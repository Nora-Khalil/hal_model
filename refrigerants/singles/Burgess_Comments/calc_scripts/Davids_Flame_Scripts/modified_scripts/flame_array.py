import subprocess
import multiprocessing
from multiprocessing.pool import Pool
import pandas as pd
import numpy as np
import time
import os
import logging
import cantera as ct
from matplotlib import pyplot as plt
from collections import defaultdict
import sys
import itertools
import scipy
import scipy.optimize

#path= '/scratch/westgroup/chlorination/git_hub/halogen-data/cantera/2-BTP/3_26_2020_gri_seed/chem_annotated.cti'
path = '/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/cantera/difftool/NIST/2-BTP_kinetics.cti'

def calculate_flame_speed(condition):
    
    def extrapolate_uncertainty(grids, speeds, plot=False):
        """
        Given a list of grid sizes and a corresponding list of flame speeds,
        extrapolate and estimate the uncertainty in the final flame speed.
        Also makes a plot, unless called with `plot=False`.
        """
        grids = list(grids)
        speeds = list(speeds)

        def speed_from_grid_size(grid_size, true_speed, error):
            """
            Given a grid size (or an array or list of grid sizes)
            return a prediction (or array of predictions)
            of the computed flame speed, based on 
            the parameters `true_speed` and `error`.
            
            It seems, from experience, that error scales roughly with
            1/grid_size, so we assume that form.
            """
            return true_speed +  error * np.array(grid_size)**-1.

        # Fit the chosen form of speed_from_grid_size, to the last four
        # speed and grid size values.
        popt, pcov = scipy.optimize.curve_fit(speed_from_grid_size, grids[-4:], speeds[-4:])

        # How bad the fit was gives you some error, `percent_error_in_true_speed`.
        perr = np.sqrt(np.diag(pcov))
        true_speed_estimate  = popt[0]
        percent_error_in_true_speed = 100.*perr[0] / popt[0]
        logging.info("Fitted true_speed is {} +/- {} cm/s ({})".format(
            popt[0]*100,
            perr[0]*100,
            percent_error_in_true_speed
            ))
        
        # How far your extrapolated infinite grid value is from your extrapolated 
        # (or interpolated) final grid value, gives you some other error, `estimated_percent_error`
        estimated_percent_error = 100. * (
                                    speed_from_grid_size(grids[-1], *popt) - true_speed_estimate
                                    ) / true_speed_estimate
        logging.info("Estimated error in final calculation {:.1f}%".format(estimated_percent_error))

        # The total estimated error is the sum of these two errors.
        total_percent_error_estimate = abs(percent_error_in_true_speed) + abs(estimated_percent_error)
        logging.info("Estimated total error {:.1f}%".format(total_percent_error_estimate))
        
        if plot:
            plt.semilogx(grids,speeds,'o-')
            plt.ylim(min( speeds[-5:]+[true_speed_estimate-perr[0]])*.95 ,
                    max( speeds[-5:]+[true_speed_estimate+perr[0]])*1.05 )
            plt.plot(grids[-4:], speeds[-4:], 'or')
            extrapolated_grids = grids + [grids[-1] * i for i in range(2,8)]
            plt.plot(extrapolated_grids,speed_from_grid_size(extrapolated_grids,*popt),':r')
            plt.xlim(*plt.xlim())
            plt.hlines(true_speed_estimate, *plt.xlim(), colors=u'r', linestyles=u'dashed')
            plt.hlines(true_speed_estimate+perr[0], *plt.xlim(), colors=u'r', linestyles=u'dashed', alpha=0.3)
            plt.hlines(true_speed_estimate-perr[0], *plt.xlim(), colors=u'r', linestyles=u'dashed', alpha=0.3)
            plt.fill_between(plt.xlim(), true_speed_estimate-perr[0],true_speed_estimate+perr[0], facecolor='red', alpha=0.1 )

            above = popt[1]/abs(popt[1]) # will be +1 if approach from above or -1 if approach from below

            plt.annotate("",
                        xy=(grids[-1], true_speed_estimate),
                        xycoords='data',
                        xytext=(grids[-1], speed_from_grid_size(grids[-1], *popt)),
                        textcoords='data',
                        arrowprops=dict(arrowstyle='|-|, widthA=0.5, widthB=0.5', linewidth=1,
                                        connectionstyle='arc3',
                                        color='black', shrinkA=0, shrinkB=0),
                        )

            plt.annotate("{:.1f}%".format(abs(estimated_percent_error)),
                        xy=(grids[-1], speed_from_grid_size(grids[-1], *popt)),
                        xycoords='data',
                        xytext=(5,15*above),
                        va='center',
                        textcoords='offset points',
                        arrowprops=dict(arrowstyle='->',
                                        connectionstyle='arc3')
                        )

            plt.annotate("",
                        xy=(grids[-1]*4, true_speed_estimate-(above*perr[0])),
                        xycoords='data',
                        xytext=(grids[-1]*4, true_speed_estimate),
                        textcoords='data',
                        arrowprops=dict(arrowstyle='|-|, widthA=0.5, widthB=0.5', linewidth=1,
                                        connectionstyle='arc3',
                                        color='black', shrinkA=0, shrinkB=0),
                        )
            plt.annotate("{:.1f}%".format(abs(percent_error_in_true_speed)),
                        xy=(grids[-1]*4, true_speed_estimate-(above*perr[0])),
                        xycoords='data',
                        xytext=(5,-15*above),
                        va='center',
                        textcoords='offset points',
                        arrowprops=dict(arrowstyle='->',
                                        connectionstyle='arc3')
                        )

            plt.ylabel("Flame speed (m/s)")
            plt.xlabel("Grid size")
            plt.show()
        
        return true_speed_estimate, total_percent_error_estimate


    def make_callback(flame):
        """
        Create and return a callback function that you will attach to 
        a flame solver. The reason we define a function to make the callback function,
        instead of just defining the callback function, is so that it can store
        a pair of lists that persist between function calls, to store the 
        values of grid size and flame speed.
        
        This factory returns the callback function, and the two lists:
        (callback, speeds, grids)
        """
        speeds = []
        grids = []

        def callback(_):
            speed = flame.u[0]
            grid = len(flame.grid)
            speeds.append(speed)
            grids.append(grid)
            logging.info("Iteration {}".format(len(grids)))
            logging.info("Current flame speed is is {:.4f} cm/s".format(speed*100.))
            if len(grids) < 5:
                return 1.0 # 
            try:
                extrapolate_uncertainty(grids, speeds)
            except Exception as e:
                logging.info("Couldn't estimate uncertainty. " + str(e))
                return 1.0 # continue anyway
            return 1.0
        return callback, speeds, grids


    sensitivity = condition['sensitivity']
    path = condition['path']
    phi = condition['phi']
    agent = condition['agent']
    suppressant_frac = condition['suppressant_frac']
    index = condition['index']
    logging.info(condition)
    logging.info(f'processing simulation for {index}')
    gas = ct.Solution(path)

    def get_suppressant_X(mol_frac,phi):
        num = (1+3.76+0.5*phi)*mol_frac
        denom = 1-mol_frac
        return num/denom

    X_suppress = get_suppressant_X(suppressant_frac,phi)
    try:
        gas.TPX = 300, ct.one_atm, {'O2(4)':1.0, 'N2':3.76, 'CH4(3)':0.5*phi, agent: X_suppress}
    except:
        gas.TPX = 300, ct.one_atm, {'O2':1.0, 'N2':3.76, 'CH4':0.5*phi, agent: X_suppress}
    
    initial_mol_fracs = gas.mole_fraction_dict()
    logging.info(gas.mole_fraction_dict())
    #width = 0.04
    width = 0.08
    flame = ct.FreeFlame(gas, width=width)
    #flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
    flame.set_refine_criteria(ratio=5, slope=0.25, curve=0.27)
    
    callback, speeds, grids = make_callback(flame)
    flame.set_steady_callback(callback)
    
    loglevel = 1
    logging.info("solving flame...")
    try:
        flame.solve(loglevel=loglevel, auto=True)
        Su0 = flame.u[0] * 100
        logging.info("Flame Speed is: {:.2f} cm/s".format(Su0))
        best_true_speed_estimate, best_total_percent_error_estimate =  extrapolate_uncertainty(grids, speeds)
    except:
        logging.info("Flame speed calculation failed :<")
        Su0 = None
        best_true_speed_estimate = best_total_percent_error_estimate = None

    flame_speeds = {}
    flame_speeds['index'] = index
    flame_speeds['agent'] = agent
    flame_speeds['agent_volume_frac'] = suppressant_frac
    flame_speeds['phi'] = phi
    flame_speeds['T(K)'] = 300
    flame_speeds['P(atm)'] = 1
    flame_speeds['mech'] = path
    flame_speeds['Su(cm/s)'] = Su0
    flame_speeds['True_Su(m/s)'] = best_true_speed_estimate
    flame_speeds['percent_error'] = best_total_percent_error_estimate
    flame_speeds['T'] = gas.T
    if 'literature' in path:
        flame_speeds['CH4'] = initial_mol_fracs.get('CH4')
        flame_speeds['O2'] = initial_mol_fracs.get('O2')
        flame_speeds['N2'] = initial_mol_fracs.get('N2')
        flame_speeds['2-BTP'] =  initial_mol_fracs.get('BTP')
        flame_speeds['CF3Br'] =  initial_mol_fracs.get('CF3BR')
    else:
        flame_speeds['CH4'] = initial_mol_fracs.get('CH4(3)')
        flame_speeds['O2'] = initial_mol_fracs.get('O2(4)')
        flame_speeds['N2'] = initial_mol_fracs.get('N2')
        flame_speeds['2-BTP'] =  initial_mol_fracs.get('2-BTP(1)')
        #flame_speeds['CF3Br'] =  initial_mol_fracs.get('CF3Br(2)')  ###################################################
    #flame_speeds.update(initial_mol_fracs)
    
    try:
        df = pd.DataFrame([flame_speeds])
        df = df[sorted(df.columns.values.tolist())]
        #out = f'./scratch/flame_speeds_rmg_{index}.csv'
        #df.to_csv(out,index=False)
        outdir = os.path.dirname(path)
        outpath = os.path.join(outdir,'flame_speeds.csv')
        if not os.path.exists(outpath):
            df.to_csv(outpath,index=False)
        else:
            df.to_csv(outpath, mode='a', header=False, index=False)
    except:
        logging.info(f'failed to save results for {index}')


    if sensitivity is True:
        sensitivities = pd.DataFrame(data=[], index=gas.reaction_equations(range(gas.n_reactions)))
        # Set the value of the perturbation
        dk = 1e-2

        # Create an empty column to store the sensitivities data
        sensitivities["baseCase"] = ""
        for m in range(gas.n_reactions):
            gas.set_multiplier(1.0) # reset all multipliers                                                                     
            gas.set_multiplier(1+dk, m) # perturb reaction m   

            # Always force loglevel=0 for this
            # Make sure the grid is not refined, otherwise it won't strictly 
            # be a small perturbation analysis
            flame.solve(loglevel=0, refine_grid=False)

            # The new flame speed
            Su = flame.u[0]

            sensitivities["baseCase"][m] = (Su-Su0)/(Su0*dk)

        # This step is essential, otherwise the mechanism will have been altered
        gas.set_multiplier(1.0)
        file_name = name + '_' +str(XO2) + 'O2_' + str(phi) + 'phi_' + str(T0) +'K_' + str(P) + 'atm'
        sensitivities.sort_values('baseCase')
        sensitivities.to_csv(os.path.join('sensitivities',"new",file_name + '.csv'))
        # Reaction mechanisms can contains thousands of elementary steps. Choose a threshold
        # to see only the top few
        # threshold = 0.03

        # firstColumn = sensitivities.columns[0]

        # # For plotting, collect only those steps that are above the threshold
        # # Otherwise, the y-axis gets crowded and illegible
        # sensitivitiesSubset = sensitivities[sensitivities[firstColumn].abs() > threshold]
        # indicesMeetingThreshold = sensitivitiesSubset[firstColumn].abs().sort_values(ascending=False).index
        # sensitivitiesSubset.loc[indicesMeetingThreshold].plot.barh(title="Sensitivities for " + file_name,
        #                                                         legend=None)
        # plt.gca().invert_yaxis()

        # plt.rcParams.update({'axes.labelsize': 20})
        # plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$');

        # # Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
        # plt.savefig(os.path.join('sensitivities','plots',file_name + '.pdf'), dpi=300)

    return True


# NIST_seed_path = os.path.join(directory,'seed-NIST/chemkin/chem.cti')
# no_seed_path = os.path.join(directory,'no-seed/chemkin/chem.cti')
# combustion_core_path = os.path.join(directory,'seed-combustion_core/chemkin/chem.cti')
# GRI_seed_path = os.path.join(directory,'GRI-MECH-seed/chemkin/chem.cti')


# O2_comps = [.21,.25,.30,.35,.40]
# phis = [0.8,0.9,1.0,1.1,1.2,1.3]
# Temps = [298.,350.,375.,400.]
# Pressures = [1.0,2.0,2.5,3.0]
# Ts_Ps = zip(Temps,Pressures)
# O2_phis = [x for x in itertools.product(O2_comps,phis)]
# conditions = [x for x in itertools.product(O2_phis,Ts_Ps)]
# conditions = [(x[0],x[1],y[0],y[1]) for x,y in conditions] #O2_comp,phi,T,P,su
# conditions = [(0.21,1.0,298.0,1.0),]

# paths = [
# ('seed-combustion_core/chemkin/chem.cti','combustion_core'),
# ('seed-NIST/chemkin/chem.cti','NIST')
# ]

# paths =[
#     ('NIST-Seed-New/cantera/chem.cti','NIST_Seed_pdep'),
# ]

# paths = [
#     ("origin.cti","origin"),
# ]
#path = '/scratch/westgroup/chlorination/git_hub/halogen-data/cantera/2-BTP/mech/reduced_chem.cti'
os.makedirs('./scratch',exist_ok=True)
array_index = int(os.getenv('SLURM_ARRAY_TASK_ID'))
print(array_index)
#path = '/scratch/westgroup/chlorination/git_hub/halogen-data/cantera/2-BTP/gri/chem_annotated.cti'
#path= '/scratch/westgroup/chlorination/git_hub/halogen-data/cantera/2-BTP/3_26_2020_gri_seed/chem_annotated.cti'
#path = '/scratch/westgroup/chlorination/git_hub/halogen-data/cantera/2-BTP/literature_mech/2-BTP_kinetics.cti'
volume_fracs = [0,0.0025,0.005,0.0075,0.01,0.015,0.02,0.025,0.03,0.035,0.040,0.045,0.05]
#phis = [0.5,0.75,1,1.25]
phis = [0.5,0.6,0.7,0.8,0.9,1,1.1,1.2]
if 'NIST' in path:
    suppressants = ['BTP'] #,'CF3BR']
else:
    suppressants = ['2-BTP(1)']#,'CF3Br(2)'] ###########################################################
conditions = []
i = 0
for frac, phi, agent in itertools.product(volume_fracs,phis,suppressants):
    if frac == 0 and agent in ('2-BTP(1)','BTP'):
        continue
    condition = {'index':i}
    condition['phi'] = phi
    condition['path'] = path
    condition['sensitivity'] = False
    condition['suppressant_frac'] = frac
    condition['agent'] = agent
    conditions.append(condition)
    i += 1
conditions.sort(key=lambda x: x['index'])
for i in conditions: 
    print(i)
condition = conditions[array_index-1]
calculate_flame_speed(condition)

# for i,v in enumerate(volume_fracs):
#     condition = {'phi':1}
#     condition['path'] = path
#     condition['sensitivity'] = False
#     condition['suppressant_frac'] = v
#     conditions.append(condition)

# pool = Pool(processes=20)
# results = pool.map(calculate_flame_speed,conditions)
# pool.close()
# pool.join()
# flame_speeds = [x for x in results]
# results_df = pd.DataFrame(flame_speeds)
# results_df.to_csv('/scratch/westgroup/chlorination/git_hub/halogen-data/cantera/2-BTP/mech/flame_speeds.csv',index=False)


# paths =[
#     ('NIST-Seed-New/cantera/chem.cti','NIST_Seed_pdep'),
# ]

#('NIST-Seed-New-wo-FAbs/cantera/chem.cti','NIST-Seed-No-FAbs')

# paths = [
#     ('NIST_2018_model.cti','NIST_2018'),
#     ('NIST_2019_model.cti','NIST_2019')
# ]

#('GRI-MECH-seed/chemkin/chem.cti','GRI_MECH'),
#('no-seed/chemkin/chem.cti','no_seed'),

# for path in paths:

# manager = Manager()
# flame_speeds = manager.list()

# processes = []
# # for i in range(len(conditions)):
# #     p = Process(target=calculate_flame_speed, 
# #     args=(path,conditions[i]))
# #     processes.append(p)

# for i in range(len(paths)):
#     p = Process(target=calculate_flame_speed, 
#     args=(paths[i],conditions[0]))
#     processes.append(p)

# active_processes = []
# for process in processes:
#     if len(active_processes) < multiprocessing.cpu_count():
#         process.start()
#         active_processes.append(process)
#         continue
#     else:
#         one_done = False
#         while not one_done:
#             for i, p in enumerate(active_processes):
#                 if not p.is_alive():
#                     one_done = True
#                     break
#         process.start()
#         active_processes[i] = process

# complete = np.zeros_like(active_processes, dtype=bool)
# while not np.all(complete):
#     for i, p in enumerate(active_processes):
#         if not p.is_alive():
#             complete[i] = True
#         time.sleep(2)

#     # time.sleep(10)

# columns = ['model','XO2','phi','T(K)','P(atm)','Su(cm/s)','True_Su(cm/s)','percent_error']
# flame_data = [x for x in flame_speeds]
# results_df = pd.DataFrame(flame_data,columns=columns)
# df_path = os.path.join('flame_speeds',path +'.csv')
# results_df.to_csv(df_path)
    

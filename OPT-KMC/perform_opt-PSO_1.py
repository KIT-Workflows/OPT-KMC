import os, sys
import yaml
import time
import shutil
import tarfile
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sp
import nevergrad as ng
from random import shuffle, random
from itertools import accumulate
from scipy import interpolate

# Extract the core code from the tar.xz file
with tarfile.open('core_code.tar.xz', 'r:xz') as f:
    f.extractall()

from core_code.support import *
from core_code.collection import *
from core_code.update import *

def kmc_code (br_change):

    barriers = [0.0227, 0.227, 0.025, 0.25, 0.2, 0.3, 0.05, 0.33, 0.343, 0.343, 0.345, 0.340]

    barriers[:8] = br_change

    for i in [10,11,12]:
        barriers[i] = br_change[8]

        if not os.path.exists('data_opt-PSO'): os.makedirs('data_opt-PSO')

    with open('data_input/args.yml') as file:
        arguments = yaml.full_load(file)

    xdim = arguments['xdim']
    ydim = arguments['ydim']
    frac_save = arguments['frac_save']
    vis_div = arguments['vis_div']
    T = arguments['T']
    concentration = br_change[9]
    # concentration = arguments['concentration']

    # set the time at 0
    t = 0
    # creating an empty list for collecting events at each step.
    pre = []
    # set the number of steps (KMC steps)
    num = 0
    # event counter and residence time
    '''to keep track of events and their consumed time, two lists will be created'''
    C = [0] * 9
    res_time = [0] * 9
    # to save output
    vis_num = int(frac_save / vis_div)
    vis_ = 0
    lists = pd.DataFrame()
    # a counter to keep track of last saved step
    keep = 0
    tracking = {}
    post_list = []

    e1 = Event(["Surf", "EC-Li+"], ["Surf", "EC-Li"], 0.24, 0) # first electron reduction wirh passivation layer rate
    e2 = Event(["EC-Li"], ["P"], 0.24, 0) # ring openning
    e3 = Event(["Surf", "P"], ["Surf", "I"], 0.376, 0)# second electron reduction  / counter for C2h4
    e4 = Event(["I", "EC-Li+"], ["F", "S"], 0.376, 0)
    e5 = Event(["P", "P"], ["O", "S"], 0.372, 0) #Li2EDC  / need a c2h4
    # Both inert and active solvents diffusion
    e6 = Event(["O", "S"], ["S", "O"], 0.3755, 0)
    e7 = Event(["P", "S"], ["S", "P"], 0.374, 0)
    e8 = Event(["EC-Li", "S"], ["S", "EC-Li"], 0.2, 0)
    e9 = Event(["EC-Li+", "S"], ["S", "EC-Li+"], 0.2, 0)
    e10 = Event(["EC-Li+", "S"], ["S", "EC-Li+"], 0.2, 0)
    e11 = Event(["EC-Li+", "S"], ["S", "EC-Li+"], 0.2, 0)
    e12 = Event(["EC-Li+", "S"], ["S", "EC-Li+"], 0.2, 0)

    events = [e1,e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12]

    # taking barriers from external list
    for i in events:
        i.barrier = barriers[events.index(i)]
        # print(i.barrier)
    # calculate rates
    for i in events:
        i.rate = bar_rate(i.barrier, T)

    lattice = Lattice(xdim, ydim)
    points = len(lattice)

    for i in lattice:
        if i.coordinate[1] == 0:
            i.species = "E"
        else:
            i.species = "S"

    all_S = [i for i in lattice if i.species == "S"]
    # shuffle the list
    shuffle(all_S)
    target = list_splitter(all_S, concentration)
    active_ec = target[0]
    #the starting number of EC-Li+
    ref_ec = len(active_ec)
    for i in all_S:
        # if i in ec_list:
        if i in active_ec:
            i.species = "EC-Li+"

    Li_con = len([i for i in lattice if i.species == "EC-Li+"])/xdim*ydim
    top_list = [i for i in lattice if i.species in ["E", "EC-Li+"]]
    EE = [i for i in lattice if i.species in ["E"]]

    start_time = time.time()
    '''
    this loop determines termination condition
    True  : up to end on no more move
    < step: at desired step
    <t    : at desired time
    '''
    while t*1e9<=9.0:
        '''starting the main part: KMC'''
        temp = []
        del temp[:]
        #get single representative of each cluster
        Li_con = Li_con_func(lattice,xdim,ydim)
        top_list = list(dict.fromkeys((top_list + EE_func(lattice))))
        temp = pre_event(top_list, Li_con, events, e1,e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12)
        pre = temp + pre
        if len(pre) == 0:
            spec_in_top = [i for i in lattice if i.species not in ["S"] and i not in top_list]
            catch = []
            spec_top = []
            for i in spec_in_top:
                if i not in catch:
                    temp = [i] + i.bonds
                    catch = catch + temp
                    spec_top.append(i)

            pre = pre_event(top_list, Li_con, events, e1,e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12)
            del spec_top[:]
            del spec_in_top[:]
            if len(pre) == 0:
                coordinate_, spec_, bonds_, status_ , color_, time_ , step_ = [], [], [], [], [], [], []
                for i in lattice:
                    coordinate_.append(i.coordinate)
                    spec_.append(i.species)
                    bonds_.append([lattice.index(ii) for ii in i.bonds])
                    status_.append(i.status[3][1])
                    color_.append(colors(i))
                    time_.append(t)
                    step_.append(num)
                object_to_write = [step_, coordinate_, spec_, bonds_, status_,color_, time_ ]
                temp_data = pd.DataFrame(object_to_write)
                lists = pd.concat([lists, temp_data.T])
                lists.reset_index(drop=True, inplace=True)
                #creating the dataframe from traj
                lists.to_pickle('data_opt-PSO/traj_' + str(keep + frac_save))
                #occ rate dataframe
                del lists
                soc_list = []
                vis_ = 0
                stat = [time.time() - start_time, t, num]
                tracking["Occurrence"] = C
                tracking["Residence_time"] = round_float_list(res_time, 30)
                tracking["Status"] = round_float_list(stat, 30)
                try:
                    with open('data_opt-PSO/tracking.yml', 'w') as out:
                        yaml.dump(tracking, out, default_flow_style=False)
                except IOError:
                    print("I/O error")
                break

        else:
            '''selection scheme'''
            R_index = []
            del R_index[:]
            # t1 = time.time()
            Events_per_type = [[i for i in pre if i.new_coord[-1] == j] for j in range(12)]
            giant_pre = [list(accumulate([i.rate for i in Events_per_type[j]])) for j in range(len(Events_per_type))]
            # print(giant_pre)
            giant_pre_index = [giant_pre.index(i) for i in giant_pre if i != []]
            giant_pre = [i for i in giant_pre if i != []]
            R_index = [0] + [i[-1] for i in giant_pre]
            R_tot = list(accumulate(R_index))
            cum_rate = R_tot[-1]
            '''Event type selection among 25 event types'''
            while True:
                rho1 = random()
                r = rho1 * (R_tot[-1])
                if (R_tot[0] < r and r <= R_tot[-1]):
                    break
            rate = []
            for i in range(1, len(R_tot)):
                temp = ratee(i, r, R_tot)
                rate.append(temp)
            which = [i for i in rate if i]
            # index
            type_index = which[0][1]
            # which one in main list
            type_index2 = giant_pre_index[type_index]
            '''the second selection same as the first one'''
            rho2 = random()
            event_at = int(len(Events_per_type[type_index2]) * rho2)
            event = Events_per_type[type_index2][event_at]
            """
            Time increment
            """
            dt = -(np.log(random()) / cum_rate)
            t = t + dt
            """
            UPDATING the list of sites using the function considering the type of selected event.
            """
            lattice = update(event, t, lattice, e1,e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12)

            """
            removing old events related to selected sites and nbrs
            """
            s1 = event.new_coord[0]
            s2 = event.new_coord[1]
            main_ = [s1] + s1.nbr + [s2] + s2.bonds + s2.nbr
            main_ = list(dict.fromkeys(main_))
            main_before = main_
            pre = [i for i in pre if i.new_coord[0] not in main_before and i.new_coord[1] not in main_before]
            counter = event.new_coord[-1]
            if counter in [11,10,9,8]:
                counter = 8
            '''prints'''
            num += 1
            # print("=============", num, "================")
            # print("selected event type:", event.reactant, "--->", event.product, event.rate)
            #  '''Keep ec-Li+ concentration constant'''
            # if event.new_coord[-1] == 0:
            #     # calculating th enumber of EC-Li+ from the concentration
            #     Li_curr = Li_con_func(lattice,xdim,ydim)*xdim*ydim
            #     #the starting number of EC-Li+
            #     if Li_curr <= ref_ec :
            #         #choice a random site in the lattice with species S and make it EC-Li+
            #         for iter in range(int(ref_ec - Li_curr)):
            #             S_new = empty_func(lattice)
            #             new_S = choice(S_new)
            #             new_S.species = "EC-Li+"
            """
            Top-list changes
            """
            del top_list[:]
            '''just kick S and A out of the  list'''
            top_list = [i for i in main_before if i.species not in ["S"]]
            res_time[counter] += dt
            C[counter] += 1
            if num % vis_div == 0:
                coordinate_, spec_, bonds_, status_ , color_, time_ , step_ = [], [], [], [], [], [], []
                for i in lattice:
                    coordinate_.append(i.coordinate)
                    spec_.append(i.species)
                    bonds_.append([lattice.index(ii) for ii in i.bonds])
                    status_.append(i.status[3][1])
                    color_.append(colors(i))
                    time_.append(t)
                    step_.append(num)
                object_to_write = [step_, coordinate_, spec_, bonds_, status_,color_, time_ ]
                temp_data = pd.DataFrame(object_to_write)
                lists = pd.concat([lists,temp_data.T])
                lists.reset_index(drop=True,inplace=True)
                vis_ += 1

            if num % frac_save == 0:
                stat = [time.time() - start_time, t, num]
                # creating the dataframe from traj
                lists.to_pickle('data_opt-PSO/traj_' + str(num))
                # occ rate dataframe
                del lists
                lists = pd.DataFrame()
                vis_ = 0
                tracking["Occurrence"] = C
                tracking["Residence_time"] = round_float_list(res_time, 30)
                tracking["Status"] = round_float_list(stat, 30)
                try:
                    with open('data_opt-PSO/tracking.yml', 'w') as out:
                        yaml.dump(tracking, out, default_flow_style=False)
                except IOError:
                    print("I/O error")
                keep = num

    dirc = [i for i in os.listdir('data_opt-PSO') if 'traj_' in i]
    it = max([int(i) for i in [s.replace("traj_", "") for s in dirc]])
    points = xdim*ydim
    df_conc = pd.DataFrame()

    df_conc = pd.DataFrame()
    for ii in range(frac_save,it+frac_save,frac_save):
        # print(ii)
        if os.path.getsize(('data_opt-PSO/traj_' + str(ii))) > 0:
            data = pd.read_pickle('data_opt-PSO/traj_' + str(ii))
            data.columns = ['step','coord','spec','bonds','status','color','time']
            EC_Lip = len(data[data['color']=='brown'])
            EC_Li = len(data[data['color']=='olive'])
            P = len(data[data['color']=='green'])
            I = len(data[data['color']=='pink'])
            F = len(data[data['color']=='red'])
            O = len(data[data['color']=='orange'])
            time_df = data.iloc[1]['time']
            fractions_ = [EC_Lip, EC_Li, P, I, F, O]
            fractions_ = [fra/points for fra in fractions_]
            fractions_.append(time_df)
            fractions = pd.DataFrame(fractions_).T
            df_conc = pd.concat([df_conc, fractions])

    df_conc.columns = [r'EC-Li$^+$', r'EC_Li', r'C$_2$H$_4$OCO$_2$Li', 'I', r'Li$_2$CO$_3$', r'Li$_2$EDC','time']
    new_order = [-1,0,1,2,3,4,5]
    df_spec = df_conc[df_conc.columns[new_order]]

    x_kmc = df_spec['time'].values.astype(float)*1e9
    x_kmc = np.concatenate(([0], x_kmc))
    name_species_kMC = [ r'EC_Li', r'C$_2$H$_4$OCO$_2$Li', r'I', r'Li$_2$EDC']

    xnew = np.linspace(x_kmc.min(),min(x_kmc.max(),9.0),200)
    y_kmc_new = np.zeros((200,len(name_species_kMC)))

    for i, species in enumerate(name_species_kMC):
        y_kmc = df_spec[f'{species}'].values.astype(float)
        y_kmc = np.concatenate(([0], y_kmc))
        f = interpolate.interp1d(x_kmc, y_kmc)
        y_kmc_new[:,i] = f(xnew)

    x_3ds = data_3DS['Time']
    name_species_3DS = ['cEC(.-)', 'oEC(.-)', 'CO3(2-)', 'EDC(2-)']
    y_3ds_new = np.zeros((200, len(name_species_3DS)))
    for i, species in enumerate(name_species_3DS):
        y_3ds = data_3DS[f'Conc_mean_{species}']
        f = interpolate.interp1d(x_3ds,y_3ds)
        y_3ds_new[:,i] = f(xnew)

    shutil.rmtree('data_opt-PSO')

    return y_kmc_new, y_3ds_new

def plot_3DS(data_3DS):
    
    plt.style.use("data_input/figure_style.mplstyle")
    fig = plt.figure(1, figsize=(4.5,3.), dpi=300)
    ax = fig.add_subplot()

    x_3ds = data_3DS['Time']
    name_species_3DS = ['cEC(.-)', 'oEC(.-)', 'CO3(2-)', 'EDC(2-)']
    color_leg = ['olive','green','pink','orange']

    for i, species in enumerate(name_species_3DS):
        y_3ds = data_3DS[f'Conc_mean_{species}']
        # f = interpolate.interp1d(x_3ds,y_3ds)
        # y_3ds_new[:,i] = f(xnew)
        ax.plot(x_3ds,y_3ds, color=f'{color_leg[i]}', label=f'{species}')

    ax.legend(loc='center', bbox_to_anchor = [0.5,1.1], ncol=3)
    plt.savefig("3DS-inputs.png")

    return "plot_3DS done"

def collect_values(data_list):
    values_dict = {}

    for data in data_list:
        for key, value in data.items():
            if key not in values_dict:
                values_dict[key] = [value]
            else:
                values_dict[key].append(value)

    return values_dict

if __name__ == "__main__":

    ''' Read the configuration file and set the 
    parameters of the algorithm and the problem to be optimized'''

    with open('rendered_wano.yml') as file:
        wano_file = yaml.full_load(file)
    

    # Input data from Simstack
    budget_var = wano_file["budget"]
    dict_br = collect_values(wano_file["Energy Barrier"])
    
    with tarfile.open('data_input.tar.xz', 'r:xz') as f:
        f.extractall()
    
    # Target data 
    data_3DS = pd.read_csv('data_input/target_data_3DS.csv')
    sys.setrecursionlimit(100000)

    #br_change = [0.0227, 0.227, 0.025, 0.25, 0.2, 0.3, 0.05, 0.33, 0.343, 0.067]
    #br_mean = np.array([0.0227, 0.227, 0.025, 0.25, 0.2, 0.3, 0.05, 0.33, 0.343, 0.2])
    br_change = dict_br['change'] 
    br_mean = np.array(dict_br['mean'])

    plot_3DS(data_3DS)

    if not os.path.exists('analysis_opt-PSO'): os.makedirs('analysis_opt-PSO')

    if wano_file["Barriers range"]:
        br_max = np.array(wano_file["br-max and br-min"]["br-max"])
        br_min = np.array(wano_file["br-max and br-min"]["br-min"])
    else: 
        br_max = br_mean + br_mean*3
        br_min = br_mean - br_mean
    print(br_max, br_min)

    def kmc_wrapper(n_randoms):

        for i in range(len(n_randoms)):
            br_change[i] = n_randoms[i]*(br_max[i]-br_min[i])+br_min[i]

        with open('analysis_opt-PSO/Energy_barrier.csv', 'a') as file_:
            str_ = [f'{br_:15.6e}' for br_ in br_change]
            file_.write(','.join(str_) + '\n')

        try:
            y_kmc_new, y_3ds_new = kmc_code(br_change)
            loss = np.sum((y_kmc_new-y_3ds_new)**2)
        except ValueError:
            loss = 99999e10

        with open('analysis_opt-PSO/Loss_function.csv', 'a') as file_:
            file_.write(f'{loss:15.6e}' + '\n')

        return loss

    instrum = ng.p.Instrumentation(ng.p.Array(shape=(len(br_mean),)).set_bounds(lower=0.0, upper=1.0))
    PSO_optimizer = ng.optimizers.PSO
    PSO_optimizer.popsize=10
    optimizer = PSO_optimizer(parametrization=instrum, budget=budget_var)

    recommendation = optimizer.minimize(kmc_wrapper)
    print(recommendation.value)
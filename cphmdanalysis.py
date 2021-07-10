"""
The Python script parses and run calculation on the lambda and log files produced by 
continuous constant pH molecular dynamics (CpHMD) simulations at independent pH conditions 
or with pH replica exchange protocol. 
For the latter case, two different classes are used to process the log files 
(:class:`log_analysis_charmm`) and (:class:`log_analysis_amber`). 
All lambda files are processed with (:class:`lambda_data`). 
Plots can be made using (:class:`plotting_lambda_data`) for lambda files and class for
log files.
"""

# Library for CpHMD Analysis 
#----------------------------
import numpy as np
from operator import add
import pandas as pd 
import re
import glob
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import cm
import math

##########################################################                
# Author: Jack Henderson                                 #
# Date: 6/14/2021                                        #
# Version: 1.2                                           #
##########################################################

##########################################################
#                                                        #
# Lambda File Processing                                 #
#                                                        #
##########################################################
#  Used for both CHARMM and AMBER lambda files           #
##########################################################

def HH_fitting(ph, hill, pka):
        return 1 / (1 + 10**(hill*(pka - ph)))

def compute_pkas(phs, ph_object, type='Macro'):
    '''
    This function is used to calculate the pKas for all titratable residues. 
    
    phs       : A list of ph values 
    ph_object : A list of 'ph_objects' generated using :class:`lambda_data`
    type      : The default 'Macro' will calculate the macro pKa for all residues the output for each
                residue will be the a list of [pKa, Hill].
                If set to 'Micro' the micro pKas for all residues with tautomeric states will 
                be calculated and the output will be a 2D list [[res#, pKa1, hill1],[res#, pKa2, hill2]].        
    '''
    data = []
    number_of_sites = 0
    if type is 'Macro':
        number_of_sites = len(ph_object[0].s)
    if type is 'Micro':
        number_of_sites = len(ph_object[0].s_x1)
    for i in range(number_of_sites):
        if type is 'Macro':
            s_values = [an_object.s[i] for an_object in ph_object]
            diffs = [abs(s_value-0.5) for s_value in s_values]
            min_diff = min(diffs)
            if min_diff < 0.4:
                initial_params = [1, phs[diffs.index(min_diff)]]
                try:
                    fit  = curve_fit(HH_fitting, phs, s_values, p0=initial_params)
                    hill = fit[0][0]
                    pka  = fit[0][1]
                    data.append([pka, hill])
                except:
                    data.append(['NaN', 'NaN'])
            else:
                data.append(['NaN', 'NaN']) 

        if type is 'Micro':
            data_sub = []
            s_values = [[an_object.s_x1[i] for an_object in ph_object],
                        [an_object.s_x2[i] for an_object in ph_object]]
            diffs    = [[abs(s_value[1]-0.5) for s_value in s_values[0]], 
                        [abs(s_value[1]-0.5) for s_value in s_values[1]]]
            min_diffs = [min(diffs[0]), min(diffs[1])]
            for n, min_diff in enumerate(min_diffs):
                if min_diff < 0.4:
                    initial_params = [1, phs[diffs[n].index(min_diff)]]
                    try:
                        fit  = curve_fit(HH_fitting, phs, [s[1] for s in s_values[n]], p0=initial_params)
                        hill = fit[0][0]
                        pka  = fit[0][1]
                        data_sub.append([s_values[n][0][0], pka, hill])
                    except:
                        data_sub.append([s_values[n][0][0], 'NaN', 'NaN'])
                else:
                    data_sub.append([s_values[n][0][0], 'NaN', 'NaN'])
            data.append(data_sub)
    return data
    
class lambda_data:
    '''
    This class will allow you to add lambda file to perform calculations on, append additional lambda 
    files (add_l_files), compute the unprotonated fraction (compute_all_s_values), and compute the 
    running average of the unprotonated fraction (compute_all_running_s).
    '''
    def __init__(self, file_path=''):
        '''
        This function is used to generate an initial ph object, where the pH objects containes all the 
        lambda values for that particular pH environment. 

        file_path : Requies a list of files lambda files arranged from lowest to highest ph environment. 

        n_itir    : Number of Titratable Sites, This includes not only the residue but the tautomers sites.
        ititr     : Indexes for all titrabtable sites.
        n_ires    : Number of Titratable Residues, current titratable residues can be Asp, Glu, His, Cys, Lys, and Arg.
        ires      : The Resid for the titratable site. 
        ires_raw  : Resids including duplicates for lambda and tautomer values.
        itauto_raw: ID for all titratable sites for Lambda and X, Cys, Lys, and Arg Type (ID: 0 Lambda); His Type (ID: 
                    1 Lambda & 2 X); Asp & Glu types (ID: 3 Lambda & 4 X).
        itauto    : ID for all titratable site, but for lambda values only.
        info_df   : Allows you to print a pandas dataframe containing all the Ititr, Ires, and Itauto values. This 
                    values will also be index according to their position in the pandas data from which is useful 
                    when want to select a specific column of lambda or x values.
        lambda_and_x_vals: These are the columns of lambda and x values from the CpHMD lambda file output. You can 
                           select a specific column of data from this list by using the index of the info_df output.
        '''
        with open(file_path, 'r') as lf:
            #### for first line of lambda file. (ititr)
            line = lf.readline()
            ititr = [int(x) for x in line.split()[2:]] 
            self.n_ititr = len(ititr) # Number of Titratable Sites 
            self.ititr = ititr        # Index for Titrabtable Sites

            #### for second line of lambda file (ires)
            line = lf.readline()
            ires_raw = [int(x) for x in line.split()[2:]]
            ires = sorted(set([int(x) for x in line.split()[2:]]))
            self.n_ires = len(ires)  # Number of Titratable Residues 
            self.ires = ires         # Unique Resids for all Titratable Residues
            self.ires_raw = ires_raw # Resids Including Duplicates for Lambda and Tautomer Values

            #### for third line of lambda file (itauto_raw)
            line   = lf.readline()
            itauto_raw = [int(x) for x in line.split()[2:]]
            itauto     = [int(x) for x in line.split()[2:] if (int(x) == 0 or int(x) == 1 or int(x) == 3)]
            self.itauto_raw = itauto_raw # ID for all titratable residues for Lambda and Chi, used for ARG LYS (ID: 0) GLU ASP (ID: 3&4) and HSP (ID: 1&2)
            self.itauto     = itauto     # ID for all titratable residues for Lambda Only

            #### for forth line of lambd file (ParK)
            line   = lf.readline()
            park   = [float(x) for x in line.split()[2:]]
            self.park = park # IDK 
            
            #### Start Making a Lambda Files DataFrame here is only the heading information  
            global info_df
            info_df = pd.DataFrame()
            info_df['Ititr']  = ititr
            info_df['Ires']   = ires_raw
            info_df['Itauto'] = itauto_raw
            self.info_table = info_df
            
            #### The rest of the lambda file 
            self.steps = []
            self.lambda_and_x_vals = [[] for x in range(self.n_ititr)] 
            for frame, line in enumerate(lf):
                self.steps.append(int(line.split()[0]))
                for j in range(self.n_ititr):
                    self.lambda_and_x_vals[j].append(float(line.split()[j+1]))

    def add_l_file(self, file_path=''):
        '''
        The add_l_file function allows you to append additional lambda files to their corresponding ph objects. 
        '''
        with open(file_path, 'r') as lf:
            for line in lf:
                if '#' not in line:
                    for j in range(len(line.split()[1:])):
                        self.lambda_and_x_vals[j].append(float(line.split()[j+1]))
    
    # Function for finding residues - Work in Progress  
    def find_residues(self, residue_list):
        '''
        I am still considering revising this.
        '''
        resids = []
        lambda_columns = []
        for n, x in enumerate(self.ires):
            for y in residue_list:
                if x == y:
                    resids.append(n)
        return resids


    def compute_all_s_values(self, cutoff=0.2, output=False, minstep = -1, maxstep = -1):
        '''
        The compute_all_s_values function allows you to calculate the unprotonated fraction for each residue for each 
        pH object. 

        cutoff: This is the lambda and x value cutoff used to denote when the proton is present or absent or 
                or on one or the other tautomer.
        output: Setting this boolean to True will provided the user with a formatted output table of the unprotonated 
                fraction information for the whole residue and if available each tautomer site as well as the pure and 
                mixed fractions for that residue.
        minstep: This is the step to start calculating the unprotonated fraction. The default value of -1 starts the 
                 calculation at 0.
        maxstep: This is the final step used for calculating the unprotonated fraction. The default value of -1 ends
                 the calculation at the final step.
        '''
        low_cut = cutoff
        up_cut  = 1 - cutoff
        self.s = []
        self.s_x1 = []
        self.s_x2 = []
        self.mixed = []
        if output:
            print('ires | itaut | S(unprot) | Pure | Mixed')
        for n, i in zip(range(len(self.lambda_and_x_vals)), self.itauto_raw):
            # Reset Time 
            if minstep == -1 or maxstep == -1:
                minstep = 0
                maxstep = len(self.lambda_and_x_vals[n])
            # ASP, GLU, and HIS like - Values are in Agreement - Two Proton Residues
            if i == 3 or i == 1:
                # Get counts for all state criterion  
                mixed_count  = 0 
                unprot_count = 0
                prot_count   = 0
                x1_count     = 0
                x2_count     = 0
                for count, x in enumerate(zip(self.lambda_and_x_vals[n], self.lambda_and_x_vals[n+1])):
                    if count >= minstep and count <= maxstep:
                        if ((x[0] < up_cut and x[0] > low_cut) or (x[1] < up_cut and x[1] > low_cut)):
                            mixed_count += 1
                        if (x[0] >= up_cut and (x[1] >= up_cut or x[1] <= low_cut)):
                            unprot_count += 1
                        if (x[0] <= low_cut and (x[1] >= up_cut or x[1] <= low_cut)):
                            prot_count += 1
                        if i == 3:
                            if x[0] <= low_cut and x[1] >= up_cut:
                                x1_count += 1
                            if x[0] <= low_cut and x[1] <= low_cut:
                                x2_count += 1
                        if i == 1: 
                            if x[0] >= up_cut and x[1] >= up_cut:
                                x1_count += 1
                            if x[0] >= up_cut and x[1] <= low_cut:
                                x2_count += 1
                # Macro and Micro State Calculations - Try and Else added not tested
                mixed = mixed_count / (maxstep - minstep)
                pure  = 1 - mixed 
                unprot_frac = unprot_count / ((maxstep - minstep) - mixed_count)
                x1_frac = 0
                x2_frac = 0
                if i == 3:
                    try: 
                        x1_frac = unprot_count / (unprot_count + x1_count)
                        x2_frac = unprot_count / (unprot_count + x2_count)
                    except:
                        x1_frac = 0
                        x2_frac = 0
                if i == 1:
                    try:
                        x1_frac  = x1_count / (x1_count + prot_count)
                        x2_frac  = x2_count / (x2_count + prot_count)
                    except:
                        x1_frac = 0
                        x2_frac = 0    
                # Results for ASP, GLU and HSP like. 
                self.s.append(unprot_frac)
                self.s_x1.append([self.ires_raw[n], x1_frac])
                self.s_x2.append([self.ires_raw[n], x2_frac])
                self.mixed.append(mixed)
                if output:
                   print('{0:4.0f}     {1:0.0f}       {2:2.2f}      {3:2.2f}   {4:2.2f}'.format(self.ires_raw[n], self.itauto_raw[n], unprot_frac, pure, mixed))
                   print('{0:4.0f}     {1:0.0f}       {2:2.2f}'.format(self.ires_raw[n], self.itauto_raw[n+1], x1_frac))
                   print('{0:4.0f}     {1:0.0f}       {2:2.2f}'.format(self.ires_raw[n], self.itauto_raw[n+1], x2_frac))

            # LYS, CYS, and ARG like - Values are in Agreement - One Proton Residues
            if i == 0:
                # Get counts for all state criterion 
                mixed_count  = 0
                unprot_count = 0 
                for count, x in enumerate(zip(self.lambda_and_x_vals[n])):
                    if count >= minstep and count <= maxstep:
                        if (x[0] < up_cut and x[0] > low_cut):
                            mixed_count += 1
                        if x[0] >= up_cut:
                            unprot_count += 1

                mixed = mixed_count / (maxstep - minstep)
                pure  = 1 - mixed
                unprot_frac = unprot_count / ((maxstep - minstep) - mixed_count)

                # Results
                self.s.append(unprot_frac)
                self.mixed.append(mixed)
                if output:
                    print('{0:4.0f}     {1:0.0f}       {2:2.2f}      {3:2.2f}   {4:2.2f}'.format(self.ires_raw[n], self.itauto_raw[n], unprot_frac, pure, mixed))


    def compute_all_running_s(self, cutoff=0.2, minstep = -1, maxstep = -1):
        '''
        The compute_all_running_s function allows you to calculate the running average of the unprotonated 
        fraction over the simulation.

        cutoff: This is the lambda and x value cutoff used to denote when the proton is present or absent or 
                or on one or the other tautomer.
        minstep: This is the step to start calculating the unprotonated fraction. The default value of -1 starts the 
                 calculation at 0.
        maxstep: This is the final step used for calculating the unprotonated fraction. The default value of -1 ends
                 the calculation at the final step.
        '''
        low_cut = cutoff
        up_cut  = 1 - cutoff
        self.running_s = [[] for x in self.ires]
        residue_count = 0
        for n, i in zip(range(len(self.lambda_and_x_vals)), self.itauto_raw):
            # Reset Time 
            if minstep == -1 or maxstep == -1:
                minstep = 0
                maxstep = len(self.lambda_and_x_vals[n])
            # ASP, GLU, and HIS like - Values are in Agreement - Two Proton Residues
            if i == 3 or i == 1:
                # Get counts for all state criterion  
                mixed_count  = 0 
                unprot_count = 0
                for count, x in enumerate(zip(self.lambda_and_x_vals[n], self.lambda_and_x_vals[n+1])):
                    if count >= minstep and count <= maxstep:
                        if ((x[0] < up_cut and x[0] > low_cut) or (x[1] < up_cut and x[1] > low_cut)):
                            mixed_count += 1
                        if (x[0] >= up_cut and (x[1] >= up_cut or x[1] <= low_cut)):
                            unprot_count += 1
                        # Calculate Running S
                        running_s_calc = 0
                        if (count + 1) <= mixed_count: 
                            #running_s_calc = 'NaN'
                            running_s_calc = 0
                        else:
                           running_s_calc = unprot_count / ((count + 1) - mixed_count)
                        self.running_s[residue_count].append(running_s_calc)
                residue_count += 1                            
    
            # LYS, CYS, and ARG like - Values are in Agreement - One Proton Residues
            if i == 0:
                # Get counts for all state criterion 
                mixed_count  = 0
                unprot_count = 0 
                for count, x in enumerate(zip(self.lambda_and_x_vals[n])):
                    if count >= minstep and count <=maxstep:
                        if (x[0] < up_cut and x[0] > low_cut):
                            mixed_count += 1
                        if x[0] >= up_cut:
                            unprot_count += 1
                        # Calculate running s
                        running_s_calc = 0
                        if (count + 1) <= mixed_count:
                            #running_s_calc = 'NaN'
                            running_s_calc = 0
                        else:
                            running_s_calc = unprot_count / ((count + 1) - mixed_count)
                        self.running_s[residue_count].append(running_s_calc)
                residue_count += 1

class plotting_lambda_data:
    '''
    This class will perform all of the standard plotting for CpHMD lambda files allowing you to 
    plot the running S values and titration curves for individual and multiple residue selections
    '''
    def plot_running_s(self, phs, resids, titles, xlabel='Steps', steps_to_time_conversion=1, save_fig=False, output='s_conv.png'):
        '''
        In order to run this function you must have first run the 'compute_all_running_s' function from the 
        'lambda_data' class.

        self: The list of pH objects
        phs: list of the pH values that corresponding to the pH objects.
        resids: The residue value determined by the 'find_residues' function to denote which running S values 
                should be plotted.
        titles: A list of titles for the plots, this allows you to decide whether you want the title to be 
                something like 'D4', 'Asp4', or 'Asp 4' ect... 
        xlabel: The label used for the X-axis. The default is set to 'Steps'
        steps_to_time_conversion: A value to multiply the X values by to convert them from steps to time. The
                                  default value is set to 1, which would results in steps.
        save_fig: A boolean when set to True allows you to save a png file. The default value is set to False.
        output: The output name of png to save. 'save_fig' must be set to True.
        '''
        colormap = cm.get_cmap('rainbow_r', len(phs))
        colors = colormap(np.linspace(0, 1, len(phs)))
        plt.close('all')

        figcols = 4
        figrows = math.ceil(len(resids)/figcols)
        fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(14,10/figcols*figrows))
        axes = axes.flatten()

        title = titles

        for index in range(figcols*figrows):
            if index < len(resids):
                for n, ph in enumerate(phs):
                    # Plot Data
                    axes[index].plot([x*steps_to_time_conversion for x in range(len(self[n].running_s[resids[index]]))], 
                                     [y for y in self[n].running_s[resids[index]]], 
                                      color=colors[n])
                    # Plot Details
                    axes[index].set_title(title[index])                  
                    axes[index].set_ylabel('S',fontsize=15)
                    axes[index].set_xlabel('{}'.format(xlabel), fontsize=15)
                    axes[index].set_ylim(0, 1)
            else: 
                axes[index].remove()
        
        legend_poses = [[0.35, 0.95], [0.56, 0.95], [0.81, 0.95], [1.06, 0.95]]
        if len(resids) == 1:
            pose = legend_poses[0]
        if len(resids) == 2:
            pose = legend_poses[1]
        if len(resids) == 3:
            pose = legend_poses[2]
        if len(resids) >= 4:
            pose = legend_poses[3]

        fig.legend(labels=phs,  # The labels for each line
           bbox_to_anchor=pose,
           borderaxespad=0.05,  # Small spacing around legend box
           title="pH Values"    # Title for the legend
           )
           
        plt.tight_layout()
        if save_fig == True:
            plt.savefig('{}'.format(output), dpi=300, transparent=False, bbox_inches='tight')

    # Note: There seems to be some strange issue where when a residue doesn't titrate the 
    # x-axis labels don't print. Others online have seen this bug too, but there doesn't 
    # appear to be a clear fix. 
    def plot_titration_curves(self, phs, resids, titles, xrange=[0, 14], save_fig=False, output='titration_curves.png'):
        '''
        In order to run this function you must have first run the 'compute_all_s_values' function from the 
        'lambda_data' class.

        self: The list of pH objects
        phs: list of the pH values that corresponding to the pH objects.
        resids: The residue value determined by the 'find_residues' function to denote which running S values 
                should be plotted.
        titles: A list of titles for the plots, this allows you to decide whether you want the title to be 
                something like 'D4', 'Asp4', or 'Asp 4' ect...
        xrange: The range to use for the X axis and is in pH units.
        save_fig: A boolean when set to True allows you to save a png file. The default value is set to False.
        output: The output name of png to save. 'save_fig' must be set to True.
        '''
        plt.close('all')
        # Get the S_Values
        s_values = []
        for res in resids:
            s_res_sub = []
            for n, ph in enumerate(phs):
                s_res_sub.append(self[n].s[res])
            s_values.append(s_res_sub)

        # Calculate the pKas
        pkas = compute_pkas(phs, self)
        # Start Plotting
        figcols = 4
        figrows = math.ceil(len(resids)/figcols)
        fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(14,10/figcols*figrows))
        axes = axes.flatten()

        title = titles

        for index in range(figcols*figrows):
            if index < len(resids):
                # Plot Data
                axes[index].plot(phs, s_values[index], 'o', color='Black')
                if pkas[resids[index]][0] is not 'NaN':
                    # Plot pKa line 
                    axes[index].axvline(pkas[resids[index]][0], color='red')
                    # Plot fit 
                    if abs(min(s_values[index])-max(s_values[index])) > 0.1:
                        local_xs = np.linspace(xrange[0], xrange[1], 100)
                        axes[index].plot(local_xs, [HH_fitting(pkas[resids[index]][0], -1*pkas[resids[index]][1], x) for x in local_xs], color='Black')
                else:
                    axes[index].hlines(np.mean(s_values[index]), 0, 14, color='black')
                # Plot Details 
                axes[index].set_title(title[index])
                axes[index].set_xlim(xrange[0], xrange[1])
                axes[index].set_xlabel('pH', fontsize=15)
                #axes[index].set_xticks([x for x in range(xrange[0], xrange[1]+1, 2)])
                axes[index].set_ylim(0, 1)
                axes[index].set_ylabel('Unprot. Fraction', fontsize=12)
            else:
                axes[index].remove()
        
        plt.tight_layout()
        if save_fig == True:
            plt.savefig('{}'.format(output), dpi=300, transparent=False)

##########################################################
#                                                        #
# Log File Processing                                    #
#                                                        #
##########################################################
#  Used for CHARMM CpHMD, CHARMM version >C38            #
#  and the Amber Async Replica Exchange                  #
##########################################################
           
class log_analysis_charmm:
    '''
    This class will allow you to process the CHARMM (version >c38) log files and calculate the exchange 
    frequencies between pH replicas and plot each replica walk.
    '''
    def __init__(self, files_path=''):
        '''
        file_path: This needs to be a 2D list of log files arranged by simulation stages/runs then 
                   by replica with the replicas in ascending order. 
        exchange_frq: can be used to print the exchange frequencies between the pH replicas.
        '''
        # Step 1: Using the multidimensional array store the data.
        storage = []

        # Gather File Info -----
        cpus_used = [] # This is the list of the #cpus used. 
        step_size = 0
        # ----------------------
        for rep in range(len(files_path[0])):
            storage_sub = []
            for stage in range(len(files_path)):
                with open(files_path[stage][rep], 'r') as f:
                    flag = 0
                    for line in f:
                        if 'PH-REX>' in line:
                            flag += 1
                            statement = str(line[9:16]).strip()
                            neighbor  = int(line[36:40])
                            step      = int(line[44:])

                            if rep % len(files_path[0]) == 0 and flag == 1: # store information
                                step_size = step
                                cpus      = neighbor
                                cpus_used.append(cpus) # This keeps a memory of the number of #cpu changes
                            if rep == 0 or rep == len(files_path[0])-1:
                                storage_sub.append([step, statement, neighbor/cpus_used[stage]])
                                storage_sub.append([step+step_size, 'N/A', 'N/A'])
                            else:
                                storage_sub.append([step, statement, neighbor/cpus_used[stage]])                            
            storage.append(storage_sub)
        # Step 2: Construct the RepWalk and Calculate the Exchange Rates
        self.repwalk = []
        repwalk_list = [rep for rep in range(len(storage))]
        self.repwalk.append([val for val in repwalk_list])
        exchange_frq_data = [0 for x in range(len(storage)-1)] 

        for step in range(len(storage[1])):
            if (step+1) % 2 != 0: # Down Swap for Odd Steps, if odd
                for i in [x for x in range(0, len(storage))]: 
                    y = storage[i][step]
                    if i > y[2]:
                        if'ACCEPT' in y[1]:
                            pose1 = 0
                            pose2 = 0
                            for pose, ph in enumerate(repwalk_list):
                                if ph == i:
                                    pose1 = pose
                                if ph == y[2]:
                                    pose2 = pose
                                    exchange_frq_data[int(y[2])] += 1
                            repwalk_list[pose1], repwalk_list[pose2] = repwalk_list[pose2], repwalk_list[pose1]
            else: # Up Swap for Even Steps, else even 
                for i in [x for x in range(1, len(storage)-1)]:
                    y = storage[i][step]
                    if i < y[2]:
                        if 'ACCEPT' in y[1]:
                            pose1 = 0
                            pose2 = 0
                            for pose, ph in enumerate(repwalk_list): # locate the position of pH for swap
                                if ph == i:
                                    pose1 = pose
                                    exchange_frq_data[i] += 1
                                if ph == y[2]:
                                    pose2 = pose
                            repwalk_list[pose1], repwalk_list[pose2] = repwalk_list[pose2], repwalk_list[pose1]
            self.repwalk.append([val for val in repwalk_list])

        self.exchange_frq = []
        for n, x in enumerate(exchange_frq_data):
            self.exchange_frq.append([n, n+1, x*2/len(storage[0])])

    def plot_replica_walk(self, xlabel='Step', steps_to_time_conversion=1, save_fig=False, output='repwalk.png'):
        '''
        Once the replica exchange files have been processed this function can be used to plot the replica 
        walks.

        xlabel: The label used for the X-axis. The default is set to 'Steps'
        steps_to_time_conversion: A value to multiply the X values by to convert them from steps to time. The
                                  default value is set to 1, which would results in steps.
        save_fig: A boolean when set to True allows you to save a png file. The default value is set to False.
        output: The output name of png to save, default name is 'repwalk.png'. 'save_fig' must be set to True.
        '''
        plt.close('all')
        figcols = 4
        figrows = math.ceil(len(self.repwalk[0])/figcols)
        fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(14,10/figcols*figrows))
        axes = axes.flatten()

        for index in range(figcols*figrows):
            if index < len(self.repwalk[0]):
                # Plot Data
                axes[index].plot([x*steps_to_time_conversion for x in range(len(self.repwalk))], [y[index] for y in self.repwalk], '-', color='Black')
                # Plot Details 
                axes[index].set_title('Rep: {}'.format(index))
                axes[index].set_ylabel('pH Replica', fontsize=15)
                axes[index].set_ylim(0, len(self.repwalk[0])-1)
                axes[index].set_xlabel('{}'.format(xlabel),fontsize=15)
            else: 
                axes[index].remove()
            plt.tight_layout()
            if save_fig == True:
                plt.savefig('{}'.format(output), dpi=300, transparent=False)

class log_analysis_amber:
    '''
    This class will allow you to process the AMBER Async CpHMD log files and calculate the exchange 
    frequencies between pH replicas and plot each replica walk.
    '''
    def __init__(self, files_path=''):
        '''
        file_path: This needs to be a 1D list of the log files for each stage of the Async CpHMD replica 
                   exchange. 
        exchange_frq: can be used to print the exchange frequencies between the pH replicas.
        '''
        self.rwalk_data = []
        all_exchange_flags = []
        prev_line  = []
        exch_flag = 0 # this flag is for testing
        for stage, file in enumerate(files_path):
            with open(file, 'r') as f:
                for count, line in enumerate(f):
                    split_line = re.split('\s+', line)
                    split_line = split_line[2:-1]

                    # Make the decoder
                    key_line = []
                    if stage > 0 and count == 0:
                        key_line = prev_line
                        # This is the heart of the matter by keeping track of the positions 
                        # of the replicas. 
                        decode = {split_line[x]:key_line[x] for x in range(len(split_line))}
                        #print('HERE {}'.format(key_line))
                        #print('HERE {}'.format(split_line))

                    save = prev_line # save the previous line so you can examine if a swap occured. 

                    # Perform the Decoding - Decoding occurs for everything after the first stage.
                    the_line = [] # This is what is actually used for the first stage it's the
                                  # same as split_line but after that it is decoded.
                    if stage == 0:
                        prev_line = split_line
                        the_line = split_line 
                    else:
                        for entry in split_line:
                            the_line.append(decode[entry])
                    
                    #print(f'cur: {count} {the_line}')
                    #print(f'sav: {count} {save}')
                    
                    # Calculate the Exchange Frequency
                    exch_flag += 1 
                    if exch_flag > 1: # The Flag here makes sure the Exch. Freq. is started on second step.
                        exchange_flags = []
                        for incr in range(len(the_line)-1):
                            a = the_line[incr:incr+2]
                            b = save[incr:incr+2]
                            if a == b[::-1]:
                                exchange_flags.append(1)
                            else:
                                exchange_flags.append(0)
                        all_exchange_flags.append(exchange_flags)
                        #print(f'Exc: {count} {exchange_flags}')

                            
                    # For decoder to remember the previously decoded lines 
                    if stage > 0:
                        prev_line = the_line
                    
                    # The Following finds where the replica is and records the walk.
                    position_collection = []
                    for n in range(len(split_line)):
                        for position, val in enumerate(the_line):
                            if n == int(val):
                                position_collection.append(position)
                    self.rwalk_data.append(position_collection)

        final_exchange_count = [0 for x in range(len(all_exchange_flags[0]))]
        for exchange_line in all_exchange_flags:
            for n_exchange, exchange_flag in enumerate(exchange_line):
                final_exchange_count[n_exchange] += exchange_flag
        #print(f'Final Exchange Count: {final_exchange_count}')
        #print(f'All Exchange Flags: {len(all_exchange_flags)}')

        self.exchange_frq = []
        for rep, val in enumerate(final_exchange_count):
            #print(f'{val} {len(all_exchange_flags)} {val/len(all_exchange_flags)}')
            self.exchange_frq.append([rep, rep+1, val/len(all_exchange_flags)])

    def plot_replica_walk(self, xlabel='Step', steps_to_time_conversion=1, save_fig=False, output='repwalk.png'):
        '''
        Once the replica exchange files have been processed this function can be used to plot the replica 
        walks.

        xlabel: The label used for the X-axis. The default is set to 'Steps'
        steps_to_time_conversion: A value to multiply the X values by to convert them from steps to time. The
                                  default value is set to 1, which would results in steps.
        save_fig: A boolean when set to True allows you to save a png file. The default value is set to False.
        output: The output name of png to save. 'save_fig' must be set to True.
        '''
        plt.close('all')
        figcols = 4
        figrows = math.ceil(len(self.rwalk_data[0])/figcols)
        fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(14,10/figcols*figrows))
        axes = axes.flatten()

        for index in range(figcols*figrows):
            if index < len(self.rwalk_data[0]):
                # Plot Data
                axes[index].plot([x*steps_to_time_conversion for x in range(len(self.rwalk_data))], [y[index] for y in self.rwalk_data], '-', color='Black')
                # Plot Details 
                axes[index].set_title('Rep: {}'.format(index))
                axes[index].set_ylabel('pH Replica', fontsize=15)
                axes[index].set_ylim(0, len(self.rwalk_data[0])-1)
                axes[index].set_xlabel('{}'.format(xlabel),fontsize=15)
            else: 
                axes[index].remove()
            plt.tight_layout()
            if save_fig == True:
                plt.savefig('{}'.format(output), dpi=300, transparent=False)

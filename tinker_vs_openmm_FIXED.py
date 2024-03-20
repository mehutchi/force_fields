#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 15:27:22 2022

@author: mehutchi
"""
import os
import matplotlib.pyplot as plt
import pickle as pickle
import numpy as np
import matplotlib
import matplotlib.colors as colors

# for openmm portion
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.unit as u
from forcebalance.nifty import printcool_dictionary
from forcebalance.molecule import Molecule
from collections import OrderedDict
import IPython

def list_splitter(raw_list):
    ''' function for splitting raw line-lists into lines of separate elements
    '''
    split_list = []
    # split a list (by spaces) into separate elements
    for i in range(len(raw_list)):
        split_list.append(raw_list[i].split())
    return split_list
# code to compute the energy of an MD trajectory frame-by-frame in OpenMM and AMOEBA and compare them

def read_coords(file, ext):    
    # splits either an xyz OR txyz MD trajectory into individual frame files
    with open(file) as raw:
        coords = raw.readlines()
    if ext == 'xyz':
        frame_start = coords[0]
        num_atoms_i = int(frame_start)
        space = 2
    elif ext == 'txyz':
        frame_start_init = coords[0].split(' ')[0:2] # get first two elements
        num_atoms_i = int(coords[0].split(' ')[0])
        space = 1
    frame_no = 0
    for i in range(len(coords)):
        if ext == 'txyz':
            frame_start = "%s %s %i\n"%(frame_start_init[0], frame_start_init[1], frame_no)
        check = coords[i]
        if coords[i] == frame_start:
            with open('frame_%i.%s'%(frame_no, ext), 'w') as out:
                out.write(coords[i])
#            if len(coords[i]) == 0:
                for j in range(1, num_atoms_i+space):
                    out.write(coords[i+j])
            frame_no += 1
    return frame_no

def extract_two_atom_dist_info(traj_info, dist_matrix):
    """Takes in a split bond order trajectory and the two atoms of interest as input and extracts BO information
    from each frame, appending a zero when the atom pair is not listed. Returns the bond orders and the element names.
    
    Example split bond order trajectory:
    
          185
    frame_0000 bond order list
        0     0    4.628  Fe  Fe
        0     1    0.400  Fe  Fe
        0     2    0.549  Fe  Fe
        0     3    0.496  Fe   C
        0     4    0.831  Fe   C
        0     5    0.842  Fe   C
        0     6    0.073  Fe   O
        0     7    0.110  Fe   O
        0     8    0.101  Fe   O
        0     9    0.091  Fe   C
        .
        .
        .
        
    Parameters
    ----------
    traj_info : np.array()
        Split list like the above example.
    BO_matrix : np.array()
        L x n x n matrix of zeros, where n is the number of atoms and L is the trajectory length.

    Returns
    -------
    BO_matrix : np.array()
        Filled in L x n x n matrix.
    """
    num_atoms = int(traj_info[0][0])
    displacement = 2
    # iterate over every entry in the list that has a length of one (each entry that contains the number of lines 
    # for that particular frame)
    frame = 0 #assuming that the frame count starts at zero
    for i in range(len(traj_info)):
        if len(traj_info[i]) == 1:
            # check the next j entries, where j is the number of bond orders in this frame, to determine whether
            # atom_one and atom_two are included in the BOs of that frame
            for j in range(num_atoms):
#                print(traj_info[i+j+displacement])
                a1 = [float(traj_info[i+j+displacement][1]), float(traj_info[i+j+displacement][2]), float(traj_info[i+j+displacement][3])]
#                for x in range(j+1, num_atoms):
#                    print(x)
#                print(range(j+1, num_atoms))
                for k in range(j+1, num_atoms):
                    a2 = [float(traj_info[i+k+displacement][1]), float(traj_info[i+k+displacement][2]), float(traj_info[i+k+displacement][3])]
#                    print("a1", a1, "a2", a2)
                    distance = dist_form(a1, a2)
                    dist_matrix[frame][j][k] = distance
#                if int(traj_info[i+j+displacement][0]) != int(traj_info[i+j+displacement][1]):
#                dist_matrix[frame][j][k] = distance                
#            print("Processed the %d data"%frame)
            frame += 1
    return dist_matrix
                
def get_dist_matrix(traj_length, num_atoms, traj_info, pickle_name):
    ''' create a list of all bond order matrices and save that into a pickle file
    if it does not already exist
    
    Parameters
    ----------
    traj_length : int
        Trajectory length.
    num_atoms : int
        Number of atoms.
    bond_order_file_name : string
        Bond order file name.
    pickle_name : string
        File name to create or open.

    Returns
    -------
    finished_BO_matrix : np.array()
        Filled in L x n x n matrix, where n is the number of atoms and L is the trajectory length.
    '''   
    """
    # if the pickle data for the BO_matrix already exists, then load it
    if os.path.exists(pickle_name):
        finished_dist_matrix = pickle.load(open(pickle_name, 'rb'))
    # if the pickle data for the BO_matrix does not already exist, create it
    else:
    """
    # create a matrix of zeros in the shape of the trajectory of dist matrices
    dist_matrix = np.zeros((traj_length, num_atoms, num_atoms))
    # create dist matrix
    finished_dist_matrix = extract_two_atom_dist_info(traj_info, dist_matrix)
    # save the numpy array using the pickle method
    pickle.dump(finished_dist_matrix, open(pickle_name, 'wb'))
    return finished_dist_matrix

def read_traj(traj_name):
    print("loading %s"%traj_name)
    with open(traj_name) as raw:
        lines = raw.readlines()
    split_lines = list_splitter(lines)
    num_atoms = int(split_lines[0][0])
    print("num_atoms: %i"%num_atoms)
    traj_length = int(len(split_lines)/(num_atoms+2)) #accounting for the atom lines +2 lines for atom count and description lines
    print("traj_length: %i"%traj_length)
    elements = []
    atomlist = []
    displacement = 2
    for i in range(num_atoms):
        atomlist.append(split_lines[i+displacement][0]+"%i"%i)
        elements.append(split_lines[i+displacement][0])
    return split_lines, num_atoms, traj_length, atomlist, elements

def dist_form(coords1, coords2):
    x1 = coords1[0]
    y1 = coords1[1]
    z1 = coords1[2]
    x2 = coords2[0]
    y2 = coords2[1]
    z2 = coords2[2]
    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def print_list(list_file):
    '''Convert a list to a string format, so lists can be included in the output file
    
    Parameters
    ----------
    list_file : list
        List to be converted into a string.
    Returns
    -------
    list_str : string
        Input list as a string (separated by commas).
    '''
    list_len = len(list_file)
    list_str = '[ '
    for i in range(list_len):
        list_str += str(list_file[i])
        # only add a comma and a space if not the last entry
        if i < (len(list_file) - 1):
            list_str += ', '
        if i == (len(list_file) - 1):
            list_str += ' ]'
    return list_str

def truncate_colormap(cmap, minval=0.0, maxval=1.0, resolution=100, n_separation=100):
    '''
    Parameters
    ----------
    cmap : cmap
        Input colormap to be adjusted.
    minval : float, default=0.0
        Minimum cutoff for colormap.
    maxval : float, default=1.0
        Maximum cutoff for colormap.
    resolution : int, default=100
        Resolution.
    n_separation : int, default=100
        Separtation.
        
    Returns
    -------
    new_cmap : cmap
        New, adjusted colormap.
    '''
    # Function to truncate a color map
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, resolution)), N=n_separation)
    return new_cmap

def plot_frame(dist_matrix, atomlist, frame, meas, tag='', title=''):
    if meas == "dist":
        title_label = 'Distance'
        dist_label = r'$\AA$ngstrom'
        jet = plt.get_cmap('jet_r')
    elif meas == "bond":
        title_label = 'Bond Order'
        dist_label = 'Bond Order'
    jet = plt.get_cmap('jet')
    print('len of shape of dist matrix = %i'%len(np.shape(dist_matrix)))
    if len(np.shape(dist_matrix)) > 2:
        num_atoms = len(dist_matrix[0][0])
        matrix = dist_matrix[frame]
    else:
        num_atoms = len(dist_matrix[0])
        matrix = dist_matrix
#    frame = 0
    fig, ax = plt.subplots()

    new_jet = truncate_colormap(jet, 0.1, 0.9)
    cmap = matplotlib.cm.get_cmap(new_jet)
    # make it so values under 0.5 are white
    cmap.set_under('white')
    
    if title == '':
        plt.title("%s Matrix for Frame: %i"%(title_label, frame))
    else:
        plt.title(title)
    im = plt.imshow(matrix, cmap = cmap, interpolation = 'nearest', aspect='auto')
    color_min = 10**-10
    plt.clim(color_min, matrix.max())
#    plt.imshow(dist_matrix[frame])
    plt.yticks(range(num_atoms), labels=atomlist)
    plt.xticks(range(num_atoms), labels=atomlist)
    ax.xaxis.tick_top()
    cbar = plt.colorbar()
    cbar.set_label(dist_label)
    plt.savefig('%s_matrix_frame_%.4i%s.jpg'%(meas, frame, tag))

def energy_diff_plot(t_list, o_list, name):
    if name == 'totalE':
        title = 'Total Energy'
    elif name =='bond':
        title = 'Bond Energy'
    elif name == 'angle':
        title = 'Angle Energy'
    elif name == 'vdw':
        title = 'VDW Energy'
    elif name == 'multi':
        title = 'Multipole Energy'
    elif name == 'strbe':
        title = 'Stretch-Bend Energy'
    elif name == 'tors':
        title = 'Torsion Energy'
    elif name == 'oopb':
        title = 'Out-of-Plane-Bend Energy'


    # AMOEBA - OpenMM energy difference plot (kJ/mol)
    diff_list = []
    for i in range(len(t_list)):
        diff_list.append(t_list[i] - o_list[i])
#        print('tinker: %.2f, openmm: %.2f'%(t_list[i], o_list[i]))
    plt.figure()
    plt.plot(range(len(diff_list)), diff_list)
    plt.title('Tinker - OpenMM %s Difference vs Frame Number'%title)
    plt.xlabel('MD frame')
    plt.ylabel('TINKER - OpenMM energy difference (kJ/mol)')
    plt.savefig('diff_vs_frame_%s.pdf'%name)
    plt.savefig('diff_vs_frame_%s.png'%name, dpi=300)

    plt.close()
    
    return diff_list


def N_comp(n, diff_list):
    # find N largest difference
#        n = 100
    diff_sort = sorted(diff_list)
    print(diff_sort[-n:])
    
    traj, num_atoms, traj_length, atomlist, elements = read_traj('all.xyz')
    dist_matrix = get_dist_matrix(traj_length, num_atoms, traj, 'all_dist_matrices.p')
    
    plot_frame(dist_matrix, atomlist, 0, 'dist')

    indices_top = np.argpartition(diff_list, -n)[-n:]
    print(indices_top)
    
    # create a collection of matrices with the largest energy differences
    top = []
    for i in indices_top:
        top.append(dist_matrix[i])
    top_m = np.array(top)
    ave = top_m.mean(axis=0)
    test2 = 0
    
    top_wo_frame0 = []
    for i in indices_top:
        top_wo_frame0.append(np.subtract(dist_matrix[i], dist_matrix[0]))
    top_adj = np.array(top_wo_frame0)
    ave_top = top_adj.mean(axis=0)
    
    plot_frame(ave_top, atomlist, 0, 'dist', '_adj_top_100', 'Ave Dist Matrix of top 100 Tink-OpMM after subtracting frame0')
    
    indices_bot = np.argpartition(diff_list, n)[:n]
    print(indices_bot)

    bot_wo_frame0 = []
    for i in indices_bot:
        bot_wo_frame0.append(np.subtract(dist_matrix[i], dist_matrix[0]))
    bot_adj = np.array(bot_wo_frame0)
    ave_bot = bot_adj.mean(axis=0)

    plot_frame(ave_bot, atomlist, 0, 'dist', '_adj_bot_100', 'Ave Dist Matrix of bottom 100 Tink-OpMM after subtracting frame0')
    
    with open('all_output.txt', 'w') as out:
        out.write('top 100 frames, AMOEBA-OpenMM\n')
        out.write(print_list(indices_top))
        out.write('\ntop 100 frames, AMOEBA-OpenMM\n')
        out.write(print_list(indices_bot))


                
def main():
    
    ### following section commented out because files aready obtained
    '''
    # read xyz coordinates
#    num_frames_xyz = read_coords('first_three.xyz', 'xyz')
    num_frames_xyz = read_coords('all.xyz', 'xyz')

    # read tinker coordinates
#    num_frames_txyz = read_coords('first_three.txyz', 'txyz')
    num_frames_txyz = read_coords('all.arc', 'txyz')
    
    num_frames_xyz = num_frames_txyz = 10001
    
    if num_frames_xyz != num_frames_txyz:
        raise ValueError('number of frames in xyz and txyz do not match')
    for i in range(num_frames_xyz):
        with open('inputA_%i.txt'%i, 'w') as ifile:
            ifile.write('frame_%i.txyz\n'%i)
            ifile.write('urea_final.prm\n')
            ifile.write('E\n')
        os.system("analyze < inputA_%i.txt > outputA_%i.txt"%(i, i))
    '''
    num_frames_xyz = num_frames_txyz = 10001

    pickle_tot_E_t = 'tinker_tot_E_dict.p'
    pickle_bond_t = 'tinker_bond_dict.p'
    pickle_angle_t = 'tinker_angle_dict.p'
    pickle_vdw_t = 'tinker_vdw_dict.p'
    pickle_multi_t = 'tinker_multi_dict.p'
    pickle_strbe_t = 'tinker_strbe_dict.p'
    pickle_tors_t = 'tinker_tors_dict.p'
    pickle_oopb_t = 'tinker_oopb_dict.p'


    if os.path.exists(pickle_tot_E_t) and os.path.exists(pickle_bond_t) and os.path.exists(pickle_angle_t) and os.path.exists(pickle_vdw_t) and os.path.exists(pickle_multi_t) and os.path.exists(pickle_strbe_t) and os.path.exists(pickle_tors_t) and os.path.exists(pickle_oopb_t):
        t_tot_E_dict = pickle.load(open(pickle_tot_E_t, 'rb'))
        t_tot_E_list = []
        t_bond_dict = pickle.load(open(pickle_bond_t, 'rb'))
        t_bond_list = []
        t_angle_dict = pickle.load(open(pickle_angle_t, 'rb'))
        t_angle_list = []
        t_vdw_dict = pickle.load(open(pickle_vdw_t, 'rb'))
        t_vdw_list = []
        t_multi_dict = pickle.load(open(pickle_multi_t, 'rb'))
        t_multi_list = []
        t_strbe_dict = pickle.load(open(pickle_strbe_t, 'rb'))
        t_strbe_list = []
        t_tors_dict = pickle.load(open(pickle_tors_t, 'rb'))
        t_tors_list = []
        t_oopb_dict = pickle.load(open(pickle_oopb_t, 'rb'))
        t_oopb_list = []

        for key, val in t_tot_E_dict.items():
            t_tot_E_list.append(val) # already multiplied by 4.184 to convert from Kcal/mol to kJ/mole
        for key, val in t_bond_dict.items():
            t_bond_list.append(val) # already multiplied by 4.184 to convert from Kcal/mol to kJ/mole
        for key, val in t_angle_dict.items():
            t_angle_list.append(val) # already multiplied by 4.184 to convert from Kcal/mol to kJ/mole
        for key, val in t_vdw_dict.items():
            t_vdw_list.append(val) # already multiplied by 4.184 to convert from Kcal/mol to kJ/mole
        for key, val in t_multi_dict.items():
            t_multi_list.append(val) # already multiplied by 4.184 to convert from Kcal/mol to kJ/mole
        for key, val in t_strbe_dict.items():
            t_strbe_list.append(val) # already multiplied by 4.184 to convert from Kcal/mol to kJ/mole
        for key, val in t_tors_dict.items():
            t_tors_list.append(val) # already multiplied by 4.184 to convert from Kcal/mol to kJ/mole
        for key, val in t_oopb_dict.items():
            t_oopb_list.append(val) # already multiplied by 4.184 to convert from Kcal/mol to kJ/mole

    else:
        t_tot_E_dict = {} # AMOEBA energy dict with frame numbers as keys
        t_tot_E_list = []
        t_bond_dict = {}
        t_bond_list = []
        t_angle_dict = {}
        t_angle_list = []
        t_vdw_dict = {}
        t_vdw_list = []
        t_multi_dict = {}
        t_multi_list = []
        t_strbe_dict = {}
        t_strbe_list = []
        t_tors_dict = {}
        t_tors_list = []
        t_oopb_dict = {}
        t_oopb_list = []

        ### run the tinker portion
        for i in range(num_frames_txyz):
            with open('outputA_%i.txt'%i) as out:
                outfile = out.readlines()
            outsplit = list_splitter(outfile)
            for j in range(len(outsplit)):
#                print('%i: %s'%(len(outsplit[j]), outsplit[j]))
                # skip all lines that are too short to contain useful info
                if len(outsplit[j]) > 1:
                    # find the "Total Potential Energy" line
                    if outsplit[j][0] == 'Total':
                        # skip lines with 7 entries to avoid wrong line
                        if len(outsplit[j]) != 7:
                            # specifically the line containing the energy value
                            if outsplit[j][3] == ':':
                                tot_e = float(outsplit[j][4])*4.184 # multiply by 4.184 to convert from Kcal/mol to kJ/mole
                                # key = krame, values = energy in Kcal/mole
                                t_tot_E_dict.update({i : tot_e})
                                print('frame %i total energy %f\n'%(i, tot_e))
                                t_tot_E_list.append(tot_e)
                    # find the "Bond Stretching" line
                    elif outsplit[j][0] == 'Bond':
                        bond_e = float(outsplit[j][2])*4.184 # multiply by 4.184 to convert from Kcal/mol to kJ/mole
                        # key = krame, values = energy in Kcal/mole
                        t_bond_dict.update({i : bond_e})
#                        print('frame %i bond energy %f\n'%(i, bond_e))
                        t_bond_list.append(bond_e)
                    elif outsplit[j][0] == 'Angle':
                        angle_e = float(outsplit[j][2])*4.184 # multiply by 4.184 to convert from Kcal/mol to kJ/mole
                        # key = krame, values = energy in Kcal/mole
                        t_angle_dict.update({i : angle_e})
#                        print('frame %i angle energy %f\n'%(i, angle_e))
                        t_angle_list.append(angle_e)
                    elif outsplit[j][0] == 'Van':
                        vdw_e = float(outsplit[j][3])*4.184 # multiply by 4.184 to convert from Kcal/mol to kJ/mole
                        # key = krame, values = energy in Kcal/mole
                        t_vdw_dict.update({i : vdw_e})
#                        print('frame %i angle energy %f\n'%(i, angle_e))
                        t_vdw_list.append(vdw_e)

                    elif outsplit[j][0] == 'Atomic':
                        multi_e = float(outsplit[j][2])*4.184 # multiply by 4.184 to convert from Kcal/mol to kJ/mole
                        # key = krame, values = energy in Kcal/mole
                        t_multi_dict.update({i : multi_e})
#                        print('frame %i multipole energy %f\n'%(i, multi_e))
                        t_multi_list.append(multi_e)
                    elif outsplit[j][0] == 'Stretch-Bend':
                        strbe_e = float(outsplit[j][1])*4.184 # multiply by 4.184 to convert from Kcal/mol to kJ/mole
                        # key = krame, values = energy in Kcal/mole
                        t_strbe_dict.update({i : strbe_e})
#                        print('frame %i str-bend energy %f\n'%(i, strbe_e))
                        t_strbe_list.append(strbe_e)
                    elif outsplit[j][0] == 'Torsional':
                        tors_e = float(outsplit[j][2])*4.184 # multiply by 4.184 to convert from Kcal/mol to kJ/mole
                        # key = krame, values = energy in Kcal/mole
                        t_tors_dict.update({i : tors_e})
#                        print('frame %i tors energy %f\n'%(i, tors_e))
                        t_tors_list.append(tors_e)
                    elif outsplit[j][0] == 'Out-of-Plane':
                        oopb_e = float(outsplit[j][2])*4.184 # multiply by 4.184 to convert from Kcal/mol to kJ/mole
                        # key = krame, values = energy in Kcal/mole
                        t_oopb_dict.update({i : oopb_e})
#                        print('frame %i tors energy %f\n'%(i, tors_e))
                        t_oopb_list.append(oopb_e)

        pickle.dump(t_tot_E_dict, open(pickle_tot_E_t, 'wb'))
        pickle.dump(t_bond_dict, open(pickle_bond_t, 'wb'))
        pickle.dump(t_angle_dict, open(pickle_angle_t, 'wb'))
        pickle.dump(t_vdw_dict, open(pickle_vdw_t, 'wb'))
        pickle.dump(t_multi_dict, open(pickle_multi_t, 'wb'))
        pickle.dump(t_strbe_dict, open(pickle_strbe_t, 'wb'))
        pickle.dump(t_tors_dict, open(pickle_tors_t, 'wb'))
        pickle.dump(t_oopb_dict, open(pickle_oopb_t, 'wb'))

    def energy_components(Sim):
        # Before using EnergyComponents, make sure each Force is set to a different group.
        EnergyTerms = OrderedDict()
        for i in range(Sim.system.getNumForces()):
            EnergyTerms[Sim.system.getForce(i).__class__.__name__] = Sim.context.getState(getEnergy=True,groups=2**i).getPotentialEnergy() / u.kilojoules_per_mole
        EnergyTerms['Potential'] = Sim.context.getState(getEnergy=True).getPotentialEnergy() / u.kilojoules_per_mole
        return EnergyTerms

    pickle_tot_E_o = 'openmm_tot_E_dict.p'
    pickle_bond_o = 'openmm_bond_dict.p'
    pickle_angle_o = 'openmm_angle_dict.p'
    pickle_vdw_o = 'openmm_vdw_dict.p'
    pickle_multi_o = 'openmm_multi_dict.p'
    pickle_strbe_o = 'openmm_strbe_dict.p'
    pickle_tors_o = 'openmm_tors_dict.p'
    pickle_oopb_o = 'openmm_oopb_dict.p'
    # AmoebaInPlaneAngleForce showing up now?
    pickle_ipaf_o = 'openmm_ipaf_dict.p'

    if os.path.exists(pickle_tot_E_o) and os.path.exists(pickle_bond_o) and os.path.exists(pickle_angle_o) and os.path.exists(pickle_vdw_o) and os.path.exists(pickle_multi_o) and os.path.exists(pickle_strbe_o) and os.path.exists(pickle_tors_o) and os.path.exists(pickle_oopb_o):
        o_tot_E_dict = pickle.load(open(pickle_tot_E_o, 'rb'))
        o_tot_E_list = []
        o_bond_dict = pickle.load(open(pickle_bond_o, 'rb'))
        o_bond_list = []
        o_angle_dict = pickle.load(open(pickle_angle_o, 'rb'))
        o_angle_list = []
        o_vdw_dict = pickle.load(open(pickle_vdw_o, 'rb'))
        o_vdw_list = []
        o_multi_dict = pickle.load(open(pickle_multi_o, 'rb'))
        o_multi_list = []
        o_strbe_dict = pickle.load(open(pickle_strbe_o, 'rb'))
        o_strbe_list = []
        o_tors_dict = pickle.load(open(pickle_tors_o, 'rb'))
        o_tors_list = []
        o_oopb_dict = pickle.load(open(pickle_oopb_o, 'rb'))
        o_oopb_list = []
        # added later
#        o_ipaf_dict = pickle.load(open(pickle_ipaf_o, 'rb'))
#        o_ipaf_list = []


        for key, val in o_tot_E_dict.items():
            o_tot_E_list.append(val) 
        for key, val in o_bond_dict.items():
            o_bond_list.append(val) 
        for key, val in o_angle_dict.items():
            o_angle_list.append(val) 
        for key, val in o_vdw_dict.items():
            o_vdw_list.append(val) 
        for key, val in o_multi_dict.items():
            o_multi_list.append(val) 
        for key, val in o_strbe_dict.items():
            o_strbe_list.append(val) 
        for key, val in o_tors_dict.items():
            o_tors_list.append(val) 
        for key, val in o_oopb_dict.items():
            o_oopb_list.append(val) 
        # added
#        for key, val in o_ipaf_dict.items():
#            o_ipaf_list.append(val) 


    else:
        o_tot_E_dict = {}
        o_tot_E_list = []
        o_bond_dict = {}
        o_bond_list = []
        o_angle_dict = {}
        o_angle_list = []
        o_vdw_dict = {}
        o_vdw_list = []
        o_multi_dict = {}
        o_multi_list = []
        o_strbe_dict = {}
        o_strbe_list = []
        o_tors_dict = {}
        o_tors_list = []
        o_oopb_dict = {}
        o_oopb_list = []
        # added
#        o_ipaf_dict = {}
#        o_ipaf_list = []


        ### run openmm portion
        for j in range(num_frames_txyz):
            
            # this portion is the old openmm section
            '''
            pdb = PDBFile('conf.pdb')
#            forcefield = ForceField('AMOEBA-urea.xml') #, 'tip3p.xml')
            forcefield = ForceField('AMOEBA-urea_NEW.xml') #, 'tip3p.xml')
            #unmatched_residues = forcefield.getUnmatchedResidues(pdb.topology)
            # use FB molecule object to set positions with xyz file
            molecule = Molecule('frame_%i.xyz'%j)
            # create an openmm list of Vec3's to update positions
            for x in range(len(molecule)):
                # FB molecule object
                xyz = molecule.xyzs[x]
                # openmm uses vec3
                xyz_vec3 = [Vec3(i[0], i[1], i[2]) for i in xyz]*angstrom
            
            # create system
            system = forcefield.createSystem(pdb.topology, nonbondedCutoff=1*nanometer)
            # give each force group its own index
            for ind, f in enumerate(system.getForces()):
                f.setForceGroup(ind)
            integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
            simulation = Simulation(pdb.topology, system, integrator)
            #simulation.context.setPositions(pdb.positions)
            simulation.context.setPositions(xyz_vec3)
            # create force component dictionary using above function
            Ecomps_OMM = energy_components(simulation)
            # print out force component dictionary
            printcool_dictionary(Ecomps_OMM, title="OpenMM energy components")
            '''
            
#            """
            # this portion is the new openmm portion (with AmoebaInPlaneAngleForce reinserted as AmoebaAngleForce, then deleted)
            pdb = PDBFile('conf.pdb')
            forcefield = ForceField('AMOEBA-urea_NEW.xml') #, 'tip3p.xml')
            #unmatched_residues = forcefield.getUnmatchedResidues(pdb.topology)
            # use FB molecule object to set positions with xyz file
            molecule = Molecule('frame_%i.xyz'%j)
            # create an openmm list of Vec3's to update positions
            for x in range(len(molecule)):
                # FB molecule object
                xyz = molecule.xyzs[x]
                # openmm uses vec3
                xyz_vec3 = [Vec3(i[0], i[1], i[2]) for i in xyz]*angstrom
            
            #oopb = AmoebaOutOfPlaneBendForce()
            #print(oopb)
            # create system
            system = forcefield.createSystem(pdb.topology, nonbondedCutoff=1*nanometer)
            
            fs = system.getForces()
            
            # determine index of AmoebaInPlaneAngleForce, which changes for each ff
            for i in range(len(fs)):
                print(fs[i])
                if "AmoebaInPlaneAngleForce" in str(fs[i]):
                    ipa_index = i
                # determine index of AmoebaAngleForce, which also changes for each ff
                elif "AmoebaAngleForce" in str(fs[i]):
                    angle_index = i
            
            '''
            print("ipa_index", ipa_index)
            
            # AmoebaAngleForce parameters (before any editing)
            print("BEFORE angle parameters")
            for i in range(fs[angle_index].getNumAngles()):
                print(fs[angle_index].getAngleParameters(i))
            '''
            
            # use the correct index to access AmoebaInPlaneAngleForce
            print(fs[ipa_index])
            print("number of in-plane angles", fs[ipa_index].getNumAngles())
            numinplane = fs[ipa_index].getNumAngles()
            for i in range(numinplane):
                print(fs[ipa_index].getAngleParameters(i))
            
            '''
            # access angle parameters in AmoebaInPlaneAngleForce
            test = fs[ipa_index].getAngleParameters(0)
            for i in range(len(test)):
                print("TEST", test[i])
            
            # angle is the 4th entry
            print("angle", str(test[4]))
            print(str(test[4]).split(' ')[0])
            # k is the 5th entry
            print("k", str(test[5]))
            '''
            
            # using the correct index for AmoebaAngleForce
            # reinsert AmoebaInPlaneAngleForce as an AmoebaAngleForce
            # example AmoebaInPlaneAngleForce: [1, 0, 2, 7, Quantity(value=152.31, unit=radian), Quantity(value=0.0235786067612, unit=kilojoule/(mole*radian**2))]
            # example AmoebaAngleForce: [1, 0, 2, Quantity(value=152.31, unit=degree), Quantity(value=0.0235786067612, unit=kilojoule/(mole*radian**2))]
            print(fs[angle_index])
            for i in range(numinplane):
                fs[angle_index].addAngle(fs[ipa_index].getAngleParameters(i)[0], fs[ipa_index].getAngleParameters(i)[1], fs[ipa_index].getAngleParameters(i)[2], float(str(fs[ipa_index].getAngleParameters(i)[4]).split(' ')[0]), float(str(fs[ipa_index].getAngleParameters(i)[5]).split(' ')[0]))
            
            print("AmoebaInPlaneAngleForces reinserted as AmoebaAngleForce")
            # remove in-plane angle force
            system.removeForce(ipa_index)
            print("AmoebaInPlaneAngleForces removed from the system")
            
            '''
            print("AFTER angle parameters")
            for i in range(fs[angle_index].getNumAngles()):
                print(fs[angle_index].getAngleParameters(i))
            '''
            # give each force group its own index
            for ind, f in enumerate(system.getForces()):
                f.setForceGroup(ind)
            integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
            simulation = Simulation(pdb.topology, system, integrator)
            #simulation.context.setPositions(pdb.positions)
            simulation.context.setPositions(xyz_vec3)
            # create force component dictionary using above function
            Ecomps_OMM = energy_components(simulation)
            # print out force component dictionary
            printcool_dictionary(Ecomps_OMM, title="OpenMM energy components")
#            """

            
            #for key, value in Ecomps_OMM.items():
            #    print(key, value)
            tot_e_o = float(Ecomps_OMM['Potential'])
            o_tot_E_dict.update({j : tot_e_o})
            o_tot_E_list.append(tot_e_o)
            bond_e_o = float(Ecomps_OMM['AmoebaBondForce'])
            o_bond_dict.update({j : bond_e_o})
            o_bond_list.append(bond_e_o)
            angle_e_o = float(Ecomps_OMM['AmoebaAngleForce'])
            o_angle_dict.update({j : angle_e_o})
            o_angle_list.append(angle_e_o)
            vdw_e_o = float(Ecomps_OMM['AmoebaVdwForce'])
            o_vdw_dict.update({j : vdw_e_o})
            o_vdw_list.append(vdw_e_o)
            multi_e_o = float(Ecomps_OMM['AmoebaMultipoleForce'])
            o_multi_dict.update({j : multi_e_o})
            o_multi_list.append(multi_e_o)
            strbe_e_o = float(Ecomps_OMM['AmoebaStretchBendForce'])
            o_strbe_dict.update({j : strbe_e_o})
            o_strbe_list.append(strbe_e_o)
            tors_e_o = float(Ecomps_OMM['PeriodicTorsionForce'])
            o_tors_dict.update({j : tors_e_o})
            o_tors_list.append(tors_e_o)
            oopb_e_o = float(Ecomps_OMM['AmoebaOutOfPlaneBendForce'])
            o_oopb_dict.update({j : oopb_e_o})
            o_oopb_list.append(oopb_e_o)
            # added
#            ipaf_e_o = float(Ecomps_OMM['AmoebaInPlaneAngleForce'])
#            o_ipaf_dict.update({j : ipaf_e_o})
#            o_ipaf_list.append(ipaf_e_o)

            
        pickle.dump(o_tot_E_dict, open(pickle_tot_E_o, 'wb'))
        pickle.dump(o_bond_dict, open(pickle_bond_o, 'wb'))
        pickle.dump(o_angle_dict, open(pickle_angle_o, 'wb'))
        pickle.dump(o_vdw_dict, open(pickle_vdw_o, 'wb'))
        pickle.dump(o_multi_dict, open(pickle_multi_o, 'wb'))
        pickle.dump(o_strbe_dict, open(pickle_strbe_o, 'wb'))
        pickle.dump(o_tors_dict, open(pickle_tors_o, 'wb'))
        pickle.dump(o_oopb_dict, open(pickle_oopb_o, 'wb'))
        # added
#        pickle.dump(o_ipaf_dict, open(pickle_ipaf_o, 'wb'))

    # plot of OpenMM vs AMOEBA energies (not too helpful)
    plt.figure()
    plt.plot(t_tot_E_list, o_tot_E_list)
    plt.title('Tinker vs OpenMM Energies')
    plt.ylabel('OpenMM kJ/mol')
    plt.xlabel('Tinker kJ/mol')
    plt.savefig('line_tinker_vs_openmm.pdf')
    plt.close()
    
    
    diff_tot_E_list = energy_diff_plot(t_tot_E_list, o_tot_E_list, 'totalE')
    # just run this for total E (for now)
    N_comp(100, diff_tot_E_list)
    
    diff_bond_list = energy_diff_plot(t_bond_list, o_bond_list, 'bond')
    diff_angle_list = energy_diff_plot(t_angle_list, o_angle_list, 'angle')
    diff_vdw_list = energy_diff_plot(t_vdw_list, o_vdw_list, 'vdw')
    diff_multi_list = energy_diff_plot(t_multi_list, o_multi_list, 'multi')
    diff_strbe_list = energy_diff_plot(t_strbe_list, o_strbe_list, 'strbe')
    diff_tors_list = energy_diff_plot(t_tors_list, o_tors_list, 'tors')
    diff_oopb_list = energy_diff_plot(t_oopb_list, o_oopb_list, 'oopb')
    

    '''
    # AMOEBA - OpenMM energy difference plot (kJ/mol)
    diff_list = []
    for i in range(len(t_list)):
        diff_list.append(t_list[i] - o_list[i])
#        print('tinker: %.2f, openmm: %.2f'%(t_list[i], o_list[i]))
    '''
    plt.figure()
    plt.plot(range(len(diff_bond_list)), diff_bond_list, label="Bond")
    plt.plot(range(len(diff_angle_list)), diff_angle_list, label="Angle")
    plt.plot(range(len(diff_strbe_list)), diff_strbe_list, label="Stretch-Bend")
    plt.plot(range(len(diff_tors_list)), diff_tors_list, label="Torsion")
    plt.plot(range(len(diff_oopb_list)), diff_oopb_list, label="Out-of-Plane-Bend")
    plt.plot(range(len(diff_multi_list)), diff_multi_list, label="Multipole")
    plt.plot(range(len(diff_vdw_list)), diff_vdw_list, label="VDW")

    # added
#    plt.plot(range(len(o_ipaf_list)), o_ipaf_list, label="In-Plane-Angle-Force")
    plt.title('Tinker - OpenMM Energy Difference vs Frame Number (Fixed)')
    plt.xlabel('MD frame')
    plt.ylabel('TINKER - OpenMM energy difference (kJ/mol)')
    plt.yticks(rotation=60)
    #plt.yscale("log")
    plt.legend()
    plt.savefig('diff_vs_frame_ALL_FIXED.pdf')
    plt.savefig('diff_vs_frame_ALL_FIXED.png', dpi=300)
    plt.close()


    plt.figure()
    plt.plot(range(len(t_angle_list)), t_angle_list, label="Tinker Angle")
    # added
#    plt.plot(range(len(o_ipaf_list)), o_ipaf_list, label="OpenMM In-Plane-Angle-Force")
    plt.plot(range(len(o_angle_list)), o_angle_list, label="OpenMM Angle")
    plt.title('Tinker and OpenMM Energies vs Frame Number')
    plt.xlabel('MD frame')
    plt.ylabel('Energy (kJ/mol)')
    #plt.yscale("log")
    plt.legend()
    plt.savefig('diff_angle-inpaf_ALL.pdf')
    
    plt.close()


    plt.figure()
    plt.plot(range(len(t_tot_E_list)), t_tot_E_list, label="Tinker Total Energy")
    # added
    plt.plot(range(len(o_tot_E_list)), o_tot_E_list, label="OpenMM Total Energy")
    plt.title('Tinker and OpenMM Energies vs Frame Number')
    plt.xlabel('MD frame')
    plt.ylabel('Energy (kJ/mol)')
    #plt.yscale("log")
    plt.legend()
    plt.savefig('diff_angle-inpaf_ALL.pdf')
    plt.close()

    plt.figure()
    plt.plot(range(len(diff_tot_E_list)), diff_tot_E_list)
    plt.title('Tinker - OpenMM Total Energy Difference vs Frame Number (Fixed)')
    plt.xlabel('MD frame')
    plt.ylabel('TINKER - OpenMM energy difference (kJ/mol)')
    plt.yticks(rotation=60)

    plt.savefig('diff_vs_frame_tot_E_FIXED.pdf')
    plt.savefig('diff_vs_frame_tot_E_FIXED.png', dpi=300)
    plt.close()



    totE_max = max(diff_tot_E_list)
    print(max(diff_tot_E_list))
    max_index = diff_tot_E_list.index(totE_max)
    per_diff = 100* abs(o_tot_E_list[max_index]-t_tot_E_list[max_index])/((o_tot_E_list[max_index]+t_tot_E_list[max_index])/2)
    print("per diff", per_diff)
        
if __name__ == '__main__':
    main()
                

import numpy as np
import argparse
import os
from os import path
import sys
import subprocess

def make_gauss_input(name, coord_name, option, charge=0, spin=1):
    ''' option #1 : create a Gaussian .com input file to run initial optimization
    option #2 : create a Gaussian .com input file to run second optimization
    option #3 : create a Gaussian .com input file to run frequency calculation
    '''
    coords_block = get_coordinates(coord_name)
    #coords_list = list_splitter(coords_block)
    #file_name = '%s.com'
    if option == 1:
        filename = '%s.com'%name
    elif option == 2:
        filename = '%s_re-opt.com'%name
    elif option == 3:
        filename = '%s_freq.com'%name
    with open(filename, 'w') as outfile:
        #outfile.write('%' + 'Chk=%s.chk\n'%name)
        if option == 1:
            outfile.write('%' + 'Chk=%s.chk\n'%name)
            outfile.write('# MP2/6-311G(1d,1p) Opt Density=MP2\n\n')
        elif option == 2:
            outfile.write('%' + 'Chk=%s_re-opt.chk\n'%name)
            outfile.write('# opt mp2/aug-cc-pVTZ density=mp2\n\n')
        elif option == 3:
            outfile.write('%' + 'Chk=%s_freq.chk\n'%name)
            outfile.write('# freq MP2/6-311G(1d,1p) density=mp2\n\n')

        outfile.write('remark line\n\n')
        outfile.write('%s %s\n'%(charge, spin))
        for i in range(2, len(coords_block)):
            outfile.write(coords_block[i])
        outfile.write('\n')

def make_gdma_input(name):
    ''' creates a GDMA input file
    '''
    with open('%s.gdmain'%name, 'w') as outfile:
        outfile.write('Title \"%s at MP2/6-311G(1d,1p)\"\n\n'%name)
        outfile.write('Verbose\n')
        outfile.write('Density MP2\n')
        outfile.write('File %s.fchk\n\n'%name)
        outfile.write('Angstrom\n')
        outfile.write('AU\n\n')
        outfile.write('Multipoles\n')
        outfile.write('Switch 0\n')
        outfile.write('Radius H 0.65\n')
        outfile.write('Limit 2\n')
        outfile.write('Start\n\n')
        outfile.write('Finish')

def make_poledit_input(name):
    ''' creates a TINKER poledit input file called input_poledit.txt
    as a reference of the inputs into that program. This program often does not
    identify equivalent atoms, and it is recommended that the user reviews the output
    parameter file name_new.key (using TINKER poledit) before moving on
    '''
    with open('input_poledit.txt', 'w') as outfile:
        outfile.write('1\n')
        outfile.write('%s.gdmaout\n'%name)
        # need this symbolic link to run everything
        outfile.write('amoeba09.prm\n\n')
        outfile.write('A\n\n\n')
        outfile.write('Y\n\n')

def make_prmedit_input(name, starting_atom_type):
    ''' creates a TINKER prmedit input file called input_prmedit.txt
    prmedit renumbers atom type indices for us. It is also capable of 
    renumbering atom class indices, but they are missing from the .key
    file at the point in time this is called. The function insert_atom_class_indices()
    inserts the atom class indices (correctly ordered).
    '''
    with open('input_prmedit.txt', 'w') as outfile:
        outfile.write('%s.key\n'%name)
        # option #3 is to edit atom types
        outfile.write('3\n')
        outfile.write('%s\n'%starting_atom_type)
        #outfile.write('%i'%starting_atom_class)

def insert_atom_class_indices(name, start_index):
    ''' this function does two main things: it copies the parameters
    from parameter.prm into name_new.key AND it inserts atom class
    indices starting from the input start_index.
    '''
    start_index = int(start_index)
    filename = "parameter.prm"
    with open(filename) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        with open(name + "_new.key", 'w') as edited:
            for i in range(len(filepieces)):
                if i == 0:
                    print("first line", filepieces[i])
                if len(filepieces[i]) != 0:
                    if filepieces[i][0] == 'atom':
                        filepieces[i].insert(2, str(start_index))
                        start_index += 1
                        editline = ""
                        for x in filepieces[i]:
                            editline += x + "    "
                        edited.write(editline + "\n")
                    else:
                        edited.write(filelines[i])
    

def edit_indices(filename, type_index, class_index):
    ''' not currently in use, but was used a starting point
    for some other functions
    '''
    """
    with open(filename) as fileblock:
        filelines = fileblock.readlines()
    # split lines by spaces
    filepieces = list_splitter(filelines)
    with open("new" + filename):
        for i in range(len(filepieces)):
            if i == 0:
                print("first line", filepieces[i])
#            if filepieces[i] == "atom":
    """
    with open(filename) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        with open("new" + filename, 'w') as edited:
            for i in range(len(filepieces)):
                if i == 0:
                    print("first line", filepieces[i])
                if filepieces[i][0] == 'atom':
                    newline = filepieces[i]
                    newline[1] = str(type_index + i)
                    newline[2] = str(class_index + i)
                    print("edited line %i"%i, newline)
                    editline = ""
                    for j in newline:
                        editline += j
                    edited.write(editline)
                else:
                    edited.write(filelines[i])

def make_potential_input(name, option):
    ''' makes input files for TINKER potential function calls (information from TINKER program)
    option #1 : creates an input file for Gaussian CUBEGEN
    option #2 : Obtain a QM Potential from a Gaussian CUBE file
    option #6 : Fits electrostatic potential parameters to a grid
    '''
    if option == 1:
        with open('input_potential_1.txt', 'w') as outfile:
            outfile.write('1\n')
            outfile.write('%s_re-opt.txyz\n'%name)
            outfile.write('%s_new.key'%name)
    elif option == 2:
        with open('input_potential_2.txt', 'w') as outfile:
            outfile.write('2\n')
            outfile.write('%s_re-opt.cube'%name)
    elif option == 6:
        with open('input_potential_6.txt', 'w') as outfile:
            outfile.write('6\n')
            outfile.write('%s_re-opt.txyz\n'%name)
            outfile.write('%s_new.key\n'%name)
            outfile.write('%s_re-opt.pot\n'%name)
            outfile.write('y\n\n')
            outfile.write('%s_new.key'%name)

def make_valence_input(name, option):
    ''' makes input files for TINKER valence function calls (information from TINKER program)
    option #1 : sets intitial values for valence parameters
    option #2 : compare QM and MM vibrational parameters (just for evaluation)
    option #3 : fits force parameters to the QM results
    '''
    if option == 1:
        with open('input_valence_1.txt', 'w') as outfile:
            outfile.write('1\n')
            outfile.write('%s_re-opt.txyz\n'%name)
            outfile.write('%s_re-opt.log\n'%name)
            outfile.write('%s_new_new.key'%name)
    elif option == 2 or option == 3:
        with open('input_valence_%i.txt'%option, 'w') as outfile:
            outfile.write('%i\n'%option)
            outfile.write('%s_re-opt.txyz\n'%name)
            outfile.write('%s_freq.log\n'%name)
            outfile.write('%s_new_new_new.key'%name)
            if option == 3:
                outfile.write('\n\n')

# function for splitting raw line-lists
def list_splitter(raw_list):
    ''' function for splitting raw line-lists into lines of separate elements
    '''
    split_list = []
    # split a list (by spaces) into separate elements
    for i in range(len(raw_list)):
        split_list.append(raw_list[i].split())
    return split_list

def get_coordinates(coord_name):
    ''' open coordinates file and split into separate lines
    '''
    with open(coord_name) as coors_file:
        coors_file_block = coors_file.readlines()
    return coors_file_block    

def line_prepender(filename, line):
    ''' function that inserts a line at the begining of a file
    '''
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)    

def combine_mp_parameters(old_params, opt_multipoles, name):
    save_file = "%s_new_new.key"%name
    # boolean to track copying of new multipole values
    copied = False
    with open(old_params) as fileblock1:
        filelines1 = fileblock1.readlines()
        # split lines by spaces
        filepieces1 = list_splitter(filelines1)
        with open(opt_multipoles) as fileblock2:
            filelines2 = fileblock2.readlines()
            # split lines by spaces
            filepieces2 = list_splitter(filelines2)

            with open(save_file, 'w') as edited:
                # iterate over each line of the original (old) parameter (.key) file
                for i in range(len(filepieces1)):
                    if i == 0:
                        print("first line", filepieces1[i])
                    if len(filepieces1[i]) != 0:
                        # identify "atom" OR "polarize" OR "FIX-MULTIPOLE" lines and copy those
                        if filepieces1[i][0] == "atom" or filepieces1[i][0] == "polarize" or filepieces1[i][0] == "FIX-MONOPOLE":
                            edited.write(filelines1[i])
                        # otherwise, skip the old multipole values and copy the new multipole values
                        elif filepieces1[i][0] == "multipole" and copied == False:
                            # iterate over each line of the optimized multipoles file and write them into the edit
                            for j in range(len(filelines2)):
                                edited.write(filelines2[j])
                            copied = True
                    # if line empty, include a line return to preserve formatting
                    else:
                        edited.write("\n")
    print("newest parameter file is %s"%save_file)

def combine_val_parameters(old_params, new_params, name):
    ''' combine valence parameter results from Tinker valence commands into a single file
    '''
    save_file = "%s_new_new_new.key"%name
    # boolean to track copying of new parameter values
    keep_copying = False
    # boolean to skip copying the Urey-Bradley parameters (which are zero for most of these systems)
    copy_urey = False
    with open(old_params) as fileblock1:
        filelines1 = fileblock1.readlines()
        # split lines by spaces
        filepieces1 = list_splitter(filelines1)
        with open(new_params) as fileblock2:
            filelines2 = fileblock2.readlines()
            # split lines by spaces
            filepieces2 = list_splitter(filelines2)

            with open(save_file, 'w') as edited:
                # iterate over and copy each line of the original (old) parameter (.key) file
                for i in range(len(filelines1)):
                    edited.write(filelines1[i])
                #edited.write("\n")
                # iterate over each line in the new parameter file (output_valence_1.txt) and copy the relevant lines
                for j in range(len(filepieces2)):
                    # skip empty lines
                    if len(filepieces2[j]) != 0:
                        # only start copying once the first "Estimated" appears, then copy every line afterward
                        if filepieces2[j][0] == "Estimated":
                            # add line returns after Estimated to aid the formatting
                            #edited.write("\n")
                            keep_copying = True
                        if keep_copying == True:
                            #print(filelines2[j])
                            # skip the ureybrad parameter lines
                            if copy_urey == False and filepieces2[j][0] != "ureybrad":
                                edited.write(filelines2[j])
#                            else:
#                                edited.write(filelines2[j])
                    # if line empty, include a line return to preserve formatting (only after encountering the first instance of "Estimated"
                    elif keep_copying == True:
                        edited.write("\n")
    print("newest parameter file is %s"%save_file)


def combine_final_parameters(old_params, new_params, name):
    ''' combine parameter results from final Tinker commands into a single, final file
    '''
    save_file = "%s_final.key"%name
    # boolean to track copying of new multipole values
    copy_new = False
    # boolean to track copying of old multipole values
    copy_old = True
    with open(old_params) as fileblock1:
        filelines1 = fileblock1.readlines()
        # split lines by spaces
        filepieces1 = list_splitter(filelines1)
        with open(new_params) as fileblock2:
            filelines2 = fileblock2.readlines()
            # split lines by spaces
            filepieces2 = list_splitter(filelines2)

            with open(save_file, 'w') as edited:
                # iterate over each line of the original (old) parameter (.key) file
                for i in range(len(filepieces1)):
                    if i == 0:
                        print("first line", filepieces1[i])
                    if len(filepieces1[i]) != 0:
                        # update copy_old boolean once the old bond and angle section is skipped
                        if filepieces1[i][0] == "Estimated" and filepieces1[i][1] == "Stretch-Bend":
                            # resume copying the old parameter lines once past the old bond and angle section
                            copy_old = True
                        # once "Estimate Bond" is reached, begin copying new bond and angle parameters
                        if filepieces1[i][0] == "Estimated" and filepieces1[i][1] == "Bond":
                            # iterate over each line of the new parameter file and write them into the edit
                            for j in range(len(filepieces2)):
                                if len(filepieces2[j]) != 0:
                                    # start copying the new parameters once "# Results" is detected
                                    if filepieces2[j][0] == "bond" and copy_new == False:
                                        edited.write("#\n")
                                        edited.write("# Results of Valence Parameter Fitting\n")
                                        edited.write("#\n")
                                        copy_new = True
                                    # copy new bond and angle parameters
                                    if copy_new == True:
                                        edited.write(filelines2[j])
                            # stop copying until "Estimated Stretch-Bend" is found in the old parameter file
                            copy_old = False
                        # by default, copy old parameter lines (except old bond and angle parameters)
                        elif copy_old == True:
                            edited.write(filelines1[i])
                    # include line returns to preserve formatting
                    else:
                        edited.write("\n")
    print("final parameter file is %s"%save_file)

def insert_params(old_params, name):
    ''' insert the parameter preamble into a force field file
    '''
    save_file = "%s_new_new_2.key"%name
    # boolean to track copying of new multipole values
    copy_new = False
    # boolean to track copying of old multipole values
    copy_old = True
    with open(old_params) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces1 = list_splitter(filelines)
        with open(save_file, 'w') as edited:
            # first line file name (needed for conversion to OpenMM later on
            edited.write("forcefield AMOEBA-%s\n\n"%name)
            # insert each of the Amoeba setup parameters
            edited.write("bond-cubic              -2.55\n")
            edited.write("bond-quartic            3.793125\n")
            edited.write("angle-cubic             -0.014\n")
            edited.write("angle-quartic           0.000056\n")
            edited.write("angle-pentic            -0.0000007\n")
            edited.write("angle-sextic            0.000000022\n")
            edited.write("opbendtype              ALLINGER\n")
            edited.write("opbend-cubic            -0.014\n")
            edited.write("opbend-quartic          0.000056\n")
            edited.write("opbend-pentic           -0.0000007\n")
            edited.write("opbend-sextic           0.000000022\n")
            edited.write("torsionunit             0.5\n")
            edited.write("vdwtype                 BUFFERED-14-7\n")
            edited.write("radiusrule              CUBIC-MEAN\n")
            edited.write("radiustype              R-MIN\n")
            edited.write("radiussize              DIAMETER\n")
            edited.write("epsilonrule             HHG\n")
            edited.write("dielectric              1.0\n")
            edited.write("polarization            MUTUAL\n")
            edited.write("vdw-12-scale            0.0\n")
            edited.write("vdw-13-scale            0.0\n")
            edited.write("vdw-14-scale            1.0\n")
            edited.write("vdw-15-scale            1.0\n")
            edited.write("mpole-12-scale          0.0\n")
            edited.write("mpole-13-scale          0.0\n")
            edited.write("mpole-14-scale          0.4\n")
            edited.write("mpole-15-scale          0.8\n")
            edited.write("polar-12-scale          0.0\n")
            edited.write("polar-13-scale          0.0\n")
            edited.write("polar-14-scale          1.0\n")
            edited.write("polar-15-scale          1.0\n")
            edited.write("polar-12-intra          0.0\n")
            edited.write("polar-13-intra          0.0\n")
            edited.write("polar-14-intra          0.5\n")
            edited.write("polar-15-intra          1.0\n")
            edited.write("direct-11-scale         0.0\n")
            edited.write("direct-12-scale         1.0\n")
            edited.write("direct-13-scale         1.0\n")
            edited.write("direct-14-scale         1.0\n")
            edited.write("mutual-11-scale         1.0\n")
            edited.write("mutual-12-scale         1.0\n")
            edited.write("mutual-13-scale         1.0\n")
            edited.write("mutual-14-scale         1.0\n\n")

            # iterate over each line of the original (old) parameter (.key) file and copy as is
            for i in range(len(filelines)):
                edited.write(filelines[i])
    print("newest parameter file is %s"%save_file)

def find_atom_types(old_params, new_params, name, starting_index):
    ''' uses an input file to find the atom types to reorder atom types the name_re-opt.txyz
    When babel is called to convert to .txyz it uses atom type assignments that do not apply to these custom setups.
    
    '''
    save_file = "%s_re-opt.txyz"%name
    # boolean to track copying of new multipole values
    copy_new = False
    # boolean to track copying of old multipole values
    copy_old = False
    # dict to hold atom index and associated atom types from the poledit output file
    atypes = {}
    with open(old_params) as fileblock1:
        filelines1 = fileblock1.readlines()
        # split lines by spaces
        filepieces1 = list_splitter(filelines1)
        with open(new_params) as fileblock2:
            filelines2 = fileblock2.readlines()
            # split lines by spaces
            filepieces2 = list_splitter(filelines2)

            with open(save_file, 'w') as edited:
                # iterate over each line of the original (old) parameter (.key) file
                for i in range(len(filepieces1)):
                    if i == 0:
                        print("first line", filepieces1[i])
                    if len(filepieces1[i]) != 0:
                        # update copy_old boolean once the old bond and angle section is skipped
                        if filelines1[i] == " Atom Type and Local Frame Definition for Each Atom :\n":
                            print("MADE IT HERE")
                            # resume copying the old parameter lines once past the old bond and angle section
                            copy_old = True
                        elif filelines1[i] == " Final Atomic Multipole Moments after Regularization :\n":
                            # halt copying (dictionary creation) once this line is encountered
                            copy_old = False

                        # skip empty lines
                        if len(filepieces1[i]) != 0 and copy_old == True:
                            # skip "Atom" line to get to index lines
                            if filepieces1[i][0] != "Atom":
                                #print("current line" + filelines1[i])
                                a_index = filepieces1[i][0]
                                a_type = filepieces1[i][1]
                                #print("atom index %s with atom type %s"%(a_index, a_type))
                                atypes.update({int(a_index) : int(a_type)})
                # iterate over each line of the new parameter file and write them into the edit
                for j in range(len(filepieces2)):
                    # just copy the first line
                    if j == 0:
                        edited.write(filelines2[j])
                    else:
                        edited_line = ""
                        for k in range(len(filepieces2[j])):
                            # if not the 5th (starting from zero) element of the list, include the original line fragment
                            if k != 5:
                                edited_line += "    %s"%filepieces2[j][k]
                            else:
                                # the new type is based on the old type which starts from 1 (hence the subtraction)
                                new_type = atypes[int(filepieces2[j][0])] + int(starting_index) -1 
                                edited_line += "    %i"%new_type
                        edited.write(edited_line + "\n")
                        # iterate over the segments of each line, updating the atom types
                        
    print("newest parameter file is %s"%save_file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--coordinates", help="starting coordinates for molecule of interest")
    parser.add_argument("-d", "--charge", default=0, help="charge of molecule")
    parser.add_argument("-s", "--spin", default=1, help="spin multiplicity of molecule")
    parser.add_argument("-n", "--name", help="name of molecule")
    parser.add_argument("-p", "--part", help="part of process to run, run part 1 before part 2")
    parser.add_argument("-t", "--type_index", help="starting atom type index")
    parser.add_argument("-i", "--class_index", help="starting atom class index")
    args = parser.parse_args()

    # Programs needed: Gaussian (v 16), Tinker (v 8.7.1), GDMA (v 2.3), Babel (2.3.1) 

    # start of code (in the urea folder)
    # PART ONE
    # start out with a traditional (not-TINKER) name.xyz file and  a symbolic link to amoeba09.prm
    if args.part == "1":
        print("***** running part 1 *****")
        #make_gauss_input(args.name=name, args.coordinates=coord_name, args.charge=charge, args.spin=spin) 
        # create the gaussian input
        make_gauss_input(args.name, args.coordinates, 1, args.charge, args.spin) 

        # run initial gaussian optimization (command_1)
        print("running initial Guassian optimization")
        os.system("g16 %s.com"%args.name)

        # run formchk (command_2)
        print("running Gaussian formchk")
        os.system("formchk %s.chk %s.fchk"%(args.name, args.name))
        # create gdma input 
        make_gdma_input(args.name)
        # run gdma program (command_3)
        print("running the GDMA program")
        os.system("gdma < %s.gdmain > %s.gdmaout"%(args.name, args.name))
        # create Tinker poledit input (command_4)
        make_poledit_input(args.name)
        ### add error check here if no amoeba09.prm file in folder!!!
        # run Tinker poledit 
        print("running TINKER poledit")
        os.system("poledit < input_poledit.txt > output_poledit.txt")
        # edit starting atom types and atom classes indices of .key file (command_4.5)
        # choline+ starts atom type indices at 500 and atom classes at 200, so the H-bond donors start from 600 and 300 respectively
        #starting_atom_type = 500
        #starting_atom_class = 200
        print("Before running part 2, review %s.key to make sure that the atom types were correctly assigned.\nSometimes TINKER doesn't correcly identify all equivalent atoms.\nUse TINKER poledit to \"Condense Symmetric Atoms to Equivalent Types\""%args.name)
        print("If changes are made, make sure to name the correct key file: %s_new.key. \nRe-running TINKER poledit manually should save the updated .key file as %s_new.key_2 (or higher), also in that file it will append the new parameters after the old ones, make sure the old ones are removed."%(args.name, args.name)) 
        print("output_poledit.txt is a required file for part 2 of this auto code. Capturing the entire screen output from manually running TINKER poledit is not compatable with the function")

    if args.part == "2":
        print("***** running part 2 *****")
        make_prmedit_input(args.name, args.type_index)
        print("running TINKER prmedit")
        # run prmedit (creates parameter.prm) !!! this removes atom classes instead of editing them!!!
        os.system("prmedit < input_prmedit.txt > output_prmedit.txt")

        # make a copy of parameter.prm
        #if path.exists("parameter.prm") == True:
        #    os.system("cp parameter.prm parameter_orig.prm")
        # manually insert missing atom class indices
        print("copying parameter.prm into %s_new.key and correcting atom type and atom class indices"%args.name)
        insert_atom_class_indices(args.name, args.class_index)

        # OUTDATED: before running part 2, run part 1, then edit parameter.prm to include atom class incides    
        #print("***** running (OLD) part 2 *****")
#        """
        # rename parameter.prm (command_5)
        #if path.exists("parameter.prm") == True:
        #    os.system("mv parameter.prm %s_new.key"%args.name)
        # get initial optimization coordinates -> name_opt.xyz
        print("running babel to convert .log to .xyz")
        os.system("babel %s.log %s_opt.xyz"%(args.name, args.name))
        # create input for next Gaussian optimization
        make_gauss_input(args.name, "%s_opt.xyz"%args.name, 2, args.charge, args.spin) 
        # run second gaussian optimization (command_6)
        print("running second Gaussian optimization")
        os.system("g16 %s_re-opt.com"%args.name)
#        """

        # run formchk on the gaussion re-optimization (command_7)
        print("running Gaussian formchk (on second opt results)")
        os.system("formchk %s_re-opt.chk %s_re-opt.fchk"%(args.name, args.name))
        # run babel to get the re-opt coordinates (command_8)
        print("running babel to convert .fchk to .txyz")
        os.system("babel -ifchk %s_re-opt.fchk -otxyz %s_re-opt_WRONG_TYPES.txyz"%(args.name, args.name))
        print("the next step before running part 3, is to edit %s_re-opt.txyz to have to the correct atom types"%args.name)
#    if args.part == "t":
        print("made it to function call")
        # trying out the function that takes care of the atom types
        find_atom_types("output_poledit.txt", "%s_re-opt_WRONG_TYPES.txyz"%args.name, args.name, args.type_index)

    # OUTDATED: before running part 3, run part 1 then part 2 and edit name_re-opt.txyz to have the correct atom type
#    if args.part == "3":
        #print("***** running (OLD) part 3 *****")
        # create input file for potential option #1
        make_potential_input(args.name, 1)
        # run tinker potential option #1 (command_9)
        print("running TINKER potential option #1")
        os.system("potential < input_potential_1.txt > output_potential_1.txt")
        # run cubegen (command_10)
        print("running Gaussian cubegen")
        os.system("cubegen 0 potential=mp2 %s_re-opt.fchk %s_re-opt.cube -5 h < %s_re-opt.grid"%(args.name, args.name, args.name))
        # create input file for potential option #2
        make_potential_input(args.name, 2)
        # run tinker potential option #2 (command_11)
        print("running TINKER potential option #2")
        os.system("potential < input_potential_2.txt > output_potential_2.txt")
        # insert "FIX-MONOPOLE" into .key file produced in step 5 (command_12)
        #???i
        line_prepender("%s_new.key"%args.name , "FIX-MONOPOLE") 
        # create input file for potential option #6
        print("running TINKER potential option #6")
        make_potential_input(args.name, 6)
        # run tinker potential option #6 (command_13)
        os.system("potential < input_potential_6.txt > output_potential_6.txt")

        # combine new fitted multipole results from name_re-opt.key with atom definitions AND polarize values from name_new.key to create name_new_new.key
        combine_mp_parameters("%s_new.key"%args.name, "%s_re-opt.key"%args.name, args.name)
        # insert amoeba parameters before running valence later on -> name_new_new_2.key
        insert_params("%s_new_new.key"%args.name, args.name)

    # OUTDATED: before running part 4, run part 1 then part 2 then part 3 and combine the multipole values from name_re-opt.key with the monopole values from name_new.key to create name_new_new.key (command_14)
        #print("***** running (OLD) part 4 *****")
        # create input file for valence option #1
        make_valence_input(args.name, 1)
        # run tinker valence option #1 (command_15)
        print("running TINKER valence option #1")
        os.system("valence < input_valence_1.txt > output_valence_1.txt")
        # combine new valence parameter values with existing parameter values into new .key file
        combine_val_parameters("%s_new_new_2.key"%args.name, "output_valence_1.txt", args.name)

    # OUTDATED: before running part 5, run parts 1-4 in sequence, then combine the parameters in output_valence_1.txt with those in name_new_new.key to form name_new_new_new.key (command_16)
        #print("***** running (OLD) part 5 *****")
        # create .xyz coordinate file to be used to create a gaussian input file
        print("running babel to convert .txyz to .xyz")
        os.system("babel %s_re-opt.txyz %s_re-opt.xyz"%(args.name, args.name))
        # make gaussian input file for frequency calculation
        make_gauss_input(args.name, "%s_re-opt.xyz"%args.name, 3, args.charge, args.spin)
        # run gaussian frequency calculation (command_17)
        print("running Gaussian frequency calculation")
        os.system("g16 %s_freq.com"%args.name)
        # create input file for valence option #2
        make_valence_input(args.name, 2)
        # run tinker valence option #2 which is just a comparison (command_18)
        print("running TINKER valence option #2")
        os.system("valence < input_valence_2.txt > output_valence_2.txt")
        # create input file for valence option #3
        make_valence_input(args.name, 3)
        # run tinker valence option #3 (functional command_19)
        print("running TINKER valence option #3")
        os.system("valence < input_valence_3.txt > output_valence_3.txt")
         
        # OUTDATED: after running parts 1-5 the new bond and angle parameters will be named after the input coordinates (name_re-opt.txyz), so these will be named name_re-opt.key_2 because we have name_re-opt.key from before
        # OUTDATED: the next step is to replace those parameters in name_new_new_new.key to form name_final.key. Note, the multipole parameters are included in name_re-opt.key_2 even though they didn't change from earlier
        combine_final_parameters("%s_new_new_new.key"%args.name, "%s_re-opt.key_2"%args.name, args.name)

if __name__ == '__main__':
    main()



import argparse
import os
from os import path

# fb/FB = ForceBalance

# function for splitting raw line-lists
def list_splitter(raw_list):
    ''' function for splitting raw line-lists into lines of separate elements
    '''
    split_list = []
    # split a list (by spaces) into separate elements
    for i in range(len(raw_list)):
        split_list.append(raw_list[i].split())
    return split_list

# testing out a different method
def insert_openmm_fb_parameters_v2(openmm_file, param):
    ''' inserts desired fb parameters into an OpenMM force field file
    '''
    intermediate_file = "intermediate_" + openmm_file
    new_file = "new_" + openmm_file
    intermediate_path = os.getcwd() + "/" + intermediate_file 
    # if there is not already an intermediate_file, create one
    if os.path.exists(intermediate_path) == False:
        # copy openmm file
        os.system("cp %s %s"%(openmm_file, intermediate_file))
    # read input openmm file
    with open(intermediate_file) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        enter_loop = False
        # create new edited file
        with open(new_file, 'w') as edited:
            bond_edit = False
            #current_line = i
            #j = 1
            if param == "Bond":
                parameters = "length, k"
            elif param == "Angle":
                parameters = "angle1, k"
                bond_edit = True
            elif param == "Vdw":
                parameters = "sigma, epsilon"
                bond_edit = True
            elif param == "Proper":
                parameters = "k1, k2, k3"
                bond_edit = True
            elif param == "StretchBend":
                parameters = "k1, k2"
                bond_edit = True
            #elif param == "Multipole":
            # this part is a bit clunky, need to be able to separate Multipole and Polarize under AmoebaMultipoleForce
            elif param == "Polarize":
                parameters = "polarizability"
                bond_edit = True
            # remove "parameterize" from every line
            #elif param == "all":
            # create line_start beginning of line check
            line_start= "<%s"%param
            for i in range(len(filepieces)):
                # special measure to dodge lines in Residues that also start with "<Bond"
                if filepieces[i][0] == "<AmoebaBondForce":
                    bond_edit = True
                # for each line that starts with line_start
                if (filepieces[i][0] == line_start and bond_edit == True):
                    #print("made it here")
                    #current_line = i+j
                    #start_line = i
                    print("current line", i)
                    # create new_line string to be added into file
                    new_line = filelines[i].split("/>")[0] + "parameterize=\"%s\" />\n"%parameters
                    # write the new_line into the edited file
                    edited.write(new_line)
                    print(new_line)
                    #j +=1
                else:
                    edited.write(filelines[i])

    # copy new_file onto intermediate file 
    os.system("cp %s %s"%("new_" + openmm_file, "intermediate_" + openmm_file))
    # delete old intermediate file
    
    return

# testing out a different method
def insert_tinker_fb_parameters_v2(tinker_file, param):
    ''' inserts desired fb parameters into an Tinker force field file
    '''
    intermediate_file = "intermediate_" + tinker_file
    new_file = "new_" + tinker_file
    intermediate_path = os.getcwd() + "/" + intermediate_file 
    # if there is not already an intermediate_file, create one
    if os.path.exists(intermediate_path) == False:
        # copy openmm file
        os.system("cp %s %s"%(tinker_file, intermediate_file))
    # read input openmm file
    with open(intermediate_file) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        enter_loop = False
        # create new edited file
        with open(new_file, 'w') as edited:
            #elif param == "all":
            for i in range(len(filepieces)):
                #print("current line", filelines[i])
                #print("len(filepieces[i])", len(filepieces[i]))
                # for each line that starts with line_start
                if (len(filepieces[i]) > 0) and (filepieces[i][0] == param):
                    radian = 57.2957795130
                    #print("made it here")
                    #current_line = i+j
                    #start_line = i
                    print("current line", i)
                    # create new_line string to be added into file
                    #new_line = filelines[i].split("/>")[0] + "parameterize=\"%s\" />\n"%parameters
                    if param == "bond":
                        new_line = filelines[i].split("\n")[0] + " # EVAL 3 PARM[\"Bond/k/%s.%s\"]/418.4 4 PARM[\"Bond/length/%s.%s\"]*10\n"%(filepieces[i][1], filepieces[i][2], filepieces[i][1],filepieces[i][2])  
                    elif param == "angle":
                        new_line = filelines[i].split("\n")[0] + " # EVAL 4 PARM[\"Angle/k/%s.%s.%s\"]/784.6085 5 PARM[\"Angle/angle1/%s.%s.%s\"]\n"%(filepieces[i][1], filepieces[i][2], filepieces[i][3], filepieces[i][1], filepieces[i][2], filepieces[i][3])  
                    elif param == "vdw":
                        new_line = filelines[i].split("\n")[0] + " # EVAL 2 PARM[\"Vdw/sigma/%s\"]*10 3 PARM[\"Vdw/epsilon/%s\"]/4.184\n"%(filepieces[i][1], filepieces[i][1])  
                    # found the conversion factors from the conversion code
                    elif param == "strbnd":
                        conversion = (4.184*10/radian)**-1
                        new_line = filelines[i].split("\n")[0] + " # EVAL 4 PARM[\"StretchBend/k1/%s.%s.%s\"%f] 5 PARM[\"StretchBend/k2/%s.%s.%s\"%f]\n"%(filepieces[i][1], filepieces[i][2], filepieces[i][3], conversion, filepieces[i][1], filepieces[i][2], filepieces[i][3], conversion)  
                    elif param == "torsion":
                        # parameterizing just the k values
                        new_line = filelines[i].split("\n")[0] + " # EVAL 5 PARM[\"Proper/k1/%s.%s.%s.%s\"/2.092] 8 PARM[\"Proper/k1/%s.%s.%s.%s\"/2.092] 11 PARM[\"Proper/k1/%s.%s.%s.%s\"/2.092]\n"%(filepieces[i][1], filepieces[i][2], filepieces[i][3], filepieces[i][4], filepieces[i][1], filepieces[i][2], filepieces[i][3], filepieces[i][4], filepieces[i][1], filepieces[i][2], filepieces[i][3], filepieces[i][4])  
                    # not quite sure this is correct
                    elif param == "polarize":
                        new_line = filelines[i].split("\n")[0] + " # EVAL 2 PARM[\"Polarize/polarizability/%s\"]*1000\n"%(filepieces[i][1])  

                    # write the new_line into the edited file
                    edited.write(new_line)
                    print("edited line:\n", new_line)
                    #j +=1
                else:
                    edited.write(filelines[i])

    # copy new_file onto intermediate file 
    os.system("cp %s %s"%("new_" + tinker_file, "intermediate_" + tinker_file))
    # delete old intermediate file
    
    return


def insert_openmm_fb_parameters(openmm_file, param):
    ### doesn't currently work properly, moving to version 2 of this function

    intermediate_file = "intermediate_" + openmm_file
    new_file = "new_" + openmm_file
    intermediate_path = os.getcwd() + "/" + intermediate_file 
    # if there is not already an intermediate_file, create one
    if os.path.exists(intermediate_path) == False:
        # copy openmm file
        os.system("cp %s %s"%(openmm_file, intermediate_file))
    # read input openmm file
    with open(intermediate_file) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        enter_loop = False
        # create new edited file
        with open(new_file, 'w') as edited:
            for i in range(len(filepieces)):
                current_line = i
                j = 1
                # check for parameter name instance in line
                if filepieces[i][0] == "<" + param:
                    print("test!!!", filelines[i])
                    edited.write(filelines[i])
                    stop = False
                    # boolean to track whether the editing loop has been entered
                    enter_loop = True
                    # increment over the lines following the parameter line while condition is True
                    while stop == False:
                        # editing specific to AmoebaBondForce
                        if param == "AmoebaBondForce":
                            line_start = "<Bond"
                            parameters = "length, k"
                        elif param == "AmoebaAngleForce":
                            line_start = "<Angle"
                            parameters = "angle1, k"
                        elif param == "AmoebaVdwForce":
                            line_start = "<Vdw"
                            parameters = "sigma, epsilon"
                        elif param == "PeriodicTorsionForce":
                            line_start = "<Proper"
                            parameters = "k1, k2, k3"
                        elif param == "AmoebaStretchBendForce":
                            line_start = "<StretchBend"
                            parameters = "k1, k2"
                        elif param == "Multipole":
                            line_start = "<Multipole"
                        # this part is a bit clunky, need to be able to separate Multipole and Polarize under AmoebaMultipoleForce
                        elif param == "AmoebaMultipoleForce":
                            line_start = "<Polarize"
                            parameters = "polarizability"
                        # remove "parameterize" from every line
                        elif param == "all":
                            line_start = "<Vdw"
                        # for each line that starts with line_start
                        if filepieces[i+j][0] == line_start:
                            #print("made it here")
                            current_line = i+j
                            start_line = i
                            print("start_line", start_line)
                            # create new_line string to be added into file
                            new_line = filelines[current_line].split("/>")[0] + "parameterize=\"%s\" />\n"%parameters
                            # write the new_line into the edited file
                            edited.write(new_line)
                            print(new_line)
                            #j +=1
                        else:
                            num_edits = j
                            stop = True
                            #break
                        j +=1
                # if no parameter instance, write line (unless it was edited in the edit loop
                else:
                    # if edit loop has not been entered, write unedited line
                    if enter_loop == False:
                        edited.write(filelines[i])
                    # if edit loop has been entered
                    else:
                        # skip over edited lines
                        exclude_line = i in range(start_line +1, start_line +num_edits)
                        if exclude_line == False:
                            edited.write(filelines[i])  
    # copy new_file onto intermediate file 
    os.system("cp %s %s"%("new_" + openmm_file, "intermediate_" + openmm_file))
    # delete old intermediate file
    
    return
    

def insert_tinker_fb_parameters(tinker_file, param):
    ### doesn't currently work properly, moving to version 2 of this function

    intermediate_file = "intermediate_" + tinker_file
    new_file = "new_" + tinker_file
    intermediate_path = os.getcwd() + "/" + intermediate_file 
    # if there is not already an intermediate_file, create one
    if os.path.exists(intermediate_path) == False:
        # copy openmm file
        os.system("cp %s %s"%(tinker_file, intermediate_file))
    # read input openmm file
    with open(intermediate_file) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        enter_loop = False
        # create new edited file
        with open(new_file, 'w') as edited:
            for i in range(len(filepieces)):
                current_line = i
                j = 0
                #print("current_line: ", filelines[i])
                #print("len(filepieces[i])", len(filepieces[i]))
                #print("filepieces[i][0]", filepieces[i][0])
                # check for parameter name instance in line
                #if filepieces[i][0] ==  param:
                if (len(filepieces[i]) > 0) and (filepieces[i][0] == param):
                    print("test!!!", filelines[i])
                    edited.write(filelines[i])
                    stop = False
                    # boolean to track whether the editing loop has been entered
                    enter_loop = True
                    # increment over the lines following the parameter line while condition is True
                    while stop == False:
                        # for each line that starts with line_start (same as param)
                        print("line no", i+j)
                        print("entry line", filelines[i+j])
                        # condition added to skip empty lines
                        if (len(filepieces[i+j]) > 0) and (filepieces[i+j][0] == param):
                            #print("made it here")
                            current_line = i+j
                            start_line = i
                            print("start_line", start_line)
                            # create the new_line depending on the param
                            if param == "bond":     
                                # create new_line string to be added into file (including conversion factors)
                                new_line = filelines[current_line].split("\n")[0] + " # EVAL 3 PARM[\"Bond/k/%s.%s\"]/418.4 4 PARM[\"Bond/length/%s.%s\"]*10\n"%(filepieces[current_line][1], filepieces[current_line][2], filepieces[current_line][1],filepieces[current_line][2])  

                                # write the new_line into the edited file
                                edited.write(new_line)
                                print(new_line)
                                #j +=1
                            else:
                                num_edits = j
                                stop = True
                                #break

                        j +=1
                # if no parameter instance, write line (unless it was edited in the edit loop
                else:
                    # if edit loop has not been entered, write unedited line
                    if enter_loop == False:
                        edited.write(filelines[i])
                    # if edit loop has been entered
                    else:
                        # skip over edited lines
                        exclude_line = i in range(start_line +1, start_line +num_edits)
                        if exclude_line == False:
                            edited.write(filelines[i])  
    # copy new_file onto intermediate file 
    os.system("cp %s %s"%("new_" + tinker_file, "intermediate_" + tinker_file))
    # delete old intermediate file
    
    return

def remove_all_fb_parameters_openmm(openmm_file):
    ''' removes all fb parameters from an OpenMM force field file
    '''
    intermediate_file = "intermediate_" + openmm_file
    new_file = "new_" + openmm_file
    intermediate_path = os.getcwd() + "/" + intermediate_file 
    # if there is not already an intermediate_file, create one
    if os.path.exists(intermediate_path) == False:
        # copy openmm file
        os.system("cp %s %s"%(openmm_file, intermediate_file))
    # read input openmm file
    with open(intermediate_file) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        enter_loop = False
        current_line = 0
        # create new edited file
        with open(new_file, 'w') as edited:
            for line in filelines:
                if "parameterize" in line:
                    print("editing line entered at line %i"%current_line)
                    new_line = filelines[current_line].split("parameterize")[0] + "/>\n"
                    # write the new_line into the edited file
                    edited.write(new_line)
                    print(new_line)
                else:
                    edited.write(filelines[current_line])
                current_line += 1 
    # copy new_file onto intermediate file 
    os.system("cp %s %s"%("new_" + openmm_file, "intermediate_" + openmm_file))
    # delete old intermediate file


def remove_all_fb_parameters_tinker(tinker_file):
    ''' removes all fb parameters from a Tinker force field file
    '''
    intermediate_file = "intermediate_" + tinker_file
    new_file = "new_" + tinker_file
    intermediate_path = os.getcwd() + "/" + intermediate_file 
    # if there is not already an intermediate_file, create one
    if os.path.exists(intermediate_path) == False:
        # copy openmm file
        os.system("cp %s %s"%(tinker_file, intermediate_file))
    # read input openmm file
    with open(intermediate_file) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        enter_loop = False
        current_line = 0
        # create new edited file
        with open(new_file, 'w') as edited:
            for line in filelines:
                if "# EVAL" in line:
                    print("editing line entered at line %i"%current_line)
                    new_line = filelines[current_line].split("#")[0] + "\n"
                    # write the new_line into the edited file
                    edited.write(new_line)
                    print(new_line)
                else:
                    edited.write(filelines[current_line])
                current_line += 1 
    # copy new_file onto intermediate file 
    os.system("cp %s %s"%("new_" + tinker_file, "intermediate_" + tinker_file))
    # delete old intermediate file



def remove_openmm_parameters(openmm_file, param):
    ''' removes a certain fb parameter from an OpenMM force field file
    '''
    intermediate_file = "intermediate_" + openmm_file
    new_file = "new_" + openmm_file
    intermediate_path = os.getcwd() + "/" + intermediate_file 
    # if there is not already an intermediate_file, create one
    if os.path.exists(intermediate_path) == False:
        # copy openmm file
        os.system("cp %s %s"%(openmm_file, intermediate_file))
    # read input openmm file
    with open(intermediate_file) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        enter_loop = False
        # create new edited file
        with open(new_file, 'w') as edited:
            for i in range(len(filepieces)):
                current_line = i
                j = 1
                # check for parameter name instance in line
                if filepieces[i][0] == "<" + param:
                    print("test!!!", filelines[i])
                    edited.write(filelines[i])
                    stop = False
                    # boolean to track whether the editing loop has been entered
                    enter_loop = True
                    # increment over the following lines while condition is True
                    while stop == False:
                        # editing specific to AmoebaBondForce
                        if param == "AmoebaBondForce":
                            line_start = "<Bond"
                        elif param == "AmoebaAngleForce":
                            line_start = "<Angle"
                        elif param == "AmoebaVdwForce":
                            line_start = "<Vdw"
                        elif param == "PeriodicTorsionForce":
                            line_start = "<Proper"
                        elif param == "AmoebaStretchBendForce":
                            line_start = "<StretchBend"
                        elif param == "Multipole":
                            line_start = "<Multipole"
                        elif param == "Polarize":
                            line_start = "<Polarize"
                        # remove "parameterize" from every line
                        elif param == "all":
                            line_start = "<Vdw"
                        # for each line that starts with line_start
                        if filepieces[i+j][0] == line_start:
                            #print("made it here")
                            current_line = i+j
                            start_line = i
                            print("start_line", start_line)
                            # create new_line string to be added into file
                            new_line = filelines[current_line].split("parameterize")[0] + "/>\n"
                            # write the new_line into the edited file
                            edited.write(new_line)
                            print(new_line)
                            j +=1
                        else:
                            num_edits = j
                            stop = True
                            #break
                # if no parameter instance, write line (unless it was edited in the edit loop
                else:
                    # if edit loop has not been entered, write unedited line
                    if enter_loop == False:
                        edited.write(filelines[i])
                    # if edit loop has been entered
                    else:
                        # skip over edited lines
                        exclude_line = i in range(start_line +1, start_line +num_edits)
                        if exclude_line == False:
                            edited.write(filelines[i])  
    # copy new_file onto intermediate file 
    os.system("cp %s %s"%("new_" + openmm_file, "intermediate_" + openmm_file))
    # delete old intermediate file


def remove_tinker_parameters(tinker_file, param):
    ''' removes a certain fb parameter from a Tinker force field file
    '''
    # this function removes fb parameter information according to param input
    # this function may not be used over the "all" version
    intermediate_file = "intermediate_" + tinker_file
    new_file = "new_" + tinker_file
    intermediate_path = os.getcwd() + "/" + intermediate_file 
    # if there is not already an intermediate_file, create one
    if os.path.exists(intermediate_path) == False:
        # copy openmm file
        os.system("cp %s %s"%(tinker_file, intermediate_file))
    # read input openmm file
    with open(intermediate_file) as fileblock:
        filelines = fileblock.readlines()
        # split lines by spaces
        filepieces = list_splitter(filelines)
        enter_loop = False
        # create new edited file
        with open(new_file, 'w') as edited:
            for i in range(len(filepieces)):
                current_line = i
                j = 1
                # check for parameter name instance in line
                if filepieces[i][0] == "<" + param:
                    print("test!!!", filelines[i])
                    edited.write(filelines[i])
                    stop = False
                    # boolean to track whether the editing loop has been entered
                    enter_loop = True
                    # increment over the following lines while condition is True
                    while stop == False:
                        # editing specific params
                        if param == "AmoebaBondForce":
                            line_start = "bond"
                        elif param == "AmoebaAngleForce":
                            line_start = "angle"
                        elif param == "AmoebaVdwForce":
                            line_start = "vdw"
                        elif param == "PeriodicTorsionForce":
                            line_start = "torsion"
                        elif param == "AmoebaStretchBendForce":
                            line_start = "strbnd"
                        elif param == "Multipole":
                            line_start = "multipole"
                        elif param == "Polarize":
                            line_start = "polarize"
                        # remove "parameterize" from every line (still needs work)
                        elif param == "all":
                            line_start = "<Vdw"
                        # for each line that starts with line_start
                        if filepieces[i+j][0] == line_start:
                            #print("made it here")
                            current_line = i+j
                            start_line = i
                            print("start_line", start_line)
                            # create new_line string to be added into file
                            # remove everything after the #
                            new_line = filelines[current_line].split("#")[0] + "/>\n"
                            # write the new_line into the edited file
                            edited.write(new_line)
                            print(new_line)
                            j +=1
                        else:
                            num_edits = j
                            stop = True
                            #break
                # if no parameter instance, write line (unless it was edited in the edit loop
                else:
                    # if edit loop has not been entered, write unedited line
                    if enter_loop == False:
                        edited.write(filelines[i])
                    # if edit loop has been entered
                    else:
                        # skip over edited lines
                        exclude_line = i in range(start_line +1, start_line +num_edits)
                        if exclude_line == False:
                            edited.write(filelines[i])  
    # copy new_file onto intermediate file 
    os.system("cp %s %s"%("new_" + openmm_file, "intermediate_" + openmm_file))
    # delete old intermediate file

    return



def main():
    # this part could use some refinement, such as using an input file to control the parameter choices
    # instead of the current method of commenting lines in and out
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--openmm", help="OpenMM forcefield file")
    parser.add_argument("-t", "--tinker", help="TINKER forcefield file")
    args = parser.parse_args()

    #remove_openmm_parameters(args.openmm, "AmoebaBondForce")
    #remove_openmm_parameters(args.openmm, "AmoebaAngleForce")
    #remove_openmm_parameters(args.openmm, "AmoebaVdwForce")
    remove_all_fb_parameters_openmm(args.openmm)
    '''
    insert_openmm_fb_parameters(args.openmm, "AmoebaBondForce")
    insert_openmm_fb_parameters(args.openmm, "AmoebaAngleForce")
    insert_openmm_fb_parameters(args.openmm, "PeriodicTorsionForce")
    insert_openmm_fb_parameters(args.openmm, "AmoebaStretchBendForce")
    insert_openmm_fb_parameters(args.openmm, "AmoebaVdwForce")
    #insert_openmm_fb_parameters(args.openmm, "Multipole")
    insert_openmm_fb_parameters(args.openmm, "AmoebaMultipoleForce")
    '''
    insert_openmm_fb_parameters_v2(args.openmm, "Bond")
    insert_openmm_fb_parameters_v2(args.openmm, "Angle")
    insert_openmm_fb_parameters_v2(args.openmm, "Proper")
    insert_openmm_fb_parameters_v2(args.openmm, "StretchBend")
    insert_openmm_fb_parameters_v2(args.openmm, "Vdw")
    #insert_openmm_fb_parameters_v2(args.openmm, "Multipole")
    insert_openmm_fb_parameters_v2(args.openmm, "Polarize")

    # test removing all the fb parameters
    remove_all_fb_parameters_tinker(args.tinker)

    # insert fb parameters
    insert_tinker_fb_parameters_v2(args.tinker, "bond")
    insert_tinker_fb_parameters_v2(args.tinker, "angle")
    insert_tinker_fb_parameters_v2(args.tinker, "vdw")
    insert_tinker_fb_parameters_v2(args.tinker, "strbnd")
    insert_tinker_fb_parameters_v2(args.tinker, "torsion")
    insert_tinker_fb_parameters_v2(args.tinker, "polarize")


    # remove openmm intermediate file
    os.system("rm %s"%("intermediate_" + args.openmm))

    # remove tinker intermediate file
    os.system("rm %s"%("intermediate_" + args.tinker))

if __name__ == '__main__':
    main()


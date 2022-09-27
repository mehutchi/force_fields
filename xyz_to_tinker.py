import argparse

# function for splitting raw line-lists
def list_splitter(raw_list):
    """function for splitting raw line-lists into lines of separate elements
    """
    split_list = []
    # split a list (by spaces) into separate elements 
    for i in range(len(raw_list)):
        split_list.append(raw_list[i].split())
        #split_list.append(raw_list[i].strip(" "))

    return split_list


def main():
    """ This code converts a .xyz file to a tinker .txyz file based on an input with the same atom ordering.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--xyz_coords_file", help="input traditional xyz file to be adjusted to tinker (t)xyz")
    parser.add_argument("-t", "--txyz_coords_file", help="input tinker (t)xyz file to be used as reference to convert traditional xyz file to (t)xyz")
    parser.add_argument("-o", "--output_file_name", default="all.arc", help="output file name, default is \"all.arc\"")
    args = parser.parse_args()

    #input_2 = 'CCU_1_2_first_frame.xyz'


    input_file = args.xyz_coords_file#'all.xyz'
    print("xyz file to be converted to tinker: %s"%args.xyz_coords_file)
    tinker_ref = args.txyz_coords_file#'CCU_1_2_first_frame.txyz'
    print("txyz file to be used as reference: %s"%args.txyz_coords_file)
    #num_atoms = 0
    # open the all.xyz coordinate file
    with open(input_file, 'r') as raw_input_file:
        input_split = raw_input_file.readlines()
        num_atoms = int(input_split[0][0])
        # open the converted and corrected TINKER xyz file
        with open(tinker_ref, 'r') as raw_ref_file:
            ref_lines = raw_ref_file.readlines()
            ref_split = list_splitter(ref_lines)
            output_file = args.output_file_name
            counter = 0
            # open the file to be written
            with open(output_file, 'w') as output:
                # iterate over each line of all.xyz (MD trajectory)
                for line in range(len(input_split)):
                    # find each line that contains a semicolon (common to terachem-produced .xyz files, might need to make this more robust)
                    if ';' in input_split[line]:
                        output.write('%i frame %i\n'%(num_atoms, counter))
                        counter += 1
                        # iterate over the next n lines, where n = num_atoms
                        for i in range(1, num_atoms +1):
                            current_line = line + i
                            #print("current line %s"%ref_split[current_line])
                            #add_on_list = ref_split[current_line][5:]
                            # this part is tricky, .xyz has an extra space that .txyz doesn't, but the index i matches the line on the reference
                            # split the current line (by whitespace)
                            add_on_list = ref_lines[i].split()
                            #print('frame %i'%counter)
                            string_add_on = ''
                            # iterate from 5th element of the split list onward and create a spaced out string (of atom class and connectivity) to convert all.xyz to all.arc (.txyz trajectory file)
                            for x in add_on_list[5:]:
                            #    print('add %s'%x)
                                x_add = '   ' + x
                                string_add_on += x_add
                            #print("add on portion %s"%string_add_on) 
                            output.write('%i '%i + input_split[current_line].rstrip('\n') + string_add_on + '\n')


if __name__ == '__main__':
    main()

import argparse

""" Frequencies.dat file atom indices start from one, but frequency indices start from zero. Documentation for vdata.txt formatting below from Lee-Ping:
"""

#==========================================#
#| File containing vibrational modes from |#
#|      experiment/QM for ForceBalance    |#
#|                                        |#
#| Octothorpes are comments               |#
#| This file should be formatted like so: |#
#| (Full XYZ file for the molecule)       |#
#| Number of atoms                        |#
#| Comment line                           |#
#| a1 x1 y1 z1 (xyz for atom 1)           |#
#| a2 x2 y2 z2 (xyz for atom 2)           |#
#|                                        |#
#| These coords will be actually used     |#
#|                                        |#
#| (Followed by vibrational modes)        |#
#| ...                                    |#
#| v (Eigenvalue in wavenumbers)          |#
#| dx1 dy1 dz1 (Eigenvector for atom 1)   |#
#| dx2 dy2 dz2 (Eigenvector for atom 2)   |#
#| ...                                    |#
#| (Empty line is optional)               |#
#| v (Eigenvalue)                         |#
#| dx1 dy1 dz1 (Eigenvector for atom 1)   |#
#| dx2 dy2 dz2 (Eigenvector for atom 2)   |#
#| ...                                    |#
#| and so on                              |#
#|                                        |#
#| Please list freqs in increasing order  |#
#==========================================#

""" Frequencies.dat stores the frequencies in columns, but vdata.txt needs them orgainized into rows.
"""

# function for splitting raw line-lists
def list_splitter(raw_list):
    ''' function for splitting raw line-lists into lines of separate elements
    '''
    split_list = []
    # split a list (by spaces) into separate elements
    for i in range(len(raw_list)):
        split_list.append(raw_list[i].split())
    return split_list

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--freq_file", default="Frequencies.dat", help="Takes a Frequencies.dat file produced by a terachem frequencies calculation and converts it to a ForceBalance input file.")
    parser.add_argument("-c", "--coords_file", default="CentralGeometry.initcond.xyz", help="coordinates at COM used to obtain the Frequencies.dat file.")
    args = parser.parse_args()

    new_file = "vdata.txt"
    with open(new_file, 'w') as edited:
        with open(args.freq_file) as fileblock:
            filelines = fileblock.readlines()
            filepieces = list_splitter(filelines)

            with open(args.coords_file) as fileblock2:
                filelines2 = fileblock2.readlines()
                for line in filelines2:
                    edited.write(line)

                num_atoms = int(filepieces[0][2])
                num_vib_modes = int(filepieces[1][3])
                mode = 0
                col = 0
                # iterate over each line in the file
                for i in range(len(filepieces)):
                    # check if the line has 4 elements
                    if len(filepieces[i]) == 4:
                        # iterate over each column
                        for col in range(4):
                            # check if col-th element of the line is the correct mode
                            if filepieces[i][col] == str(mode):
                                # write the magnitude of the vib mode
                                edited.write("\n" + filepieces[i+1][col] + "\n")
                                # counter to track current row (skipping mode index, mode magnitude, and ----- lines)
                                row_count = 3
                                # iterate over the next 3 lines (coordinates)
                                for j in range(num_atoms):
                                    edited.write("    %s    %s    %s\n"%(filepieces[i+row_count][1+col], filepieces[i+row_count+1][col], filepieces[i+row_count+2][col]))
                                    row_count += 3
                                mode += 1



if __name__ == '__main__':
    main()

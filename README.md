# force_fields

## automatic_generate_tinker_v3.py
This code runs in two parts to create an AMOEBA force field from a small molecule input.

## fb_parameter_adjust.py
This code prepares AMOEBA and OpenMM force field files for usage within ForceBalance (FB). First any exisiting FB keywords are removed from the files, then the desired keywords are added in. Choice of keywords depends on the desired optimization targets for FB, ex: bond, angle, VDW, etc. 

## vib_mode_format_converter.py
This short code takes a Frequencies.dat file calculated by Terachem vibrational frequencies calculation and converts it into vdata.txt, a required ForceBalance input file.

## xyz_to_tinker.py
This short code converts a .xyz input file into a Tinker .txyz file using a reference with the same atom ordering. The .txyz file format uses the same coordinates as .xyz, except that it explicitly lists atom connectivity and requires labeling of atom types. Exisiting molecular file conversion software can easily convert .txyz to .xyz, but it is not as suited for converting in the opposite direction when using custom AMOEBA force fields with unique atom types and classes.

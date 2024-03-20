# force_fields

## automatic_generate_tinker_v3.py
This code creates an AMOEBA force field (using Python to run Tinker, Gaussian, and Open Babel) from a small molecule input.

The Ponder group has a tutorial which instructs how to use their Tinker software to create a custom AMOEBA force field for any desired small molecule. This code automates that tutorial into two parts, with only one small intervention needed to verify atoms of equivalent symmetry. 

## fb_parameter_adjust.py
This code prepares AMOEBA and OpenMM force field files for usage within ForceBalance (FB). FB is a force field optimizer and requires placement of certain keywords and information to identify targets.

To prepare force fields for usage in FB, this code first removes any existing FB keywords from the files, then the desired keywords are added in. Choice of keywords depends on the desired optimization targets for FB, ex: bond, angle, VDW, etc. 

## vib_mode_format_converter.py
This short code takes a Frequencies.dat file (calculated by Terachem vibrational frequencies calculation) and converts it into vdata.txt, a required ForceBalance input file.

The main obstacle is that Frequencies.dat stores the frequencies in columns, but vdata.txt needs them orgainized into rows.

## xyz_to_tinker.py
This short code converts a .xyz input file into a Tinker .txyz file using a reference with the same atom ordering. The .txyz file format uses the same coordinates as .xyz, except that it explicitly lists atom connectivity and requires labeling of atom types. Existing molecular file conversion software can easily convert .txyz to .xyz, but it is not as suited for converting in the opposite direction when using custom AMOEBA force fields with unique atom types and classes.

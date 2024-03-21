# force_fields
The goal of this project is to efficiently and with low user error create custom small-molecule AMOEBA force fields using the TINKER software.

## automatic_generate_tinker_v3.py
This code creates an AMOEBA force field (using Python to run Tinker, Gaussian, and Open Babel) from a small molecule input.

The Ponder group has a tutorial which instructs how to use their Tinker software to create a custom AMOEBA force field for any desired small molecule. This code automates that approximately 25-step tutorial into two "one-click" parts, with only one small intervention needed in-between to verify atoms of equivalent symmetry. 

## tinker_vs_openmm_FIXED.py

The goal of this project was to take the AMOEBA force fields generated by automatic_generate_tinker_v3.py and optimize their parameters using ForceBalance (FB). However, it is convenient for FB to carry the parameter opimization forward using funcationally identical AMOEBA and OpenMM force fields. This is due to the fact that some calculations are easier or just not present in one engine verses the other.

Peter Eastman from the OpenMM Project kindly provided a AMOEBA-OpenMM force field conversion code which I adjusted to our needs. Unfortunately, I encounted many pitfalls attempting to create an exact copy of the custom AMOEBA force fields within OpenMM. This code was created to troubleshoot the observed energy differences when comparing "identical" force fields in AMOEBA and OpenMM. The following python plots illustrate part of that troubleshooting process:

The first troubleshooting measure was to compare the atomic distance matrices in a trial molecular dynamics simulation using a Urea force field in both AMOEBA and OpenMM. The purpose was to check if any extreme geometry contributed to the energy differences displayed by the two force fields which should have appeared identical.

![urea_dist_matrix](https://github.com/mehutchi/force_fields/assets/20996215/db751d71-7303-4232-8ae2-cccb95724efc)

In this figure, the upper left is the distance matrix for the 0th frame in the trial MD simulation
for urea. Lower left is an average distance matrix of the 100 least energetically different MD
frames between AMOEBA-Tinker and AMOEBA-OpenMM. Before averaging, each frame
had the 0th difference matrix subtracted (to highlight deviation from the planar starting
point). Lower right is the same as lower left, except using the 100 most energetically different
MD frames. All distances in ˚Angstrom, redder color indicates a larger distance.

## fb_parameter_adjust.py
This code prepares AMOEBA and OpenMM force field files for usage within ForceBalance (FB). FB is a force field optimizer and requires placement of certain keywords and information to identify targets.

To prepare force fields for usage in FB, this code first removes any existing FB keywords from the files, then the desired keywords are added in. Choice of keywords depends on the desired optimization targets for FB, ex: bond, angle, VDW, etc. 

## vib_mode_format_converter.py
This short code takes a Frequencies.dat file (calculated by Terachem vibrational frequencies calculation) and converts it into vdata.txt, a required ForceBalance input file.

The main obstacle is that Frequencies.dat stores the frequencies in columns, but vdata.txt needs them orgainized into rows.

## xyz_to_tinker.py
This short code converts a .xyz input file into a Tinker .txyz file using a reference with the same atom ordering. (Note that Tinker files often use the .xyz extension even though they differ from traditional XYZ files, but I've found it useful to include a 't' in the extension when using both to help distinguish them)

The .txyz file format uses the same coordinates as .xyz, except that it explicitly lists atom connectivity and requires labeling of atom types. Existing molecular file conversion software can easily convert .txyz to .xyz, but it is not as suited for converting in the opposite direction when using custom AMOEBA force fields with uniquely indexed atom types and classes.

Aside from a desired .xyz file to convert, a reference .txyz file with the correct atom classes and attachments is needed. Atom ordering must be respected for this to work, but is not an issue if the molecules are not manipulated than re-saved.

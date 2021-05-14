# StructureAnalysis  
Collection of Scripts to Analyze 2º and 3º RNA Structure with a focus on internucleotide measurements.  
Scripts were designed to work with outputs of ShapeJumper https://github.com/Weeks-UNC/ShapeJumper_V1  
However, they can be applied to work with other internucleotide pair data sets.  

[3DDistanceHistogram.py](https://github.com/Weeks-UNC/StructureAnalysis/blob/master/README.md#3DDistanceHistogrampy)  
[draw_3D_RNA_contacts.py](https://github.com/Weeks-UNC/StructureAnalysis/blob/master/README.md#draw_3D_RNA_contactspy)  
[AnnotateSecondaryStructure.py](https://github.com/Weeks-UNC/StructureAnalysis/blob/master/README.md#AnnotateSecondaryStructurepy)  
[contactDistanceHistogram.py](https://github.com/Weeks-UNC/StructureAnalysis/blob/master/README.md#contactDistanceHistogrampy)  
[deletionLengthDistribution.py](https://github.com/Weeks-UNC/StructureAnalysis/blob/master/README.md#deletionLengthDistributionpy)  
[Common input files types used](https://github.com/Weeks-UNC/StructureAnalysis/blob/master/README.md#file-descriptions)  

### General Requirements
All scripts require python 2.7  

### Installation
clone this package to your home directory with this command:  
`git clone https://github.com/Weeks-UNC/StructureAnalysis.git`  

# 3DDistanceHistogram.py  
Maps internucleotide pairs onto provided pdb structure. Returns a histogram of the 3D distances between these pairs.  
The pairs mapped are limited to those with a frequency above a given percentile.  
The distance histrogram is plotted against a histogram of all pairwise distances in the given pdb structure.  


### Requirements
This script requires the python 2 modules matplotlib, numpy, and Biopython.  
The python script RNAtools2 is also required, available here: https://github.com/Weeks-UNC/RNATools  

### Required Inputs
**PDB File** - Structure file in pdb format. It is critical that the residue numbering matches that used in the deletion file.    
**Deletion File** - list of nucleotide pairs and their frequency. See *File Description* section at bottom of readme for more detail.  
**Gene** - Input a gene or sequence name to limit nucleotide pairs plotted to just those from a specific sequence name contained in the deletion file.
**Percentile** - Decimal value between 0 and 100. Limits pairs plotted to those with a frequency above this given percentile.  

### Execution Instructions
`python 3DDistanceHistogram.py TertiaryStructure.pdb DeletionFile.txt geneName PercentileNumber`  
Example: `python 3DDistanceHistogram.py 3DHS.pdb RNasePDeletions.txt rnasep 97.5`  

### Optional Inputs
`-u,--upperPercentile FLOAT` - Apply an upper cutoff as well, so all deletions are below this percentile. For example if percentile is set at 90 and upperPercentile is set to 95 only deletions between the 90th and 95th percentiles will be included in the histogram.  
`-c,--color PYPLOTCOLOR` - Change the color of the histogram from the default of black. Requires a color name compatible with pyplot. Acceptable color names here: https://matplotlib.org/stable/gallery/color/named_colors.html  
`--experimentalOnly` - Only plot the histogram of distances from the deletion file, not the set of all possible distances in the PDB.    
`--inclusivePercentile` - By default the percentile cutoff is calculated from only the set of nucleotide pairs present in the given PDB file. This option changes the percentile cutoff calculation to include all nucleotide pairs. This option will only affect the output if nucleotides are missing from the PDB.
  `--topDels INT` - Limits pairs plotted to given number rather than percentile eg. most frequent 100 pairs.    
`--lengthAdjust SEQUENCELENGTH` - Create background histogram from random Nt pairs that reflect the 1D length of deletions. Requires length of sequence.  
`--shapeFile FILE.shape` - Used with other options to limit pairs plotted by shape reactivities provided in file.  
`--shapeCutoff FLOAT` - Limits pairs plotted to those where both nucleotides have a shape reactivity above cutoff input. Use with --shapeFile.  
`--oneNTSelection ` - Changes function of --shapeCutoff so that only one of two nucleotides in pair must exceed cutoff. Use with --shapeCutoff and --shapeFile.  
`--shift5nt INT` - Shift all deletion start sites, the 5' end, by the input value. Can be either positive or negative.  
`--shift3nt INT` - Shift all deletion stop sites, the 3' end, by the input value. Can be either positive or negative.  
`--centralPoint` - Rather than calculating distances between 2' hydroxyl groups, calculate between central point of the nucleobase. Useful when plotting UV induced or psoralen crosslinking interactions.    
`--corrChi` - Allows deletion file input to be in chi square correlation file format. See https://github.com/Weeks-UNC/RingMapper  
`--StrucFilt STRUCTURE.ct CONTACTDISTANCE` - Limits nucleotide pairs plotted to those with a contact distance above given cutoff. See https://rna.urmc.rochester.edu/Text/File_Formats.html for description of CT file format.  

### Output Description
Generates a pdf containing a histogram of the 3D distances between nucleotide pairs. This is in a black outline and overlaid onto a solid light purple histogram showing all possible 3D pairwise distances in the provided PDB structure.   

# draw_3D_RNA_contacts.py
Draws connections between set of nucleotide pairs on a pdb file using the pymol distance command.  

### Requirements
PyMol - Should be executable on the command line.  

### Required Inputs
**Contact File** - Text file of nucleotide pairs and frequencies. See *File Description* at bottom of readme for more detail.  
**PDB File** - Structure file in pdb format. It is critical that the residue numbering matches that used in the deletion file.  

### Execution Instructions
This script is designed to run via PyMol on the command line.  
`pymol -qcr draw_3D_RNA_contacts.py -- contactFile.txt structure.pdb `  
The `-qcr` flag executes the script via pymol without opening the GUI or printing a startup message.  

### Output Description
Generates a pymol session file of the pdb with contacts drawn. The ouput file name is the contact file name with a .pse suffix. For example if the input contact file was RNasePContacts.txt the output would be RNasePContacts.pse  

# AnnotateSecondaryStructure.py  
Annotates a varna or xrna secondary structure file and generates a SVG file ouput. Available annotations include coloring nucleotides, changing font size and drawing contacts between nucleotides, among others.  
While some of these annotations can be accomplished in the VARNA and/or xRNA programs, their implementation here is meant to be quicker and more user friendly.  

### Requirements
This script requires the python 2 modules matplotlib, numpy, and Biopython.  

### Required Inputs
**Input** - Varna or xRNA file containing a secondary structure. To obtain these files use the programs[VARNA](http://varna.lri.fr/) or [xRNA](http://rna.ucsc.edu/rnacenter/xrna/). Make sure to save your structure in the .varna or .xrna format. If the input is an xRNA file you must use the -x option. 
**Output** - Filename for output file. Should end in .svg.  

### Execution Instructions
`python AnnotateSecondaryStructure.py inputStructure.VARNA outputFileName.svg`  

### Optional Inputs
`-x,--xrna` - Use this flag if the input file is in the xRNA format.  
`--ct SecondaryStructure.ct` - Colors structure by its agreement with provided reference structure. Base pairs present in both structures are colored green. Those present in only the CT file are colored red and those only in the varna are colored purple.    
`-e,--enzyme`  
`-d,--diff` - Changes chemical probing type to differential, coloring cutoffs +=0.3.  
`-c,--colors ColorFile` - Colors circle behind nucleotide letter according to numbers in input file. 0 for white, 1 for red, 2 for grey, 3 for blue, 4 for orange, 5 for purple, 6 for black. Color file should be a text file with a number for each nucleotide on its own line. The number in line 1 colors nt 1, line 2 nt 2 and so on.  
`-l,--lines "ContactFile.txt LowCutoff HighCutoff"` - Draws lines between nucleotides specified in contact file. Contacts with rates above the low cutoff value are colored orange, values above the high cutoff are colored red. Contacts with rates below the low cutoff are not plotted. See *File Description* section at bottom of readme for in depth description of contact file format.
`-p,--percentileLines "DeletionFile.txt Percentile"` - Draws lines between nucleotide pairs specified in deletion file. Only draws lines for deletions with rates above given percentile. Multiple percentiles, up to 5, may be entered in ascending order. The order of line coloring, corresponding to percentiles entered, is cyan, read, orange, green, blue. See *File Description* section at bottom of readme for in depth description of deletion file format.  
`-r,--rawCountLines "DeletionFile.txt CountX"`- Similar to --percentLines but limits to the top X most frequent deletions. Up to 5 cutoffs may be entered, in ascending order.  
`-cl,--categoryLines "DeletionFile.txt Percentile"` - Similar to --percentileLines except the deletion file should have a fifth column where every column is assigned a category number, 0-2. Deletions in category 0 are colored green, 1 are violet and 2 are orange. Unlike --percentileLines, only one percentile may be entered.  
`--distanceLines "DeletionFile.txt 3DStructure.pdb Percentile` - Similar to --percentileLines except lines are colored by their 3D distance, as determined from the input pdb file. Lines with distances less than 10 angstroms are colored pure green, longer than 50 are pure blue. Distances between 10 and 50 are colored on a proportional gradient from green to blue. If one or both nts in a pair are missing from the pdb, the line is colored grey. Unlike --percentileLines, only one percentile may be entered.  
`--centralPoint` - Use with --distanceLines. Calculates 3D distances between nucleotide pairs by the central point of the nucleobase rather than the default 2' hydroxyl group.  
`-s,--shape reactivityData.shape` - Colors nucleotides by shape reactivity. Reactivites above 0.8 are colored red, above 0.4 are colored orange, below -0.4 are colored black and -999 values are colored grey.  
`--dms DMSReactivityData.shape` - Similar to the --shape option. Sets color limits to 0.2 and 0.4  
`-o,--offset OffsetNum` - Adds the input offset number to nucleotide numbering. Eg. if the offset number is 9, nucleotide 1 would be labeled as 10.  
`--offsetStart StartNum ` - Only add the offset to the start number nucleotide and beyond. Allows a gap to be introduced to nucleotide numbering. Use with --offset option.  
`-n,--nonCanonicalPairs pairsFile.txt` - Plots a circle between noncanonical base pairs. Input is a 2 column text file where each line contains the nucleotide numbers of a noncanonical base pair. If the provided varna or xrna input already displays these nucleotides as paired, the circle is overlaid on to the base pairing line.  
`--switch` - reverse the pairing coloring scheme used by --ct.  
`-f,--ntFontSize FontSizeNum` - Sets the font size for all nucleotides. The default font size is 16.  
`--hideNtNums` - Don't display nucleotide numbering.  

### Output Description
The script will generate a visualization of the secondary structure, with annotations, in the SVG format. SVG files can be open with dedicated SVG viewers like Gapplin or Adobe SVG viewer. Most viewers will allow the user to export svg files in png or pdf formats.  
SVG files can be edited with Adobe Illustrator.

# contactDistanceHistogram.py  
Maps internucleotide pairs onto provided secondary structure in ct format. Returns a histogram of the contact distances between these pairs.  
The contact distance histrogram is plotted against a histogram of all pairwise distances in the given secondary structure.  

*Contact Distance* is a measurement of the shortest path between two nucleotides where the RNA secondary structure is represented as a graph with nucleotides as vertices and with base pair and backbone connections as edges. The contact distance is calculated using a "breadth first search" algorithm, implemented in RNAtools2.  

### Requirements
This script requires the python 2 modules matplotlib and numpy.  
The python script RNAtools2 is also required, available here: https://github.com/Weeks-UNC/RNATools   

### Required Inputs
**Deletion File** - list of nucleotide pairs and their frequency. See *File Description* section at bottom of readme for more detail.  
**Gene** - Input a gene or sequence name to limit nucleotide pairs plotted to just those from a specific sequence name contained in the deletion file.  
**CT File** - Text File in ct format describing secondary structure. See *File Descriptions* at bottom of README for more details.  

### Execution Instructions
`python contactDistanceHistogram.py DeletionFile.txt geneName SecondaryStructure.ct`  
Example: `python contactDistanceHistogram.py RNasePDeletions.txt rnasep RNasePSecondaryStructure.ct`  

### Optional Inputs
`--lengthAdjust` - Create background histogram from random Nt pairs that reflect the 1D length of deletions.   
`--seqLength INT` - Provide length of sequence in ct file to limit input deletion pairs to just those within the ct file. Useful for preventing "out of bounds" errors.  

### Output Description
Generates a pdf with 3 contact distance histograms, one of all contact distances in the deletion file, one of the 97th percentile and one of the 98th. All histograms are plotted over a histogram of all possible contact distances in the secondary structure.  

# deletionLengthDistribution.py  

### Requirements
### Required Inputs
### Execution Instructions
### Optional Inputs
### Output Description

# File Descriptions  
### Deletion File
Example:  
Total Reads Aligned:11450       Total Deletions:728.0  
rnasep  123     184     0.0017623086  
rnasep  81      94      0.0013324450  

The first line is a header, denoting total reads aligned and deletions observed. In most cases a line must be present but the actual contents are unimportant.  
Subsequent lines follow the same 4 column format:  
- Column 1 = Reference name of gene or sequence contact is found in.
- Column 2 = Deletion start site.
- Column 3 = Deletion stop site.
- Column 4 = Deletion rate frequency. This value is often normalized by read depth.

### Contact File
Example:  
nt1 nt2 rate  
123 184 0.0017623086  
81  94  0.0013324450  

The first line is a header. The header line must be present but the contents are unimportant.  
- Column 1 = Contact start site.  
- Column 2 = Contact stop site.  
- Column 3 = Contact frequency.  

**CT File**
Format to describe RNA secondary structure.
See https://rna.urmc.rochester.edu/Text/File_Formats.html for detailed description.

# StructureAnalysis
Collection of Scripts to Analyze 2º and 3º RNA Structure with a focus on internucleotide measurements.  
Scripts were designed to work with outputs of ShapeJumper https://github.com/Weeks-UNC/ShapeJumper_V1  
However, they can be applied to work with other internucleotide pair data sets.  

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
`--inclusivePercentile` - By default the percentile cutoff is calculated from only the set of nucleotide pairs present in the given PDB file. This option changes the percentile cutoff calculation to include all nucleotide pairs. This option will only affect the output if nucleotides are missing from the PDB. `--topDels INT` - Limits pairs plotted to given number rather than percentile eg. most frequent 100 pairs.    
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

# draw_3D_RNA_contacts.py  

### Requirements  
### Required Inputs  
### Execution Instructions  
### Optional Inputs  
### Output Description 

# AnnotateSecondaryStructure.py  

### Requirements  
### Required Inputs  
### Execution Instructions  
### Optional Inputs  
### Output Description  

# contactDistanceHistogram.py  

### Requirements  
### Required Inputs  
### Execution Instructions  
### Optional Inputs  
### Output Description  

# deletionLengthDistribution.py  

### Requirements  
### Required Inputs  
### Execution Instructions  
### Optional Inputs  
### Output Description  

# File Descriptions  
**Deletion File**  
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

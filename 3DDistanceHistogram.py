#3DDistanceHistogram.py
#DEPENDENCIES: RNAtools2 package from https://github.com/Weeks-UNC/RNATools
#Created By Tom Christy
#Major Edit, January 16, 2019
#Changed this script to only print one histogram that is larger and more legible. This should make it more figure friendly.
#You Have to give it a percentile now.

#This script takes in a pdb and a deletions file from ShapeJumper. It then compares
#all possible 3D distances to those in the deletion file. A

import matplotlib
#matplotlib.use("PDF")
from matplotlib import pyplot as plt
import numpy as np
from operator import itemgetter
from Bio.PDB import *
import re
import argparse
import sys
from RNAtools2 import CT
from matplotlib.patches import Rectangle

def parseArgs():
	prs = argparse.ArgumentParser()
	prs.add_argument("pdb",type=str,help ="PDB input file to determine distances from")
	prs.add_argument("deletionFile",type=str,help ="deletions input file. 4 Column File.")
	prs.add_argument("gene",type=str,help ="target gene to analyze from deletion file")
	prs.add_argument("percentile",type=float,help="Percentile cutoff to plot")
	prs.add_argument("-u","--upperPercentile",type=float,default=100.0,help="Apply an upper cutoff as well, so all deletions are below this percentile")
	prs.add_argument('-c',"--color",type=str,default="black",help="input a color for the experimental histogram")
	prs.add_argument("--experimentalOnly",action='store_true',help="Plot only the experimental, no background lengths.")
	prs.add_argument("--inclusivePercentile",action='store_true',help="this option includes deletions not found in the crystal in the calculation of the percentile cutoffs")
	prs.add_argument("--topDels",action="store",type=int,help="Input an integer. The ouput histograms will then the top deletions of this count rather than percentiles. ie 100 will get the top 100 deletions.")
	prs.add_argument("--filledLoop", action='store_true',help='labels histograms that the full PDB structure is being used (for RNaseP-Cat)')
	prs.add_argument("--allOptimal", action='store_true',help='writes out a list of ALL deletions between 10 and 25 Angstroms')
	prs.add_argument("--lengthAdjust", type=int, default=0, help='Create background histogram from random Nt pairs that reflect the 1D length of deletions. Requires length of sequence. Assumes structure cassettes included in length.')
	prs.add_argument("--shapeFile",type=str,default="",help ="A shape file, in .shape format, for each nucleotide")
	prs.add_argument("--shapeCutoff", type = float, default=0.4, help = "minimum SHAPE reactivity to allow nts with, assumes shapeFile has been given")
	prs.add_argument("--oneNTSelection",action='store_true',default=False,help="If true, only one nt has to be reactive in deletion, else both do. Assumes a shape file has been given")
	prs.add_argument("--shift5nt",type=int,default=0,help="Shift all deletion start sites, the 5 prime end, by the input value, either positive or negative")
	prs.add_argument("--shift3nt",type=int,default=0,help="Shift all deletion start sites, the 3 prime end, by the input value, either positive or negative")
	prs.add_argument("--centralPoint",action='store_true',default=False,help="Rather than calculating distances between hydroxyl groups, calculate between central point of the base")
	prs.add_argument("--corrChi",action='store_true',help="Set this option if the input file is chi square correlation file.")
	prs.add_argument("--StrucFilt",action='store',nargs=2,default = [],help="Takes in a CT file and an integer. After Percentile cutoffs, filters deletions by primary and secondary structure distance. Modeled off of DMD filtering procedure.")
	o = prs.parse_args()
	return o
	
def format(value):
    return "%.10f" % value
    
def loadPDB(filename):
	parser = PDBParser()
	structure = parser.get_structure(filename, filename)
	return structure
	
def getResNums(structure):
	# input must be a structure object generated with PDBparser
	# as in the loadPDB function or on its own.
	resnums = []
	residues = structure.get_residues()
	counter = 0
	for i in residues:
		counter += 1
		atoms = []
		for a in i:
			atoms.append(str(a))
		if '<Atom O2\'>' in atoms:
			#print i
			# split residue object string into fields and
			# then extract resnum from resseq=N field. This is hacky and I know it.
			resnums.append(int(str(i).split()[3].split('=')[1]))
	return resnums

def getDistances(pairs, structure, chain):
	# pairs = list of i,j tuples (can be i,j,k,...; only i,j considered)
	# calculates and returns the distance in angstroms
	# between the 2' oxygen atoms of residues i and j.
	# If you want ALL possible pairwise distances, use getALLdistances()
	distances = []
	residues = getResNums(structure)
	for i in pairs:
		j = i[0]
		k = i[1]
		#set to I for 23S
		#set to V for 16S
		if abs(j-k) > 10 and j in residues and k in residues:
			#handle case for central point of base calculation
			if(args.centralPoint):
				Base1Coords = getCentralCoord(structure, chain, j)
				Base2Coords = getCentralCoord(structure, chain, k)
				#prevent residues without base coordinates from causing issues
				if("NO" not in Base1Coords and "NO" not in Base2Coords):
					dist = calcDistance(Base1Coords, Base2Coords)
					distances.append([j,k,i[2],dist])
			else:		
				atom1 = structure[0][chain][j]['O2\'']
				atom2 = structure[0][chain][k]['O2\'']
				dist = atom1 - atom2
				distances.append([j,k,i[2],dist])
#		if(j not in residues or k not in residues):
#			print "Deletion "+ str(j) +" to "+ str(k)+ " not in set of residues"
#			if(j not in residues):
#				print str(j) + " not in set of residues"
#			if(k not in residues):
#				print str(k)+ " not in set of residues"
	return distances

def getALLdistances(structure, chain):
	# calculates ALL pairwise thru-space distances for a given structure.
	all_dists = []
	residues = getResNums(structure)
	for i in residues:
		for j in residues:		
			if abs(i-j) > 10:
				if(args.centralPoint):
					Base1Coords = getCentralCoord(structure, chain, i)
					Base2Coords = getCentralCoord(structure, chain, j)
					if("NO" not in Base1Coords and "NO" not in Base2Coords):
						dist = calcDistance(Base1Coords, Base2Coords)
						all_dists.append([i,j,dist])
				else:
					atom1 = structure[0][chain][i]['O2\'']
					atom2 = structure[0][chain][j]['O2\'']
					dist = atom1 - atom2
					all_dists.append([i,j,dist])
	return all_dists

def getDelDistance1D(deletions):
	#This function takes in a list of deletions and returns a list of all the 1D distances of
	#the deletions
	lengths = []
	for line in deletions:
		delLength = line[1] - line[0]
		lengths.append(delLength)
	return lengths

def createDeletions(delLengths, seqLen, pdbDist):
	#This function takes in a list of deletion lengths and a sequence length. These are used to create fake deletions
	#that span the whole structure but have the same deletions length distribution.
	
	#prep, take in the list of all PDBDistances and turn it into a dictionary.
	#This will allow me to easily check if fake deletions exist in the crystal structure and then 
	#pull their distances.
	pdbDict = {}
	for line in pdbDist:
		pdbDict[(line[0], line[1])] = line[2]
	#determine the total number of deletions to create. This is derived from all possible NT to NT pairs, corrected for 
	#structure cassettes and 10NT del Length limit.
	totalDeletions = (seqLen - 10)**2
	count = 0
	iterated = 0
	deletions = []
	
	while count < totalDeletions:
		iterated += 1
		#pull out a deletion length
		delLength = np.random.choice(delLengths)
		delStop = np.random.randint(seqLen-42)
		delStart = delStop - delLength
		#make sure deletions land within actual sequence ie. not in a Structure cassette
		if (delStart in range(15,seqLen-43) and delStop in range(15,seqLen-43)):
			#now check to make sure the sites are in the crystal
			delStart = delStart - 14
			delStop = delStop - 14
			if(pdbDict.has_key((delStart,delStop))):
				deletions.append([delStart, delStop, pdbDict[(delStart,delStop)]])
				count += 1
	return deletions

def calcDistance(coords1, coords2):
	#This script takes in 2 lists of xyz coordinates and calculates the distance between them
	x1 = coords1[0]
	y1 = coords1[1]
	z1 = coords1[2]
	x2 = coords2[0]
	y2 = coords2[1]
	z2 = coords2[2]
	dist = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
	return dist


def getCentralCoord(structure, chain, residueNum):
	#This function takes in residue number and returns the central xyz coordinates of the base of that residue
	
	#list of atoms to not use
	sp = ["P", "OP1", "OP2", "O5\'", "C5\'", "C4\'", "O4\'", "C3\'", "O3\'", "C2\'", "O2\'", "C1\'"]
	#instantiate coordinates and atomCount
	x = 0.0
	y = 0.0
	z = 0.0
	atomCount = 0
	#get atoms in residue and run through them
	residue = structure[0][chain][residueNum]
	for atom in residue:
		atomName = atom.get_name()
		if(atomName not in sp):
			atomCount += 1
			coords = atom.get_coord()
			x += coords[0]
			y += coords[1]
			z += coords[2]
	
	#now that you've add up all the xyz values, get their average and return it
	if(atomCount != 0):
		x = x/atomCount
		y = y/atomCount
		z = z/atomCount
	else:
		x = "NO"
	return [x, y, z]

def filterByStruc(delList, bps, cutOff):
	#Takes in the list of deletion, list of base pairs and distance cutoff
	#Filters deletions that are below the cutoff in primary or secondary space, not contact distance, 
	#works the same way as DMD's filtering method for RINGS

	new_pairs = []
	for pair in delList:
		add = True
		#pull the nucleotide numbers of the deletion pair
		nt1 = pair[0]
		nt2 = pair[1]
		#curate based on primary
		if abs(nt1 - nt2) <= cutOff:
			add = False
		#curate based on Secondary
		for bp in bps:
			if abs(nt1 - bp[0]) <= cutOff:
				if abs(nt1 - bp[0]) + abs(nt2 - bp[1]) <= cutOff:
					add = False
		if add:
			new_pairs.append(pair)
	delList = new_pairs
	return delList

args = parseArgs()
pdb = loadPDB(args.pdb)
#determine chain of molecule
chain = ''
for model in pdb:
	for c in model:
		chain = str(c)

p = re.compile('Chain id=(\w)')
m = p.search(chain[1:-1])
if m:
	chain = m.group(1)

print chain

pdb_dist = getALLdistances(pdb, chain)

#find max distance
max = 0
for line in pdb_dist:
	if(line[2]>max):
		max = line[2]
max = int(max)
#if selected, load in a .shape files reactivities into a dict
if(args.shapeFile != ""):
	shapeReactivities = {}
	inShapeFile = open(args.shapeFile, 'r')
	shapeLines = inShapeFile.readlines()
	inShapeFile.close()

	for line in shapeLines:
		cols = line.split()
		shapeReactivities[int(cols[0])] = float(cols[1])

#load in deletion file
inF = open(args.deletionFile,'r')
inLines = inF.readlines()
inF.close()
if(args.corrChi):
	#remove correlation file header
	inLines.pop(0)
targetGene = args.gene
#load deletions in multi-dimensional array
deletions = list()
for i in range(1,len(inLines)):
	line = inLines[i].strip().split()
	if(line[0] == targetGene or args.corrChi):
		#if SHAPE Filtering has been selected do it here
		if(args.corrChi):
			nt1 = int(line[0])
			nt2 = int(line[1])
		else:
			nt1 = int(line[1])
			nt2 = int(line[2])
		#shift nts if selected, if not the default is 0 so they won't change
		nt1 = nt1 + args.shift5nt
		nt2 = nt2 + args.shift3nt
		if(args.shapeFile != ""):
			aboveSHAPE = False		
			#make sure nucleotides exist within range of nucleotides
			if(shapeReactivities.has_key(nt1) and shapeReactivities.has_key(nt2)):
				if(args.oneNTSelection):
					if(shapeReactivities[nt1] > args.shapeCutoff or shapeReactivities[nt2] > args.shapeCutoff):
						aboveSHAPE = True
				#if not selected, make sure both nts have shape reactivities above the cutoff
				else:
					if(shapeReactivities[nt1] > args.shapeCutoff and shapeReactivities[nt2] > args.shapeCutoff):
						aboveSHAPE = True
				#if the selected number of NTS are above the cutoff, put deletion into list, otherwise remove it from consideration
				if(aboveSHAPE):
					newLine = [0,0,0.0]
					newLine[0] = nt1
					newLine[1] = nt2 
					newLine[2] = float(line[3])
					deletions.append(newLine)
		else:
			newLine = [0,0,0.0]
			newLine[0] = nt1
			newLine[1] = nt2
			if args.corrChi:
				newLine[2] = float(line[4])
			else:
				newLine[2] = float(line[3])
			deletions.append(newLine)
#make sure that data has been loaded into deletions (ie. make sure there are actually deletions)
try:
	x = deletions[0]
except IndexError:
	sys.exit("There are no deletions in this file!")

rxn_dist = getDistances(deletions, pdb, chain)

#if asked for, write out all deletions that span between 10 and 25 angstroms
if(args.allOptimal):
	outF = open(args.deletionFile[:-4]+'_optimalDels.txt','w')
	outF.write('Optimal 10 - 25 Angstrom Deletions\n')
	for delLine in rxn_dist:
		if(delLine[3] > 10 and delLine[3] < 25):
			outF.write(str(delLine[0])+'\t'+str(delLine[1])+'\t'+str(format(delLine[2]))+'\n')

# Calculate median and stdev of deletion frequencies, use as master cutoff for all data.
median = np.median([i[2] for i in rxn_dist])
stdev = np.std([i[2] for i in rxn_dist])

if args.topDels:
	percA = int(args.topDels[0])
else:
	upperPerc = args.upperPercentile
	percA = args.percentile
	

	
if(args.inclusivePercentile):
	if args.topDels:
		#sort the deletion rates and pull the values of top A and B so we can pull out all deletions with a
		#rate higher than A or B
		deletions = sorted(deletions, key=lambda x: x[2], reverse=True)
		pA = deletions[percA][2]
	else:
		upperP = np.percentile([i[2] for i in deletions],upperPerc)
		pA = np.percentile([i[2] for i in deletions],percA)
		
else:
	if args.topDels:
		#sort the deletion rates and pull the values of top A and B so we can pull out all deletions with a
		#rate higher than A or B
		rxn_dist = sorted(rxn_dist, key=lambda x: x[2], reverse=True)
		pA = rxn_dist[percA][2]
	else:
		upperP = np.percentile([i[2] for i in rxn_dist],upperPerc) 
		pA = np.percentile([i[2] for i in rxn_dist],percA)

rxn_dist_pA = [i for i in rxn_dist if (i[2] > pA and i[2]<= upperP)]

print "percentile "+str(percA)
print pA

#if requested, remove any of deletions that don't pass the DMD type secondary structure filtering 
if(len(args.StrucFilt) == 2 ):
	if(isinstance(args.StrucFilt[0],str) and isinstance(int(args.StrucFilt[1]),int)):
		rnaCT = CT(args.StrucFilt[0])
		helices = rnaCT.extractHelices()
		basePairs = []
		#extract base pairs for helices and load them as tuples into a list
		for key in helices.keys():
			helix = helices[key]
			for bp in helix:
				basePairs.append(bp)
			
		rxn_dist_pA = filterByStruc(rxn_dist_pA, basePairs, int(args.StrucFilt[1]))
		rxn_dist = filterByStruc(rxn_dist, basePairs, int(args.StrucFilt[1]))
	else:
		print "ERROR! If you want to filter by secondary structure you need to input a ct file and then an integer!"
		sys.exit()

#print "del 99 is"
#print delp99
#print "rxn 99 is "
#print rxn_dist_p99

###testing percentiles!!!!!!####
#print "median = "+str(median)


#Steve's Idea: Plot deletions against a set of random deletions that mimic the 1D distance of 
#the deletions in each set. 

#get 1D distances of all deletion sets
#using these deletions lengths, create a list of 3D distances
if(args.lengthAdjust != 0):
	rxn_dist_pA_1DBackground = createDeletions(getDelDistance1D(rxn_dist_pA), args.lengthAdjust, pdb_dist)

#create title for saved figure. This will also be the plot's title
if args.corrChi:
	title = args.deletionFile[:-5]
else:
	title = args.deletionFile[:-4]
if args.filledLoop:
	title += "_FilledLoop"
if args.lengthAdjust != 0:
	title += "_lengthAdjustedBackground"
if args.shapeFile != "":
	title += "_SHAPE"
	title += "_"+str(args.shapeCutoff)[2:]+"Cutoff"
	if(args.oneNTSelection):
		title += "_1ntFilter"
	else:
		title += "_2ntFilter"
if args.shift5nt != 0:
	title += "_shift5nt"+str(args.shift5nt)
if args.shift3nt != 0:
	title += "_shift3nt"+str(args.shift3nt)
if(args.inclusivePercentile):
	title += "_inclusivePercentile"
if(args.centralPoint):
	title += "_BaseToBase"
if len(args.StrucFilt) == 2:
	title += "_StrucFilt"+str(args.StrucFilt[1])
	
if args.topDels:
	title += "_3DDistancehist_Top"+str(percA)+"Dels.pdf"
elif args.experimentalOnly:
	title +="_3DDistancehist_"+str(percA)+"PercNoBG.pdf"
else:
	if upperPerc == 100.0:
		title += "_3DDistancehist_"+str(percA)+"Perc.pdf"
	else:
		title += "_3DDistancehist_"+str(percA)+"_to_"+str(upperPerc)+"Perc.pdf"

#change the font size
plt.rcParams.update({'font.size': 10})
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['axes.linewidth'] = 1.5

### Plot the distance disributions ###
fig = plt.figure("Distance Distributions", figsize=(8,4))
ax1 = plt.gca()
ax1.tick_params(axis='both',direction='out',length=4,width=1.5)

#bins = np.linspace(0,max,num=25)
bins = np.arange(0,max+5,5)
#set up bins for deletions, same spacing just limited by the highest deletion distance
maxDel = np.max(rxn_dist_pA, axis=0)[3]
minDel = np.min(rxn_dist_pA, axis=0)[3]
expBins = np.arange(int(minDel) / 5 * 5, maxDel + 5, 5)

#Could uncomment this to draw sample background of all deletions#
#if(args.lengthAdjust != 0):
#	ax1.hist([i[2] for i in rxn_dist_1DBackground], bins=bins, color='purple', alpha=0.2)
#	plt.text(0.8, 0.81,' Background\nMean Distance:\n'+str(round(np.mean([i[2] for i in rxn_dist_1DBackground]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
#else:
########Plot Percentile A Cutoff Deviation Hist ########################

print "Experimental mean 3d distance "+str(round(np.mean([i[3] for i in rxn_dist_pA]),2))

if args.experimentalOnly == False:
	if (args.lengthAdjust != 0):
		hist1, b1, p1 = ax1.hist([i[2] for i in rxn_dist_pA_1DBackground], bins=bins, label="Pairwise", color='purple',
								 alpha=0.2)
		
		#plt.text(0.77, 0.72,' Pairwise\nMean Distance:\n'+str(round(np.mean([i[2] for i in rxn_dist_pA_1DBackground]),2))+ " "+ u'\u00C5', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	else:
		ax1.hist([i[2] for i in pdb_dist], bins=bins, color='purple', label="Pairwise", alpha=0.2)
		#plt.text(0.77, 0.72,' Crystal\nMean Distance:\n'+str(round(np.mean([i[2] for i in pdb_dist]),2))+ " "+ u'\u00C5', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	#plt.text(0.72, 0.95,'Total Deletions:\n'+str(len(deletions)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	
	ax1.set_ylabel('Count (reference)',fontsize=12)
	plt.xlabel("Distance (" + u'\u00C5' + ")\n",fontsize=12)
	lines, labels = ax1.get_legend_handles_labels()
	plt.ylim(0,5200)

	#plt.xlim(0,120)
	#ax1.set_xbound(0,120)
	ax2 = ax1.twinx()
	ax1.yaxis.tick_right()
	ax1.yaxis.set_label_position("right")
	hist2 = ax2.hist([i[3] for i in rxn_dist_pA], bins=expBins, histtype='step', color=args.color, linewidth=2, label='Experimental')
	#ax2.set_xbound(0,120)
	ax2.yaxis.tick_left()
	ax2.yaxis.set_label_position("left")
	
	#plt.ylabel('Count (Experimental)')

	##add text labels on mean distances
	#if args.topDels:
	#	plt.title('Experimental vs. All Pairwise\n3D Distances (Top ' +str(percA)+' Dels)')
	#	plt.text(0.8, 0.94,'Top '+str(percA)+' Deletions:\n'+str(len(rxn_dist_pA)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	#else:
		#plt.title('Experimental vs. All Pairwise\n3D Distances ('+str(percA)+' Perc)')
	#plt.text(0.8, 0.85,str(percA)+' Perc Deletions:\n'+str(len(rxn_dist_pA)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	#plt.text(0.77, 0.85,'Deletions:\n'+str(len(rxn_dist_pA)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	
	#determine max count
	bottom, top = plt.ylim()
	print "Max count is near: "+str(top)
	
	plt.xlim(1,max)
	#99th percentile
	#plt.ylim(0,7.5)
	#97th percentile
	#plt.ylim(0,12.5)
	
	#plt.ylim(0,26)
	
	#plt.text(0.72, 0.94,'Mean Distance:\n'+str(round(np.mean([i[3] for i in rxn_dist_pA]),2)) + " "+ u'\u00C5', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	ax2.set_xticks(range(0,126,25))
	plt.xlim(0,115)
	lines2, labels2 = ax2.get_legend_handles_labels()
	
	plt.ylabel('Count (experimental)',fontsize=12)
	ax2.tick_params(axis='both',direction='out',length=5,width=1.5)
	ax2.legend(lines + lines2, labels + labels2, loc=0)
	#ax1.legend(loc=0)
else:
	hist2 = ax1.hist([i[3] for i in rxn_dist_pA], bins=expBins, histtype='step', color=args.color, linewidth=2, label='Experimental')
	#plt.ylim(0,26)
	plt.xlabel("Distance (" + u'\u00C5' + ")",fontsize=12)
	plt.ylabel('Count (Experimental)',fontsize=12)
	ax1.set_xticks(range(0,126,25))
	plt.xlim(0,115)
	plt.legend()

if(args.lengthAdjust == 0):
	plt.tight_layout()
plt.tight_layout()
print title
#plt.show()
plt.savefig(title)
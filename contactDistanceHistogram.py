#contactDistanceHistogram.py

#Created by Tom Christy
#November 10, 2016
#DEPENDENCIES: RNAtools package from https://github.com/Weeks-UNC/RNATools

#This Script takes in a set of deletions and a CT File. It plots how the contact distance 
#of the deletions compares to all possible contact distances. It includes an option to adjust the background
#plotted distances to reflect those possible with the same distribution of deletion lengths. Also this script
#plots 4 cutoffs based on percentile of deletion rate.

#import stuffs
import matplotlib
matplotlib.use("PDF")
from RNAtools2 import CT
import argparse
from matplotlib import pyplot as plt
import numpy as np
from operator import itemgetter
import re

def parseArgs():
	prs = argparse.ArgumentParser()
	prs.add_argument("deletionFile",type=str,help="Normalized deletions input file. 4 Columns.")
	prs.add_argument("targetGene",type=str,help = "name of gene to analyze in deletions file")
	prs.add_argument("CTFile",type=str,help = "Reference CT File to assess contact distance with")
	prs.add_argument("--lengthAdjust",action="store_true",help="Background histogram now samples from contact distances with the same 1D length")
	prs.add_argument("--seqLength",type=int, default=0,help="this option takes in the sequence length (without Structure Cassettes) and uses that to remove bad deletions")
	o = prs.parse_args()
	return o

def getDelDistance1D(deletions):
	#This function takes in a list of deletions and returns a list of all the 1D distances of
	#the deletions
	lengths = []
	for line in deletions:
		delLength = line[1] - line[0]
		lengths.append(delLength)
	return lengths

def createDeletions(delLengths, seqLen, ctDist):
	#This function takes in a list of deletion lengths and a sequence length. These are used to create fake deletions
	#that span the whole structure but have the same deletions length distribution.
	
	#prep, take in the list of all CT file derived contact distances and turn it into a dictionary.
	#This will allow me to easily check if fake deletions exist in the ct file and then 
	#pull their distances.
	ctDict = {}
	for line in ctDist:
		ctDict[(line[0], line[1])] = line[2]
	#determine the total number of deletions to create. This is derived from all possible NT to NT pairs, corrected for 
	#10NT del Length limit.
	totalDeletions = (seqLen - 10)**2
	count = 0
	iterated = 0
	deletions = []
	
	while count < totalDeletions:
		iterated += 1
		#pull out a deletion length
		delLength = np.random.choice(delLengths)
		delStop = np.random.randint(seqLen+1)
		delStart = delStop - delLength
		#make sure deletions land within actual sequence ie. not in a Structure cassette
		if (delStart in range(1,seqLen+1) and delStop in range(1,seqLen+1)):
			#now check to make sure the sites are in the crystal
			if(ctDict.has_key((delStart,delStop))):
				deletions.append([delStart, delStop, ctDict[(delStart,delStop)]])
				count += 1
	return deletions

args = parseArgs()
#load deletion file
inF1 = open(args.deletionFile,'r')
delLines = inF1.readlines()
inF1.close()

#load deletions into a 2D array
targetGene = args.targetGene
#load deletions in multi-dimensional array
deletions = list()
for i in range(1,len(delLines)):
	line = delLines[i].strip().split()
	if(line[0] == targetGene):
		newLine = [0,0,0.0]
		newLine[0] = int(line[1])
		newLine[1] = int(line[2]) 
		newLine[2] = float(line[3])
		if(args.seqLength != 0):
			if(newLine[0] < args.seqLength and newLine[1] < args.seqLength):
				deletions.append(newLine)
		else:
			deletions.append(newLine)

#load CT File
rnaCT = CT(args.CTFile)
rnaLen = len(rnaCT.num)
print "Length of RNA: "+str(rnaLen)
#determine all possible contact distances
totalContactDistance = list()
for i in range(1,len(rnaCT.num)+1):
	for j in range(i,len(rnaCT.num)+1):
		totalContactDistance.append([i,j,int(rnaCT.contactDistance(i,j))])

#determine samples contact distances at different percentile cutoffs
median = np.median([i[2] for i in deletions])
stdev = np.std([i[2] for i in deletions])
p97 = np.percentile([i[2] for i in deletions],97)
p98 = np.percentile([i[2] for i in deletions],98)
p99 = np.percentile([i[2] for i in deletions],99)
p99_5 = np.percentile([i[2] for i in deletions],99.5)

sampleContactDistance = list()
ContactDistance_p97 = list()
ContactDistance_p98 = list()
ContactDistance_p99 = list()
ContactDistance_p99_5 = list()
for line in deletions:
	if(line[0] <= rnaLen and line[1] <= rnaLen):
		CD = rnaCT.contactDistance(line[0],line[1])
		sampleContactDistance.append([line[0],line[1],CD])
		if(line[2] > p97):
			ContactDistance_p97.append([line[0],line[1],CD])
			if(line[2] > p98):
				ContactDistance_p98.append([line[0],line[1],CD])
				if(line[2] > p99):
					ContactDistance_p99.append([line[0],line[1],CD])
					if(line[2] > p99_5):
						ContactDistance_p99_5.append([line[0],line[1],CD])



#if requested create length adjusted background contact distance distributions for each percentile
if(args.lengthAdjust == True):
	#rxn_dist_1DBackground = createDeletions(getDelDistance1D(rxn_dist), args.lengthAdjust, pdb_dist)
	CD_p97_1DBackground = createDeletions(getDelDistance1D(ContactDistance_p97), len(rnaCT.num), totalContactDistance)
	CD_p98_1DBackground = createDeletions(getDelDistance1D(ContactDistance_p98), len(rnaCT.num), totalContactDistance)
	CD_p99_1DBackground = createDeletions(getDelDistance1D(ContactDistance_p99), len(rnaCT.num), totalContactDistance)
	CD_p99_5_1DBackground = createDeletions(getDelDistance1D(ContactDistance_p99_5), len(rnaCT.num), totalContactDistance)

# determine max contact distance
max = 0
for line in totalContactDistance:
	if(line[2]>max):
		max = line[2]

if(args.lengthAdjust):
	titleHigh = args.deletionFile[:-4]+'_LengthAdjusted_ContactDistanceHist_99and99_5Perc.pdf'
	titleLow = args.deletionFile[:-4]+'_LengthAdjusted_ContactDistanceHist_97and98Perc.pdf'
else:
	titleHigh = args.deletionFile[:-4]+'_ContactDistanceHist_99and99_5Perc.pdf'
	titleLow = args.deletionFile[:-4]+'_ContactDistanceHist_97and98Perc.pdf'

### Plot the distance disributions ###
fig1 = plt.figure("High Percentile Contact Distance Distributions", figsize=(20,10))
ax1 = plt.subplot(131)
#bins = np.linspace(0,max,num=25)
bins = np.arange(0,max+5,5)
ax1.hist([i[2] for i in totalContactDistance], bins=bins, color='blue', alpha=0.2)
plt.text(0.8, 0.81,' Crystal\nMean Distance:\n'+str(round(np.mean([i[2] for i in totalContactDistance]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))

plt.ylabel('Count (Pairwise)')
plt.xlabel('Contact Distance')
ax2 = ax1.twinx()
ax2.hist([i[2] for i in sampleContactDistance], bins=bins, histtype='step', color='red', linewidth=2)
plt.title('Experimental vs. All Pairwise\nContact Distances (no cutoff)')
plt.text(0.8, 0.95,'Total Deletions:\n'+str(len(deletions)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean([i[2] for i in sampleContactDistance]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))

ax1 = plt.subplot(132)
if(args.lengthAdjust):
	ax1.hist([i[2] for i in CD_p99_1DBackground], bins=bins, color='purple', alpha=0.2)
	plt.text(0.8, 0.81,' Background\nMean Distance:\n'+str(round(np.mean([i[2] for i in CD_p99_1DBackground]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
else:
	ax1.hist([i[2] for i in totalContactDistance], bins=bins, color='blue', alpha=0.2)
plt.xlabel(titleHigh[:-18])
ax2 = ax1.twinx()
ax2.hist([i[2] for i in ContactDistance_p99], bins=bins, histtype='step', color='orange', linewidth=2, label='99th Percentile')
plt.text(0.8, 0.95,'99 Perc Deletions:\n'+str(len(ContactDistance_p99)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean([i[2] for i in ContactDistance_p99]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
plt.ylabel('Count (Experimental)')
plt.title('Experimental vs. All Pairwise\nContact Distances (99 Perc)')
plt.xlim(1,max)

ax1 = plt.subplot(133)
if(args.lengthAdjust):
	ax1.hist([i[2] for i in CD_p99_5_1DBackground], bins=bins, color='purple', alpha=0.2)
	plt.text(0.8, 0.81,' Background\nMean Distance:\n'+str(round(np.mean([i[2] for i in CD_p99_5_1DBackground]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
else:
	ax1.hist([i[2] for i in totalContactDistance], bins=bins, color='blue', alpha=0.2)
plt.xlabel('Contact Distance')
ax2 = ax1.twinx()
ax2.hist([i[2] for i in ContactDistance_p99_5], bins=bins, histtype='step', color='green', linewidth=2, label='99.5th Percentile')
plt.ylabel('Count (Experimental)')
plt.title('Experimental vs. All Pairwise\nContact Distances (99.5 Perc)')
plt.xlim(1,max)
plt.text(0.8, 0.95,'99.5 Perc Deletions:\n'+str(len(ContactDistance_p99_5)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean([i[2] for i in ContactDistance_p99_5]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))


plt.tight_layout()
#plt.show()
plt.savefig(titleHigh)
print titleHigh

###################################################################
### Plot the distance disributions of lower Percentile cutoffs ###
fig2 = plt.figure("Low Percentile Contact Distance Distributions", figsize=(20,10))
ax1 = plt.subplot(131)
#bins = np.linspace(0,max,num=25)
bins = np.arange(0,max+5,5)
ax1.hist([i[2] for i in totalContactDistance], bins=bins, color='blue', alpha=0.2)
plt.text(0.8, 0.81,' Crystal\nMean Distance:\n'+str(round(np.mean([i[2] for i in totalContactDistance]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))

plt.ylabel('Count (Pairwise)')
plt.xlabel('Contact Distance')
ax2 = ax1.twinx()
ax2.hist([i[2] for i in sampleContactDistance], bins=bins, histtype='step', color='red', linewidth=2)
plt.title('Experimental vs. All Pairwise\nContact Distances (no cutoff)')
plt.text(0.8, 0.95,'Total Deletions:\n'+str(len(deletions)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean([i[2] for i in sampleContactDistance]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))

ax1 = plt.subplot(132)
if(args.lengthAdjust):
	ax1.hist([i[2] for i in CD_p97_1DBackground], bins=bins, color='purple', alpha=0.2)
	plt.text(0.8, 0.81,' Background\nMean Distance:\n'+str(round(np.mean([i[2] for i in CD_p97_1DBackground]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
else:
	ax1.hist([i[2] for i in totalContactDistance], bins=bins, color='blue', alpha=0.2)
plt.xlabel(titleLow[:-16])
ax2 = ax1.twinx()
ax2.hist([i[2] for i in ContactDistance_p97], bins=bins, histtype='step', color='orange', linewidth=2, label='97th Percentile')
plt.ylabel('Count (Experimental)')
plt.title('Experimental vs. All Pairwise\nContact Distances (97 Perc)')
plt.xlim(1,max)
plt.text(0.8, 0.95,'97 Perc Deletions:\n'+str(len(ContactDistance_p97)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean([i[2] for i in ContactDistance_p97]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))


ax1 = plt.subplot(133)
if(args.lengthAdjust):
	ax1.hist([i[2] for i in CD_p98_1DBackground], bins=bins, color='purple', alpha=0.2)
	plt.text(0.8, 0.81,' Background\nMean Distance:\n'+str(round(np.mean([i[2] for i in CD_p98_1DBackground]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
else:
	ax1.hist([i[2] for i in totalContactDistance], bins=bins, color='blue', alpha=0.2)
plt.xlabel('Contact Distance')
ax2 = ax1.twinx()
ax2.hist([i[2] for i in ContactDistance_p98], bins=bins, histtype='step', color='green', linewidth=2, label='98th Percentile')
plt.ylabel('Count (Experimental)')
plt.title('Experimental vs. All Pairwise\nContact Distances (98 Perc)')
plt.xlim(1,max)
plt.text(0.8, 0.95,'98 Perc Deletions:\n'+str(len(ContactDistance_p98)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean([i[2] for i in ContactDistance_p98]),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))

plt.tight_layout()
#plt.show()
plt.savefig(titleLow)
print titleLow
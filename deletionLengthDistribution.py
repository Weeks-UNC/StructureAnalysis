#deletionLengthDistribution.py

#Created by Tom Christy
#November 15, 2016

#This script takes in a list of nucleotide pairs, calculates the distance between them and plots a distribution.
#The distribution of the 99.5, 99 and 98th percentile by rate frequency are also plotted.
import argparse
import matplotlib
matplotlib.use("PDF")
import math
import matplotlib.pyplot as plt
import numpy as np
def parseArgs():
	prs = argparse.ArgumentParser()
	prs.add_argument("input", type=str, help="input a 4 column tab delimeted deletions.txt file. col1=geneName col2=nt1 col3=nt2, col4=rate")
	prs.add_argument('-w','--weighted',action='store_true',default=False,help='factor deletion count into output')
	prs.add_argument("--custom",type=int,nargs=3, default = [], action="store", help="if you select this option, input 3 custom percentile cutoffs to use rather than the default")
	prs.add_argument("-n","--noBackground",action='store_true',help="if selected, purple background of all deletions will not display")
	prs.add_argument("-l","--length",type=int,default=0,help ="Length of the sequence, without structure casettes. Used to define longest possible deletion")
	o = prs.parse_args()
	return o

def getMax1DDistance(deletions):
	#This function takes in a list of deletions and returns a list of all the 1D distances of
	#the deletions
	lengths = []
	for line in deletions:
		delLength = abs(int(line[1]) - int(line[0]))
		lengths.append(delLength)
	return max(lengths)


args = parseArgs()
#read in deletions file
inF = open(args.input, 'r')
inLines = inF.readlines()
inF.close
#load deletions into a matrix
inputDels = np.zeros((len(inLines)-1,3))
for i in range(1, len(inLines)):
	cols = inLines[i].strip().split()
	start = int(cols[1])
	stop = int(cols[2])
	delCount = float(cols[3])
	inputDels[i-1]=np.array([start,stop,delCount])

#if the weighted option is selected, first determine the total sum of all deletions' frequencies
#then use that to convert all deletion frequencies into percents of the total

if(args.weighted):
	print inputDels[:,2].sum()
	print np.absolute(inputDels[:,2]).sum()
	print inputDels[:,2].sum()
	tot = np.absolute(inputDels[:,2]).sum()
	inputDels[:,2] = inputDels[:,2] / tot * 100


if(args.length == 0):
	seqLen = getMax1DDistance(inputDels)
else:
	seqLen = args.length

#load an array with the distribution of all possible del lengths >10 nts
allPosDels = np.zeros(seqLen+1)
for i in range(1,seqLen-10):
	for j in range(i,seqLen-10):
		dist = j-i
		if(dist > 10):
			allPosDels[dist] = allPosDels[dist] + 1

#instantiate array of zeros that will hold all possible deletion lengths
dels = np.zeros(seqLen+1)
delsC = np.zeros(seqLen+1)
delsB = np.zeros(seqLen+1)
delsA = np.zeros(seqLen+1)

if args.weighted == False:
	dels = dels.astype(int)
	delsA = delsA.astype(int)
	delsB = delsB.astype(int)
	delsC = delsC.astype(int)

if(len(args.custom) == 3):
	C = args.custom[0]
	B = args.custom[1]
	A = args.custom[2]
else:
	C = 98
	B = 99
	A = 99.5


#set up percentile cutoffs
pC = np.percentile([i[2] for i in inputDels],C)
pB = np.percentile([i[2] for i in inputDels],B)
pA = np.percentile([i[2] for i in inputDels],A)

#determine the lowest deletion frequency, this is useful for determining means on the weighted data
minDelFreq = inputDels[-1][2]
#now handle the subtracted case where the minimum del frequency is negative. In this case we want to get the 
#deletion frequency closest to 0.
if(minDelFreq < 0):
	for i in range(len(inputDels)):
		if(inputDels[i][2]<0):
			minDelFreq = inputDels[i-1][2]
			break
#if the weighted argument is given
#run through deletions
#this will hold each deletion length with repeats for the sake of mean and std
cumDels = list()
#these lists will hold the lengths in each percentile so I can pull a mean
delLensC = list()
delLensB = list()
delLensA = list()
for line in inputDels:
	start = line[0]
	stop = line[1]
	delCount = line[2]
	delLen = abs(int(stop - start))
	#make sure the deletion actually occurs in the sequence
	if(delLen <= seqLen):
		#adjust del length frequency for how often this particular deletion occurs
		if(args.weighted):
			dels[delLen] = dels[delLen] + delCount
			if(delCount > pC):
				delsC[delLen] = delsC[delLen] + delCount
				for i in range(int(delCount/minDelFreq)):
					delLensC.append(delLen)
			if(delCount > pB):
				delsB[delLen] = delsB[delLen] + delCount
				for i in range(int(delCount/minDelFreq)):
					delLensB.append(delLen)
			if(delCount > pA):
				delsA[delLen] = delsA[delLen] + delCount
				for i in range(int(delCount/minDelFreq)):
					delLensA.append(delLen)
			for j in range(int(delCount/minDelFreq)):
				cumDels.append(delLen)
		#handle unadjusted case
		else:
			dels[delLen] = dels[delLen] + 1
			if(delCount > pC):
				delsC[delLen] = delsC[delLen] + 1
				delLensC.append(delLen)
			if(delCount > pB):
				delsB[delLen] = delsB[delLen] + 1
				delLensB.append(delLen)
			if(delCount > pA):
				delsA[delLen] = delsA[delLen] +1
				delLensA.append(delLen)
			cumDels.append(delLen)
#create a custom title based on the options selected
if(args.weighted):
	title = str(args.input)[:-4]+"_Weighted"
else:
	title = str(args.input)[:-4]+"_Unweighted"
if(len(args.custom) == 3):
	title += "_Perc"+str(C)+"_"+str(B)+"_"+str(A)
if args.length != 0:
	title+="_MaxDelLenPossible"+str(seqLen)
#plot deletion frequencies
fig = plt.figure(title, figsize=(20,10))
matplotlib.rcParams.update({'font.size': 18})
if args.noBackground:
	#plot no cutoff
	ax1 = plt.subplot(221)
	ax1.bar(range(seqLen+1),dels, edgecolor = 'r', fill = False, linewidth=0.5)
	plt.xlabel('Sequence Deletion Lengths')
	if(args.weighted):
		[ymin, ymax] = plt.ylim()
		plt.ylim(ymin=0)
		plt.ylabel('Experimental Deletion Frequency')
	else:
		plt.ylabel('Experimental Deletion Count')
	plt.title('Experimental vs. All\nPairwise Distances (no cutoff)')
	plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean(cumDels),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	
	#plot Cth Percentile cutoff
	ax1 = plt.subplot(222)
	ax1.bar(range(seqLen+1),delsC, edgecolor = 'r', fill = False, linewidth=0.5)
	plt.xlabel('Sequence Deletion Lengths')
	if(args.weighted):
		[ymin, ymax] = plt.ylim()
		plt.ylim(ymin=0)
		plt.ylabel(str(C)+'th Percentile Deletion Frequency')
	else:
		plt.ylabel(str(C)+'th Percentile Deletion Count')
		[ymin, ymax] = plt.ylim()
		if ymax < 5:
			ticks = range(int(ymin),int(ymax)+1)
			plt.yticks(ticks)
	plt.title('All Experimental vs. '+str(C)+'th Percentile Deletions')
	plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean(delLensC),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	
	#plot Bth Percentile cutoff
	ax1 = plt.subplot(223)
	ax1.bar(range(seqLen+1),delsB, edgecolor = 'r', fill = False, linewidth=0.5)
	plt.xlabel('Sequence Deletion Lengths')
	if(args.weighted):
		[ymin, ymax] = plt.ylim()
		plt.ylim(ymin=0)
		plt.ylabel(str(B)+'th Percentile Deletion Frequency')
	else:
		plt.ylabel(str(B)+'th Percentile Deletion Count')
		[ymin, ymax] = plt.ylim()
		if ymax < 5:
			ticks = range(int(ymin),int(ymax)+1)
			plt.yticks(ticks)
	plt.title('All Experimental vs. '+str(B)+'th Percentile Deletions')
	plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean(delLensB),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))
	
	#plot Ath Percentile cutoff
	ax1 = plt.subplot(224)
	ax1.bar(range(seqLen+1),delsA, edgecolor = 'r', fill = False, linewidth=0.5)
	plt.xlabel('Sequence Deletion Lengths')
	if(args.weighted):
		[ymin, ymax] = plt.ylim()
		plt.ylim(ymin=0)
		plt.ylabel(str(A)+'th Percentile Deletion Frequency')
	else:
		plt.ylabel(str(A)+'th Percentile Deletion Count')
		[ymin, ymax] = plt.ylim()
		if ymax < 5:
			ticks = range(int(ymin),int(ymax)+1)
			plt.yticks(ticks)
	plt.title('All Experimental vs. '+str(A)+' Percentile Deletions')
	plt.text(0.8, 0.88,'Mean Distance:\n'+str(round(np.mean(delLensA),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))

else:
	#plot no cutoff
	ax1 = plt.subplot(221)
	ax1.bar(range(seqLen+1), allPosDels, alpha =0.2, color = 'grey')
	plt.xlabel('Sequence Deletion Lengths')
	plt.ylabel('All Possible Deletion Frequency')
	ax2 = ax1.twinx()
	ax2.bar(range(seqLen+1),dels, edgecolor = 'r', fill = False, linewidth=0.5)
	if(args.weighted):
		[ymin, ymax] = plt.ylim()
		plt.ylim(ymin=0)
		plt.ylabel('Experimental Deletion Frequency')
	else:
		plt.ylabel('Experimental Deletion Count')
	plt.title('Experimental vs. All\nPairwise Distances (no cutoff)')
	plt.text(0.9, 0.88,'Mean Distance:\n'+str(round(np.mean(cumDels),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))

	#plot Cth Percentile cutoff
	ax1 = plt.subplot(222)
	ax1.bar(range(seqLen+1), dels, alpha =0.2, color = 'b')
	plt.xlabel('Sequence Deletion Lengths')
	plt.ylabel('Total Deletion Frequency')
	ax2 = ax1.twinx()
	if(args.weighted):
		plt.ylim(ymin,ymax)
	ax2.bar(range(seqLen+1),delsC, edgecolor = 'r', fill = False, linewidth=0.5)
	plt.ylabel(str(C)+'th Percentile Deletion Frequency')
	plt.title('All Experimental vs. '+str(C)+'th Percentile Deletions')
	plt.text(0.9, 0.88,'Mean Distance:\n'+str(round(np.mean(delLensC),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))


	#plot Bth Percentile cutoff
	ax1 = plt.subplot(223)
	ax1.bar(range(seqLen+1), dels, alpha =0.2, color = 'b')
	plt.xlabel('Sequence Deletion Lengths')
	plt.ylabel('Total Deletion Frequency')
	ax2 = ax1.twinx()
	if(args.weighted):
		plt.ylim(ymin,ymax)
	ax2.bar(range(seqLen+1),delsB, edgecolor = 'r', fill = False, linewidth=0.5)
	plt.ylabel(str(B)+'th Percentile Deletion Frequency')
	plt.title('All Experimental vs. '+str(B)+'th Percentile Deletions')
	plt.text(0.9, 0.88,'Mean Distance:\n'+str(round(np.mean(delLensB),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))


	#plot B.5th Percentile cutoff
	ax1 = plt.subplot(224)
	ax1.bar(range(seqLen+1), dels, alpha =0.2, color = 'b')
	plt.xlabel('Sequence Deletion Lengths')
	plt.ylabel('Total Deletion Frequency')
	ax2 = ax1.twinx()
	if(args.weighted):
		plt.ylim(ymin,ymax)
	ax2.bar(range(seqLen+1),delsA, edgecolor = 'r', fill = False, linewidth=0.5)
	plt.ylabel(str(A)+'th Percentile Deletion Frequency')
	plt.title('All Experimental vs. '+str(A)+' Percentile Deletions')
	plt.text(0.9, 0.88,'Mean Distance:\n'+str(round(np.mean(delLensA),2)), horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, bbox=dict(boxstyle="square",facecolor='white'))

	
fig.suptitle(title)
plt.tight_layout()

#append file name, mean and standard deviation to output file
#outF = open('output.txt','a')
#this will later be put into and excel spreadsheet
#outF.write(title+"\t"+str(round(np.mean(cumDels),4))+"\t"+str(round(np.std(cumDels),4))+"\n")
#outF.close()
#plt.show()
plt.savefig(title+'_DelLengthDistribution.pdf')

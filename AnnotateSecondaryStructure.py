##########################################
#	 Varna SVG color script
#
#	   Author: Gregg Rice
#			   gmr@unc.edu
#				Steve modified to plot T1 cleavage data
#				Steve modified to plot additional custom nuc colors and lines between nucs
#				Tom modified to plot lines between nucs colored by their 3D distances
#
# Affiliation: Weeks Lab, UNC Chemistry
#
#		 Date: Oct 6, 2014
#	  Version: 0.93beta2
#
# released under GPL 2.0
##########################################

import math,sys
import numpy as np
from operator import itemgetter
from Bio.PDB import *
import re

class RNA:
	def __init__(self,x):
		self.name = x
		self.num,self.seq,self.ct = self.readCT(x)
		#print len(self.num), len(self.seq), len(self.ct),len(self.dif),len(self.shape)
	
	def __str__(self):
		a = '{ Name= %s, CT= %s }' % (self.name, str(len(self.ct)))
		return a
	
	def readCT(self,z):
		num,seq,bp = [],[],[]
		for i in open(z).readlines()[1:]:
			a = i.rstrip().split()
			num.append(int(a[0])),seq.append(str(a[1])),bp.append(int(a[4]))
		return num,seq,bp
	def readSHAPE(self,z):
		self.shape = []
		for i in open(z).readlines():
			a = i.rstrip().split()
			self.shape.append(float(a[1]))

class structureCoords:
	def __init__(self,x,xrna=False):
		self.name = x
		if xrna:
			self.num, self.pairMap, self.baseMap,self.shape,self.ct,self.scale, self.period = self.parseXRNA(x)
		else:
			print len(self.readVarna(x))
			self.num, self.pairMap, self.baseMap,self.center, self.shape,self.ct,self.scale, self.period = self.readVarna(x)
			
		
		# find the max X and Y coords
		self.maxX = 0.0
		self.minX = 0.0
		self.maxY = 0.0
		self.minY = 0.0
		for i,j in self.pairMap:
			if i < self.minX:
				self.minX = i
			if i > self.maxX:
				self.maxX = i
			if j < self.minY:
				self.minY = j
			if j > self.maxY:
				self.maxY = j
		
		#print self.center
		# initialize color id array
		self.colorIDs = [0]*len(self.num)

		
		self.numMap = self.baseMap
		#self.center = self.pairMap #bad idea here, center stores the location of the center of the helix
		self.diff = False
		self.dms = False
		self.customColors = False
		self.ntFont = 16
		
		#object variables to keep track of which lines are being drawn and 
		#the cutoffs to use
		self.extraLines = False
		self.extraLinesPerc = False
		self.extraLinesDist = False
		self.extraLineData = []
		self.threshList = []
		self.threshA = 0.001
		self.threshB = 0.03
		self.xrna = xrna
		self.pdb = ""
	
	def readSHAPE(self,z):
		self.shape = []
		for i in open(z).readlines():
			a = i.rstrip().split()
			self.shape.append(float(a[1]))
			
	def readColors(self,z):
		self.customColors = True
		#empty color IDs to be filled with custom colors provided by the file
		self.colorIDs = []
		lines = open(z).read().replace("\r","\n").replace("\n\n","\n").split("\n")
		#print len(lines)
		for line in lines:
			ID = line.strip().split()
			if ID:
				ID = ID[0]
				self.colorIDs.append(int(ID))

	def readExtraLines(self,z):
		self.extraLines = True
		splitZ = z.split()
		filename=splitZ[0]
		self.threshA = float(splitZ[1])
		self.threshB = float(splitZ[2])
		self.extraLineData = open(filename).read().replace("\r","\n").replace("\n\n","\n").split("\n")[1:]
		self.extraLineData = [x.split() for x in self.extraLineData]
		#print self.extraLineData
	 
	#used with --percentileLines option.
	#reads in the argument consisting of a 4 column deletion file and list of percentiles
	#these arguments are stored in a string that is parsed  
	def readExtraPercLines(self,z):
		self.extraLinesPerc = True
		#split the string
		splitZ = z.split()
		#grab the file name holding the deletion data (nts to draw lines b/w)
		filename=splitZ[0]
		#loop through rest of elements and load into threshold list
		for element in splitZ[1:]:
			self.threshList.append(float(element))
		#read deletions in object
		self.extraLineData = open(filename).read().replace("\r","\n").replace("\n\n","\n").split("\n")[1:]
		self.extraLineData = [x.split() for x in self.extraLineData]
		#print self.extraLineData
	
	#used with --categoryLines option.
	#reads in the argument consisting of a 5 column deletion file and a percentile
	#these arguments are stored in a string that is parsed  
	def readExtraCatLines(self,z):
		self.extraLinesPerc = True
		#split the string
		splitZ = z.split()
		#grab the file name holding the deletion data (nts to draw lines b/w)
		filename=splitZ[0]
		#loop through rest of elements and load into threshold list
		for element in splitZ[1:]:
			self.threshList.append(float(element))
		#read deletions in object
		self.extraLineData = open(filename).read().replace("\r","\n").replace("\n\n","\n").split("\n")[1:]
		self.extraLineData = [x.split() for x in self.extraLineData]
		#print self.extraLineData
	
	#used with --rawCountLines option.
	#reads in the argument consisting of a 4 column deletion file and list of counts
	#these arguments are stored in a string that is parsed  
	def readExtraCountLines(self,z):
		self.extraLinesPerc = True
		#split the string
		splitZ = z.split()
		#grab the file name holding the deletion data (nts to draw lines b/w)
		filename=splitZ[0]
		#loop through rest of elements and load into threshold list
		for element in splitZ[1:]:
			self.threshList.append(int(element))
		#read deletions in object
		self.extraLineData = open(filename).read().replace("\r","\n").replace("\n\n","\n").split("\n")[1:]
		self.extraLineData = [x.split() for x in self.extraLineData]
		#print self.extraLineData
	
	#used with --distanceLines option.
	#reads in the argument consisting of a 4 column deletion file, a pdb file and list of percentiles
	#these arguments are stored in a string that is parsed   
	def readExtraDistLines(self,z):
		self.extraLinesDist = True
		splitZ = z.split()
		filename=splitZ[0]
		pdbFile=splitZ[1]
		#grab percentile cutoff
		self.threshA = float(splitZ[2])
		self.threshList.append(float(splitZ[2]))
		self.extraLineData = open(filename).read().replace("\r","\n").replace("\n\n","\n").split("\n")[1:]
		self.extraLineData = [x.split() for x in self.extraLineData]
		#store pdb
		self.pdb = loadPDB(pdbFile)
	
	def parseXRNA(self,x):
		import xml.etree.ElementTree as ET
		tree = ET.parse(x)
		root = tree.getroot()

		nucList = root.findall('./Complex/RNAMolecule/')
		nucLists = []
		for i in nucList:
			#print i.tag
			if i.tag == 'NucListData':nucLists.append(i)
			
		#print nucLists
		startNT = int(nucLists[0].get('StartNucID'))
		#print startNT
		num = []
		bases,x_pos,y_pos,x_cen,y_cen = [],[],[],[],[]
		for nt in nucLists[0].text.split('\n'):
			if nt == '':continue
			line = nt.split()
			num.append(startNT)
			startNT+=1
			bases.append(line[0]),x_pos.append(float(line[1])),y_pos.append(float(line[2]))
		
		#print x_pos
		#determine offset
		x_pos,y_pos = np.array(x_pos),-1*np.array(y_pos)
		
		vec = ((x_pos[0]-x_pos[1])**2+(y_pos[0]-y_pos[1])**2)**0.5
		offset = vec/1.5

		#transpose coords to positive space
		xoff = abs(min(x_pos))+offset
		yoff = abs(min(y_pos))+offset
		x_pos = x_pos + xoff
		y_pos = y_pos + yoff
		y_pos=2.3*y_pos
		x_pos=2.3*x_pos
		
		x_cen += xoff
		y_cen += yoff
		center = zip(x_cen,y_cen)
		
		#make expected arrays
		coord = zip(x_pos,y_pos)
		basemap = zip(bases,zip(x_pos,y_pos+offset))
		#print basemap
		
		shape = np.zeros(len(num))
		clist = {'bbbbbb':-999,'999999':-999,'ff0000':1.0, '0':0.0, '000000':0.0, '1c75bc':-0.45,'00ff00':0.45, 'ff9900':0.45, 'f57e1f':0.45}
		for shapeLine in root.findall('./Complex/RNAMolecule/Nuc'):
			nucRange = shapeLine.get('RefIDs')
			preColor = shapeLine.get('Color')
			if not preColor:continue
			try:nucColor = clist[shapeLine.get('Color')]
			except:nucColor = 0.0
			if not nucRange:continue
			for i in nucRange.split(','):
				if len(i.split('-'))==1:
					try:shape[int(i)-1]=nucColor
					except:pass
				else:
					line = i.split('-')
					line = map(int,line)
					for j in range(line[0],line[1]+1):
						try:shape[j-1]=nucColor
						except:pass
		
		shape = map(float,shape)
		period = 20
		#get pairing informationo
		ct = np.zeros(len(num))
		for pair in root.findall('./Complex/RNAMolecule/BasePairs'):
			pairStart,pairLength,pairEnd = pair.get('nucID'),pair.get('length'),pair.get('bpNucID')
			for nt in range(int(pairLength)):
				p5 = int(pairStart)+nt
				p3 = int(pairEnd)-nt
				try:
					ct[p5-1] = p3
					ct[p3-1] = p5
				except:pass
		ct = map(int,ct)
		#print num
		return num,coord,basemap,shape,ct,offset,period
	
	def readVarna(self,x):
		import xml.etree.ElementTree as ET
		tree = ET.parse(x)
		root = tree.getroot()
		
		#initialize some arrays
		offset = 15/2.
		num,bases,x_pos,y_pos,shape,x_cen,y_cen = [],[],[],[],[],[],[]
		
		#read the varna xml file, nucleotides are under bases
		for nt in root.findall('./RNA/bases/nt'):
			num.append(int(nt.get('num')))
			shape.append(float(nt.get('val')))
			base = nt.find('base').text
			for i in nt:
				if i.get('r')=='pos':
					x_pos.append(float(i.get('x')))
					y_pos.append(float(i.get('y')))
					bases.append(base)
				if i.get('r')=='center':
					x_cen.append(float(i.get('x')))
					y_cen.append(float(i.get('y')))
		#determine offset
		x_pos,y_pos = np.array(x_pos),np.array(y_pos)
		
		vec = ((x_pos[0]-x_pos[1])**2+(y_pos[0]-y_pos[1])**2)**0.5
		offset = vec/4
		
		#print x_pos, x_cen

		#transpose coords to positive space
		xoff = abs(min(x_pos))+offset
		yoff = abs(min(y_pos))+offset
		x_pos = x_pos + xoff
		y_pos = y_pos + yoff
		
		x_cen += xoff
		y_cen += yoff
		center = zip(x_cen,y_cen)
		
		#make expected arrays
		coord = zip(x_pos,y_pos)
		basemap = zip(bases,zip(x_pos,y_pos+offset))
		
		#read varna xml file for pairing information
		ct = np.zeros(len(num))
		for pair in root.findall('./RNA/BPs/bp'):
			p5,p3 = int(pair.get('part5')),int(pair.get('part3'))
			#print p5+1,p3+1
			ct[p5] = p3+1
			ct[p3] = p5+1
		ct = map(int,ct)
		
		#get the number period
		period = int(root.find('config').get('numperiod'))
		#print num
		return num,coord,basemap,center,shape,ct,offset,period

class Enzyme:
	
	def __init__(self, filePath):
		self.arrowSizes = []
		self.nums = []
		if filePath != "":
			self.nums, self.arrowSizes = self.parseEnzymeFile(filePath)
	
	def parseEnzymeFile(self,filePath):
		fileIn = open(filePath, "rU")
		fileIn.readline() # skip header line
		nums = []
		arrowSizes = [] 
		for line in fileIn:
			splitLine = line.strip().split()
			nums.append(int(splitLine[0]))
			arrowSizes.append(int(splitLine[5]))
		fileIn.close()
		return nums, arrowSizes
			
	def rotatePoint(self, coords, theta):
		x = coords[0]
		y = coords[1]
		xRot = x*math.cos(-theta)+ y*math.sin(-theta)
		yRot = -x*math.sin(-theta) + y*math.cos(-theta)		   
		#print "x,y= %f, %f\ttheta= %f\t xrot,yrot= %f, %f"%(x,y,theta,xRot,yRot) 
		return [xRot,yRot]
	
	def translatePoint(self, coords, offset):
		return [coords[i] + offset[i] for i in range(len(coords))]
	
	def midpoint(self, c1, c2):
		return [c1[0]+(c2[0]-c1[0])/2, c1[1]+(c2[1]-c1[1])/2]
	
	
	def calcArrowCoords(self, origin, h, w, theta):
		# define initial wedge V
		left = [-w/2,h]
		right = [w/2,h]
		# rotate initial wedge about 0,0
		leftRot = self.rotatePoint(left,theta)
		rightRot = self.rotatePoint(right,theta)
		# translate to given origin
		left = self.translatePoint(leftRot, origin)
		right = self.translatePoint(rightRot, origin)
		# return three coords
		return origin, left, right
		
	def arrowToString(self, pt1, pt2, pt3):
		svgString = '<polygon fill="rgb(80,100,255)" stroke="green" stroke-width="0" points="%f,%f	%f,%f %f,%f" />'%(pt1[0],pt1[1],pt2[0],pt2[1],pt3[0],pt3[1])  
		return svgString  
		
	def drawArrows(self, varna):
		heights = [0,20,40,80]
		width = 10
		arrowString = ""
		for i in range(len(self.nums)):
			num = self.nums[i]
			arrowSize = self.arrowSizes[i]
			loc1 = [varna.baseMap[num-1][1][0],varna.baseMap[num-1][1][1]]
			loc2 = [varna.baseMap[num][1][0],varna.baseMap[num][1][1]]
			loc = self.midpoint(loc1, loc2) # put arrow 3prime of cleaved nuc
			#print loc
			# 0=no arrow, 1=lo, 2=mid, 3=high
			if arrowSize != 0:
				height = heights[arrowSize]
				
				#pt = findFarPoint(num-1,varna)
				#xDiff = pt[0] - loc[0]
				#yDiff = pt[1] - loc[1]
				#theta = math.atan2(yDiff, xDiff)
				
				# assuming clockwise nuc order, find normal angle
				diff = [loc2[n]-loc1[n] for n in [0,1]]
				theta = math.atan2(diff[1], diff[0]) + math.pi
				
				coords = self.calcArrowCoords(loc,height,width,theta)
				arrowString += self.arrowToString(coords[0],coords[1],coords[2])
		
		return arrowString
		
def evalStructures(rna1,rna2):
	# Returns shared basepairs, those only in rna1, those only in rna2
	# n number of accepted pairs x2, p number of predicted pairs x2
	n,p = 0,0
	shared,acceptedOnly,predictedOnly = [],[],[]
	for i in range(len(rna1.ct)):
		#clean out duplicates and nonbasepairs
		if rna1.ct[i] == 0 and rna2.ct[i] ==0:continue
		#count for Sens,PPV
		if rna1.ct[i] != 0: n+=1
		if rna2.ct[i] != 0: p+=1
		#shared bps
		if rna1.ct[i] == rna2.ct[i]:
			if rna1.ct[i] < rna1.num[i]:continue
			if rna2.ct[i] < rna2.num[i]:continue
			shared.append((rna1.num[i],rna1.ct[i]))
			continue
		if rna1.ct[i] != 0 and rna1.ct[i] < rna1.num[i]:
			acceptedOnly.append((rna1.num[i],rna1.ct[i]))
		if rna2.ct[i] != 0 and rna2.ct[i] < rna2.num[i]:
			predictedOnly.append((rna2.num[i],rna2.ct[i]))
	return acceptedOnly,predictedOnly,shared,len(shared)/(float(n)/2),len(shared)/(float(p)/2)	  

def offPoint(p,q,r):
	p,q = np.array(p),np.array(q)
	v_u = (q-p)/(np.sum((q-p)**2)**0.5)
	return v_u*r+p

def newLines(pointSet,locMap,r,struct=False):
	a = []
	#print "locMap"
	#print locMap
	for i,j in pointSet:
		p1 = offPoint(locMap[i-1],locMap[j-1],r)
		p2 = offPoint(locMap[j-1],locMap[i-1],r)
		#check noncanonical
		#if struct:
			#canonical
			#print struct.baseMap[i-1][0],struct.baseMap[j-1][0]
		a.append((p1,p2))
	return a

def getSHAPEcolor(x):
	#if x < -4: return '160,160,160'
	if x < -4: return '180,180,180'
	if x > 0.85: return '255,0,0'
	#if 0.85 >= x >0.4:return '255,210,0'
	if 0.85 >= x >0.4:return '255,164,26'
	if 0.4 >= x > -4: return '1,0,0'

def getDMScolor(x):
	if x < -4: return '180,180,180'
	if x > 0.4: return '255,0,0'
	#if 0.85 >= x >0.4:return '255,210,0'
	if 0.4>= x >0.2:return '255,164,26'
	if 0.2 >= x > -4: return '1,0,0'

def getDiffcolor(x):
	if x < -500:return '160,160,160'
	elif x <= -0.3:return '41,171,226'
	elif x >= 0.3:return '0,210,0'
	else: return '0,0,0'
	
def drawBases(varna):
	bases = varna.baseMap
	shape = varna.shape
	line = ''
	if varna.diff:
		for i in range(len(bases)):
			if abs(shape[i])>0.3 and shape[i] >-500:
				#line += '<text x="%s" y="%s" text-anchor="middle" font-family="Sans-Serif" font-weight="bold" font-size="19" stroke="rgb(%s)" fill="rgb(%s)" >%s</text>' % (bases[i][1][0],bases[i][1][1],getDiffcolor(shape[i]),getDiffcolor(shape[i]),bases[i][0])
				line += '<text x="%s" y="%s" text-anchor="middle" font-weight="bold" font-size="%s" stroke="rgb(%s)" fill="rgb(%s)" >%s</text>' % (bases[i][1][0],bases[i][1][1],varna.ntFont,getDiffcolor(shape[i]),getDiffcolor(shape[i]),bases[i][0])
			else:
				#line += '<text x="%s" y="%s" text-anchor="middle" font-family="Sans-Serif" font-weight="bold" font-size="18" fill="rgb(%s)" >%s</text>' % (bases[i][1][0],bases[i][1][1],getDiffcolor(shape[i]),bases[i][0])
				line += '<text x="%s" y="%s" text-anchor="middle" font-weight="bold" font-size="%s" fill="rgb(%s)" >%s</text>' % (bases[i][1][0],bases[i][1][1],varna.ntFont,getDiffcolor(shape[i]),bases[i][0])
	elif varna.dms:
		for i in range(len(bases)):
			line += '<text x="%s" y="%s" text-anchor="middle" font-weight="bold" font-size="%s" fill="rgb(%s)" >%s</text>' % (bases[i][1][0],bases[i][1][1],varna.ntFont,getDMScolor(shape[i]),bases[i][0])		
	else:
		for i in range(len(bases)):
			#line += '<text x="%s" y="%s" text-anchor="middle" font-family="Sans-Serif" font-size="18" fill="rgb(%s)" >%s</text>\n' % (bases[i][1][0],bases[i][1][1],getSHAPEcolor(shape[i]),bases[i][0])
			if varna.colorIDs[i] == 6 or varna.colorIDs[i] == 3:
				line += '<text x="%s" y="%s" text-anchor="middle" font-family="Arial" font-size="%s" fill="rgb(%s)" >%s</text>\n' % (bases[i][1][0],bases[i][1][1],varna.ntFont,'256,256,256',bases[i][0])
			elif varna.colorIDs[i] == -1:
				#color the nucleotide letter Grey
				line += '<text x="%s" y="%s" text-anchor="middle" font-family="Arial" font-size="%s" fill="rgb(%s)" >%s</text>\n' % (bases[i][1][0],bases[i][1][1],varna.ntFont,'128,128,128',bases[i][0])				
			elif varna.colorIDs[i] == -2:
				#color the nucleotide letter Brown
				line += '<text x="%s" y="%s" text-anchor="middle" font-family="Arial" font-size="%s" fill="rgb(%s)" >%s</text>\n' % (bases[i][1][0],bases[i][1][1],varna.ntFont,'153,76,0',bases[i][0])							
			else:
				#line += '<text x="%s" y="%s" text-anchor="middle" font-family="Arial" font-size="18" fill="rgb(%s)" >%s</text>\n' % (bases[i][1][0],bases[i][1][1],'0,0,0',bases[i][0])
				line += '<text x="%s" y="%s" text-anchor="middle" font-family="Arial" font-size="%s" fill="rgb(%s)" >%s</text>\n' % (bases[i][1][0],bases[i][1][1],varna.ntFont,getSHAPEcolor(shape[i]),bases[i][0])
	return line
			

def findFarPoint(pt,varna):
	#finds a good direction to draw the number label
	# goes through the distance pairs, finds all nts within ~ 63pts and finds the center of mass
	x,y = [],[]
	for i,j in varna.pairMap:
		x.append(i),y.append(j)
	x,y = np.array(x),np.array(y)
	point = np.array((x[pt],y[pt]))
	#print point
	dist = np.sum((point-np.transpose((x,y)))**2,axis=1)
	#print len(dist[dist<5000]),'#'
	#print (x,y)
	cutoff=4000
	length = len(x[dist<cutoff])
	centerMass = np.sum(x[dist<cutoff])/length, np.sum(y[dist<cutoff])/length
	#print str(np.array(centerMass))
	return np.array(centerMass)

def findCWNormalPoint(pt,varna):
	x,y = [],[]
	for i,j in varna.pairMap:
		x.append(i),y.append(j)
	x,y = np.array(x),np.array(y)
	point = np.array((x[pt],y[pt]))
	#print "point: "+str(point)
	try:
		pointAfter = np.array((x[pt+1],y[pt+1]))
	except IndexError:
		pointAfter = np.array((x[pt]+1,y[pt]))
	# assuming clockwise nuc order, find normal angle
	diff = [pointAfter[n]-point[n] for n in [0,1]]
	#print "diff: "+str(diff)
	theta = math.atan2(diff[1], diff[0]) + math.pi/2
	#print "theta: "+str(theta)
	distance = 20
	newX = point[0]+math.cos(theta)*distance
	newY = point[1]+math.sin(theta)*distance
	newPoint = np.array((newX,newY))
	#print "newPoint: "+str(newPoint)+"\n"
	return newPoint 

def drawNums(varna,offset,startNt):
	period = varna.period
	#draw over i
	nums = []
	lines = []
	for i in [1]+range(0,len(varna.num),period)[1:]+[len(varna.num)]:
		key = i-1
		if varna.xrna:
			center = findFarPoint(key,varna)
			#center = findCWNormalPoint(key,varna)
		if not varna.xrna:
			center = np.array(varna.center[key])
		a = np.array(varna.pairMap[key])
		base = np.array(varna.pairMap[key])
		#print "base: "+str(base)
		#print "center: "+str(center)
		norm = np.sum((base-center)**2)**0.5
		#print base, center, norm, varna.scale
		u_vec = (base-center)/norm*varna.scale*7 + base
		nums.append((str(i),map(float,u_vec)))
		
		p1 = offPoint(map(float,u_vec),map(float,base),varna.scale*2)
		p2 = offPoint(map(float,base),map(float,u_vec),varna.scale*2)
		lines.append((p1,p2))
	#add lines connecting base and letters
	line = processConnect(lines,(3,0,0),lineMap=True)
	#add numbering
	for i in range(len(nums)):
		#line += '<text x="%s" y="%s" text-anchor="middle" font-family="Sans-Serif" font-weight="bold" font-size="16" fill="rgb(0,1,0)" >%s</text>' % (nums[i][1][0],nums[i][1][1]+varna.scale,str(int(nums[i][0])+offset))
		#only apply the offset to numbers greater than the starting nt. Otherwise, don't apply offset
		if int(nums[i][0]) >= startNt:
			line += '<text x="%s" y="%s" text-anchor="middle" font-weight="bold" font-size="24" fill="rgb(0,1,0)" >%s</text>\n' % (nums[i][1][0],nums[i][1][1]+varna.scale,str(int(nums[i][0])+offset))
		else:
			line += '<text x="%s" y="%s" text-anchor="middle" font-weight="bold" font-size="24" fill="rgb(0,1,0)" >%s</text>\n' % (nums[i][1][0],nums[i][1][1]+varna.scale,str(int(nums[i][0])))			
	return line

def drawOutline(varna):
	outlineString = '<polyline points= "'
	for nuc in varna.baseMap:
		pt = nuc[1]
		outlineString += '{0},{1} \n'.format(pt[0],pt[1]-4)
	outlineString += '" stroke="rgb(200,200,200)" stroke-width="2" fill="none"/>\n'
	return outlineString

#draws a colored circle behind the nucleotide letters.
#color is white by default or picked by 0,1,2 code from color file --colors
def drawCircles(varna):
	def pct(x):
		return str(x)+'%'
	#0=white,1=red,2=grey,3=blue,4=orange,5=purple,6=black
	circleColors=["white","red","grey","blue","orange","purple","black"]
	outlineString = ''
	#print "customColors = "+str(varna.customColors)
	#print "colorIDs = "+str(varna.colorIDs)
	#print range(len(varna.baseMap))
	#loop through every nucleotide
	for i in range(len(varna.baseMap)):
		nuc = varna.baseMap[i]
		pt = nuc[1]
		#these color selections aren't really used in the code anymore
		col = circleColors[0]
		if varna.customColors==True and i < len(varna.colorIDs):
			col = circleColors[varna.colorIDs[i]]
			#draw an white circle (ie nothing)
			#if varna.colorIDs[i] == 0: 
			#	outlineString += '<circle cx="%f" cy="%f" r="10" stroke="white" stroke-width="2" fill="%s"/>\n'%(pt[0],pt[1]-4,"white")
			#draw a red circle
			if varna.colorIDs[i] == 1: 
				#outlineString += '<circle cx="%f" cy="%f" r="10" stroke="{3}" stroke-width="2" fill="%s"/>\n'%(pt[0],pt[1]-4,col)
				outlineString += '<circle cx="%f" cy="%f" r="11" stroke="white" stroke-width="2" fill="rgb(%s)"/>\n'%(pt[0],pt[1]-6,', '.join(map(pct,(100,0,0))))
			#draw a grey circle
			elif varna.colorIDs[i] == 2:
				#outlineString += '<circle cx="%f" cy="%f" r="10" stroke="white" stroke-width="2" fill="%s"/>\n'%(pt[0],pt[1]-4,"white")
				outlineString += '<circle cx="%f" cy="%f" r="11" fill="rgb(%s)"/>\n'%(pt[0],pt[1]-6,', '.join(map(pct,(65,65,65))))
			#draw a blue circle
			elif varna.colorIDs[i] == 3:
				outlineString += '<circle cx="%f" cy="%f" r="11" stroke="white" stroke-width="2" fill="rgb(%s)"/>\n'%(pt[0],pt[1]-6,', '.join(map(pct,(0,0,100))))
			#draw a orange circle
			elif varna.colorIDs[i] == 4:
				outlineString += '<circle cx="%f" cy="%f" r="11" stroke="white" stroke-width="2" fill="rgb(%s)"/>\n'%(pt[0],pt[1]-8,', '.join(map(pct,(100,65,0))))
			#draw a purple circle
			elif varna.colorIDs[i] == 5:
				outlineString += '<circle cx="%f" cy="%f" r="11" stroke="white" stroke-width="2" fill="rgb(%s)"/>\n'%(pt[0],pt[1]-8,', '.join(map(pct,(100,0,100))))
			#draw a black circle
			elif varna.colorIDs[i] == 6:
				outlineString += '<circle cx="%f" cy="%f" r="11" stroke="white" stroke-width="2" fill="rgb(%s)"/>\n'%(pt[0],pt[1]-6,', '.join(map(pct,(0,0,0))))
	return outlineString
	

def processConnect(pairlist,color,dashed=False, extraLine=False, nonCanon=False, lineMap=False,strokewidth=1.7,opacity=0.95):
	def pct(x):
		return str(x)+'%'
	out = ''
	# print "pairlist: "+str(pairlist)
	for i,j in pairlist:
		if nonCanon == True:
			#first determine if the pair is horizontal or vertical by comparing x and y coordintates
			x1 = i[0]
			x2 = j[0]
			y1 = i[1]
			y2 = j[1]
			cx = (x1+x2)/2
			cy = (y1+y2)/2
			#print "i: "+str(i)
			#print "j: "+str(j)
			
			line = '<circle cx="%s" cy="%s" r="%s" />' % (cx,cy,5)
		else:
			line = '<line x1="%s" y1="%s" x2="%s" y2="%s" ' % (i[0],i[1],j[0],j[1])
			#add a rounded butt for extra lines ie 3D distance lines
			if extraLine == True:
				line += 'stroke-linecap="round" '
			if lineMap == True:
				line += 'stroke="rgb(%s)" stroke-width="0.9" opacity="%.2f" />' % (', '.join(map(pct,color)),opacity)
				out+=line
				continue
			if dashed==False:line += 'stroke="rgb(%s)" stroke-width="%.1f" opacity="%.2f" />' % (', '.join(map(pct,color)),strokewidth,opacity)
			if dashed==True:line += 'stroke="rgb(%s)" stroke-width="1.2" opacity="%.2f" />' %( ', '.join(map(pct,color)),opacity)
		out+=line
	return out

#added functions to allow for 3D Distance calculations
def loadPDB(filename):
	parser = PDBParser()
	structure = parser.get_structure(filename, filename)
	return structure

#returns all residue numbers in a pdb file
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

#returns the 3D distance between 2 nucleotides.
#if one of the nucleotides doesn't exist in the pdb, return a distance of -1
def calculate3DDistance(nt1, nt2, structure, chain, central):
	# input must be a structure object generated with PDBparser
	# as in the loadPDB function or on its own.
	
	#get a list of all residues
	residues = getResNums(structure)
	#make sure both nucleotides exist in pdb
	if nt1 in residues and nt2 in residues:
		
		#if the contacts are from UV or AMT data, calculate distance b/w central point of bases
		if central:
			Base1Coords = getCentralCoord(structure, chain, nt1)
			Base2Coords = getCentralCoord(structure, chain, nt2)
			#prevent residues without base coordinates from causing issues
			if("NO" not in Base1Coords and "NO" not in Base2Coords):
				dist = calcDistance(Base1Coords, Base2Coords)
		else:
			#grab their 2'O atoms and calculate a distance
			atom1 = structure[0][chain][nt1]['O2\'']
			atom2 = structure[0][chain][nt2]['O2\'']
			dist = atom1 - atom2
		return dist
	else:
		#return -1 if one or both of the nucleotides are not in the pdb
		return -1

def drawExtraLines(varna):
	#lineColors = ["orange","red"]
	orange,red = (98,68.6,25), (98.5,14,0)
	out = ""
	if varna.extraLines == True:
		print varna.threshA
		print varna.threshB
		#print str(varna.extraLineData)
		for fields in varna.extraLineData:
			if fields == []: continue
			fromNuc = int(fields[0])-1
			toNuc = int(fields[1])-1
			corrCoeff = float(fields[2])
			if corrCoeff > varna.threshA:
				col = orange
				if corrCoeff > varna.threshB:
					col = red
				fromNucCoord = varna.baseMap[fromNuc][1]
				#print "corr: "+str(corrCoeff)+", color: "+str(col)
				try:
					toNucCoord = varna.baseMap[toNuc][1]
				except IndexError:
					print('Error because it cant find '+fields[1])
					print('The from Nuc is '+fields[0])
					print('With a rate of '+str(corrCoeff))
					raise
				out += processConnect([( ([fromNucCoord[0],fromNucCoord[1]]),
										([toNucCoord[0],toNucCoord[1]])	 )],
										 col,strokewidth=5)
	return out

def drawExtraLinesDist(varna, colorGradient=False, centralPoint=False):
	#color names
	lineColors = ["red","orange","green","grey"]
	#rgb values of colors
	lineColors=[(255,0,0),(255,165,0),(0,255,0),(128,128,128)]
	#svg form (as percents)
	lineColors=[(100,0,0),(100,65,0),(0,100,0),(50,50,50)]
	
	#default opacity
	o = 0.25
	#determine the total number of deletions that will be drawn
	lowestPerc = varna.threshList[-1]
	delsToDraw = (100 - lowestPerc)/100 * len(varna.extraLineData)
	#make lines more transparent if we're going to draw more than 100 deletions
	#if delsToDraw > 100:
	#	o = 0.4
		
	out = ""
	#determine chain of molecule
	chain = ''
	for model in varna.pdb:
		for c in model:
			chain = str(c)
	p = re.compile('Chain id=(\w)')
	m = p.search(chain[1:-1])
	if m:
		chain = m.group(1)
	
	if varna.extraLinesDist == True:
		#print str(varna.extraLineData)
		#make a list of all the rates from the given input file
		rates=[]
		for deletion in varna.extraLineData:
			if deletion == []: continue
			rate = float(deletion[3])
			rates.append(rate)
		#determine rate cutoffs from given percentile
		p = varna.threshA
		varna.threshA = np.percentile(rates,p)
		#print varna.threshList
		#plot deletions above cutoff
		for fields in varna.extraLineData:
			if fields == []: continue
			fromNuc = int(fields[1])-1
			toNuc = int(fields[2])-1
			rate = float(fields[3])
			if rate >= varna.threshA:
				#calculate distance
				dist = calculate3DDistance(fromNuc+1,toNuc+1,varna.pdb, chain, centralPoint)
				#select color based on distance:
				#option 1: color by distance category
				if colorGradient== False:
					if dist > 30:
						#green
						col = lineColors[2]
					elif dist > 15:
						#orange
						col = lineColors[1]
					elif dist > 0:
						#red
						col = lineColors[0]
					else:
						#grey for no distance
						col = lineColors[3]
				#option 2: color as a gradient: close = more red, far = more green
				if colorGradient == True:
					#first handle negative distances ie no distance available
					if dist < 0:
						#set color to grey
						col = (50,50,50)
						#make these lines more transparent
						o = 0.4
					else:
						#calculate distance as a percent of 50, starting at 10.
						#So at 10A the line should be 100% Red and at 50 or more it should be 100% Blue
						#print str(fromNuc)+" "+str(toNuc)+" "+str(dist)
						if dist < 10:
							dist = 10
						dist = dist - 10
						
						percBlue = dist/40.0 * 100.0
						if percBlue > 100:
							percBlue = 100
						#percent Red is the inverse of the percent Blue
						percGreen = 100 - percBlue
						#Kevin addition: divide blue in 1/2 so blue will be darker
						percBlue = percBlue/2
						col = (0,percGreen,percBlue)
						
						o = 0.75
						#print "percent Red: "+str(percRed) + " Percent Blue: "+str(percBlue)+" \n"
						#print col
					
				fromNucCoord = varna.baseMap[fromNuc][1]
				#print "corr: "+str(rate)+", color: "+str(col)
				try:
					toNucCoord = varna.baseMap[toNuc][1]
				except IndexError:
					print('Error because it cant find '+fields[1])
					print('The from Nuc is '+fields[0])
					print('With a rate of '+str(corrCoeff))
					raise
				out = processConnect([( ([fromNucCoord[0],fromNucCoord[1]]),
										([toNucCoord[0],toNucCoord[1]])	 )],
										 col,extraLine=True,strokewidth=4,opacity=o) + out
		return out

def drawExtraLinesPerc(varna):
	#color names
	lineColors = ["cyan","red","orange","green","blue"]
	#rgb values of colors
	lineColors=[(0,255,255),(255,0,0),(255,165,0),(0,255,0),(0,0,255)]
	#svg form (as percents)
	lineColors=[(0,100,100),(100,0,0),(100,65,0),(0,100,0),(0,0,100)]

	#default opacity
	o = 0.6
	#determine the total number of deletions that will be drawn
	lowestPerc = varna.threshList[-1]
	delsToDraw = (100 - lowestPerc)/100 * len(varna.extraLineData)
	
	#make lines more transparent if we're going to draw more than 100 deletions
	#if delsToDraw > 100:
	#	o = 0.4
	out = ""
	#print len(varna.baseMap)
	if varna.extraLinesPerc == True:
		#print str(varna.extraLineData)
		#make a list of all the rates from the given input file
		rates=[]
		for deletion in varna.extraLineData:
			if deletion == []: continue
			rate = float(deletion[3])
			rates.append(rate)
		#determine rate cutoffs from given percentiles
		#print varna.threshList
		for i in range(len(varna.threshList)):
			p = varna.threshList[i]
			varna.threshList[i] = np.percentile(rates,p)
		#print varna.threshList
		#plot deletions above top cutoff
		cutoff = varna.threshList[0]	
		for fields in varna.extraLineData:
			if fields == []: continue
			fromNuc = int(fields[1])-1
			toNuc = int(fields[2])-1
			rate = float(fields[3])
			if rate >= cutoff:
				col = lineColors[0]
				#make sure nt's exist in structure, otherwise, skip
				if fromNuc < len(varna.baseMap) and toNuc < len(varna.baseMap):
					fromNucCoord = varna.baseMap[fromNuc][1]
					#print "corr: "+str(rate)+", color: "+str(col)
					try:
						toNucCoord = varna.baseMap[toNuc][1]
					except IndexError:
						print('Error because it cant find '+fields[1])
						print('The from Nuc is '+fields[0])
						raise
					out = processConnect([( ([fromNucCoord[0],fromNucCoord[1]]),
											([toNucCoord[0],toNucCoord[1]])	 )],
											 col,strokewidth=4,opacity=o) + out
				else:
					continue
		#loop through each cutoff and plot line above
		for i in range(1,len(varna.threshList)):
			upperCutoff = varna.threshList[i-1]
			lowerCutoff = varna.threshList[i]
			for fields in varna.extraLineData:
				if fields == []: continue
				fromNuc = int(fields[1])-1
				toNuc = int(fields[2])-1
				rate = float(fields[3])
				if rate >= lowerCutoff and rate < upperCutoff:
					#pick the color matched to the current percentile cutoff
					col = lineColors[i]
					if fromNuc < len(varna.baseMap) and toNuc < len(varna.baseMap):
						fromNucCoord = varna.baseMap[fromNuc][1]
						#print "corr: "+str(rate)+", color: "+str(col)
						try:
							toNucCoord = varna.baseMap[toNuc][1]
						except IndexError:
							print('Error because it cant find '+fields[1])
							print('The from Nuc is '+fields[0])
							print('With a rate of '+str(corrCoeff))
							raise
						out = processConnect([( ([fromNucCoord[0],fromNucCoord[1]]),
												([toNucCoord[0],toNucCoord[1]])	 )],
												 col,strokewidth=5,opacity=o) + out
					else:
						continue
		return out

#draw lines for the categorically colored lines
def drawExtraLinesCat(varna):
	print "Category Breakdown:"
	print "Green = mutual b/w two sets"
	print "violet = just set 1 (DMSO)"
	print "orange = just set 2 (IA)"
	#color names
	#lineColors = ["green","blue","orange"]
	#rgb values of colors
	#lineColors=[(0,255,0),(0,0,255),(255,165,0)]
	#svg form (as percents)
	#lineColors=[(0,100,0),(0,0,100),(100,65,0)]
	
	#purple alt
	lineColors = ["green","violet","orange"]
	#rgb values of colors
	lineColors=[(0,255,0),(237,102,219),(255,165,0)]
	#svg form (as percents)
	lineColors=[(0,100,0),(93,40,86),(100,65,0)]

	#default opacity
	o = 0.95
	#determine the total number of deletions that will be drawn
	lowestPerc = varna.threshList[-1]
	delsToDraw = (100 - lowestPerc)/100 * len(varna.extraLineData)
	
	#make lines more transparent if we're going to draw more than 100 deletions
	#if delsToDraw > 100:
	#	o = 0.4
	out = ""
	#print len(varna.baseMap)
	if varna.extraLinesPerc == True:
		#print str(varna.extraLineData)
		#make a list of all the rates from the given input file
		rates=[]
		for deletion in varna.extraLineData:
			if deletion == []: continue
			rate = float(deletion[3])
			rates.append(rate)
		#determine rate cutoffs from given percentiles
		#print varna.threshList

		p = varna.threshList[0]
		varna.threshList[0] = np.percentile(rates,p)
		#plot deletions above top cutoff
		cutoff = varna.threshList[0]	
		for fields in varna.extraLineData:
			if fields == []: continue
			fromNuc = int(fields[1])-1
			toNuc = int(fields[2])-1
			rate = float(fields[3])
			category = int(fields[4])
			if rate >= cutoff:
				#print fields
				col = lineColors[category]
				#set line color based on the deletions category

				#make sure nt's exist in structure, otherwise, skip
				if fromNuc < len(varna.baseMap) and toNuc < len(varna.baseMap):
					fromNucCoord = varna.baseMap[fromNuc][1]
					#print "corr: "+str(rate)+", color: "+str(col)
					try:
						toNucCoord = varna.baseMap[toNuc][1]
					except IndexError:
						print('Error because it cant find '+fields[1])
						print('The from Nuc is '+fields[0])
						raise
					out = processConnect([( ([fromNucCoord[0],fromNucCoord[1]]),
											([toNucCoord[0],toNucCoord[1]])	 )],
											 col,strokewidth=3,opacity=o) + out
				else:
					continue
	return out


#draws the lines for the raw count option
def drawExtraLinesCount(varna):
	#color names
	lineColors = ["cyan","red","orange","green","blue"]
	#rgb values of colors
	lineColors=[(0,255,255),(255,0,0),(255,165,0),(0,255,0),(0,0,255)]
	#svg form (as percents)
	lineColors=[(0,100,100),(100,0,0),(100,65,0),(0,100,0),(0,0,100)]

	#default opacity
	o = 0.6
	#determine the total number of deletions that will be drawn
	biggestCount = varna.threshList[-1]	
	#make lines more transparent if we're going to draw more than 100 deletions
	if biggestCount > 100:
		o = 0.4
	out = ""
	#if the user is asking for more lines than there are data, yell about it and return nada
	if biggestCount > len(varna.extraLineData):
		sys.exit("YOU ARE TRYING TO PLOT MORE LINES THAN THERE ARE DATA FOR")

	#print varna.baseMap
	if varna.extraLinesPerc == True:
		#print str(varna.extraLineData)
		#remove any empty data from extra lines list
		while varna.extraLineData[-1] == []:
			varna.extraLineData.pop()
		#make sure that the extraLineData is sorted by rate
		#convert all rates to numbers first so sort workd
		for i in range(len(varna.extraLineData)):
			varna.extraLineData[i][3] = float(varna.extraLineData[i][3])
		varna.extraLineData = sorted(varna.extraLineData, key=lambda x: x[3], reverse=True)
		#print str(varna.extraLineData)
		#print varna.threshList
		#plot deletions above top count cutoff
		#basically just loop through top x lines	
		for i in range(varna.threshList[0]):
			fields = varna.extraLineData[i]
			if fields == []: continue
			fromNuc = int(fields[1])-1
			toNuc = int(fields[2])-1
			col = lineColors[0]
			#make sure nt's exist in structure, otherwise, skip
			if fromNuc < len(varna.baseMap) and toNuc < len(varna.baseMap):
				fromNucCoord = varna.baseMap[fromNuc][1]
				#print "corr: "+str(rate)+", color: "+str(col)
				try:
					toNucCoord = varna.baseMap[toNuc][1]
				except IndexError:
					print('Error because it cant find '+fields[1])
					print('The from Nuc is '+fields[0])
					raise
				out = processConnect([( ([fromNucCoord[0],fromNucCoord[1]]),
										([toNucCoord[0],toNucCoord[1]])	 )],
										 col,strokewidth=5,opacity=o) + out
			else:
				continue
		#loop through each all lines between each set of x to y count cutoffs
		for i in range(1,len(varna.threshList)):
			upperCutoff = varna.threshList[i-1]
			lowerCutoff = varna.threshList[i]
			for j in range(upperCutoff,lowerCutoff):
				fields = varna.extraLineData[j]
				if fields == []: continue
				fromNuc = int(fields[1])-1
				toNuc = int(fields[2])-1
				#pick the color matched to the current percentile cutoff
				col = lineColors[i]
				if fromNuc < len(varna.baseMap) and toNuc < len(varna.baseMap):
					fromNucCoord = varna.baseMap[fromNuc][1]
					#print "corr: "+str(rate)+", color: "+str(col)
					try:
						toNucCoord = varna.baseMap[toNuc][1]
					except IndexError:
						print('Error because it cant find '+fields[1])
						print('The from Nuc is '+fields[0])
						print('With a rate of '+str(corrCoeff))
						raise
					out = processConnect([( ([fromNucCoord[0],fromNucCoord[1]]),
											([toNucCoord[0],toNucCoord[1]])	 )],
											 col,strokewidth=5,opacity=o) + out
				else:
					continue
		return out


def parseArgs():
	import argparse
	prs = argparse.ArgumentParser(description='Colors and optionally compares a VARNA or XRNA file with a reference ct and .SHAPE file')
	prs.add_argument('input',action='store',type=str,help='input file, if xrna add -x flag')
	prs.add_argument('output',action='store',type=str,help='output .svg file')
	prs.add_argument('--ct',action='store',type=str,help='compare structure with a reference ct file')
	prs.add_argument('-x','--xrna',action='store_true',default=False,help='changes input file type to XRNA')
	prs.add_argument('-e','--enzyme',action='store',type=str,help='draw enzymatic cleavage data from file')
	prs.add_argument('-d','--diff',action='store_true',default=False,help='changes chemical probing type to differential, coloring cutoffs +=0.3')
	prs.add_argument('-c','--colors',action='store',type=str,help='color behind nucs with custom colors. 0 for white, 1 for red, 2 for grey, 3 for blue, 4 for orange, 5 for purple, 6 for black.')
	prs.add_argument('-l','--lines',action='store',type=str,help='draw additional lines between certain nucs')
	prs.add_argument('-p','--percentileLines',action='store',type=str,help='similar to lines, but only draws lines above the given percentiles. Takes in a 4 column format file.')
	prs.add_argument('-r','--rawCountLines',action='store',type=str,help='similar to lines, but only draws lines for the top x rates. Can take in multiple counts and color them separately. Takes in a 4 column format file.')	
	prs.add_argument('-cl','--categoryLines',action='store',type=str,help='similar to lines, but only draws lines above a given percentile. Also takes in a 5 column file where the 5th column is numbered 0 to 2 and corresponds to how each line shall be colored.')
	prs.add_argument('--distanceLines',action='store',type=str,help="similar to lines but takes in: 4 column deletion file, pdb, and percentile cutoff. Draws lines of all deletions above the percentile and colors them by 3D distance.")
	prs.add_argument("--centralPoint",action='store_true',default=False,help="Rather than calculating distances between hydroxyl groups, calculate between central point of the base. Use in conjunction with distanceLines.")	
	prs.add_argument('-s','--shape',action='store',type=str,help='overide stored chemical probing values from varna')
	prs.add_argument('--dms',action='store',type=str,help='similar to the shape option. But sets color limits to 0.2 and 0.4')
	prs.add_argument('-o','--offset', action='store',type=int,default=0,help='numbering ofset, adds this to the numbering in the file')
	prs.add_argument('--offsetStart',action='store',type=int,default=1,help='Only add the offset to this nucleotide and beyond')
	prs.add_argument('-n','--nonCanonicalPairs',action='store',type=str,help='input a 2 column list of non canonical base pairs. Will plot circles on top of pair lines to denote non canonical pairing.')
	prs.add_argument('--switch',action='store_true',default=False,help='reverse the pairing coloring scheme')
	prs.add_argument('-f','--ntFontSize',action='store',type=int,default=16,help="Set the font size for all nucleotides. Default value is 16.")
	prs.add_argument('--hideNtNums',action='store_true',help="Hide the nucleotide numbering.")
	o=prs.parse_args()
	return o

def drawNonCanonPairs(nonCanonicalPairsFile,varna):
	out = ""
	#read in pairs
	inFnonCanon = open(nonCanonicalPairsFile,'r')
	lines = inFnonCanon.readlines()
	inFnonCanon.close()
	nonCanonPairNts = []
	for line in lines:
		cols = line.strip().split()
		nonCanonPairNts.append([int(cols[0]),int(cols[1])])
	#get coordinates of pairs
	setScale = varna.scale*2
	print nonCanonPairNts
	nonCanon_pts = newLines(nonCanonPairNts,varna.pairMap,setScale)
	#pass coordinates to process connect
	black = (100,100,100)
	out = processConnect(nonCanon_pts,black,nonCanon = True)
	return out

def hasCT(correct,varna,switch=False):
	ac,pred,both,s,p = evalStructures(correct,varna)
	print 'PPV: ', round(p,2),'SENS: ',round(s,2)
	
	setScale = varna.scale*2
	if varna.xrna:
		setScale = varna.scale*2
	#make lines
	both_pt = newLines(both,varna.pairMap,setScale)
	pred_pt = newLines(pred,varna.pairMap,setScale)
	ac_pt = newLines(ac,varna.pairMap,setScale)
	
	#define colors
	green,red,purple = (0,50,0),(100,0,0),(39,17,56)
	if switch:
		# switch the predicted and accepted colors
		red,purple=(39,17,56),(100,0,0)
	
	#draw lines
	drawnLines = processConnect(pred_pt,purple) + processConnect(ac_pt,red,dashed=True) + processConnect(both_pt,green)
	return drawnLines

def noCT(varna):
	ac,pred,both,s,p = evalStructures(varna,varna)
	
	setScale = varna.scale*2
	if varna.xrna:
		setScale = varna.scale*2
	both_pt = newLines(both,varna.pairMap,setScale,struct=varna)
	black = (0,0,0)
	drawnLines = processConnect(both_pt,black)
	return drawnLines

if __name__ == '__main__':

	arg = parseArgs()

	#read input file, correct ct file
	svgStruct = structureCoords(arg.input,xrna=arg.xrna)
	#if arg.xrna:
	#	 svgStruct = xRNA(arg.input)
	#else:
	#	 svgStruct = varnaSVG(arg.input)

	#overwrite stored chemical probing values
	if arg.shape:
		svgStruct.readSHAPE(arg.shape)
	if arg.dms:
		svgStruct.readSHAPE(arg.dms)
		svgStruct.dms = True
	if arg.diff:
		svgStruct.diff = True

	# do custom colors if given
	if arg.colors:
		svgStruct.readColors(arg.colors)
	#set a custom font size for nucleotides
	svgStruct.ntFont = arg.ntFontSize

	#draw non canonical pairs 
	if arg.nonCanonicalPairs:
		nonCanon = drawNonCanonPairs(arg.nonCanonicalPairs,svgStruct)
	else:
		nonCanon = ""


	# do extra lines
	if arg.lines:
		svgStruct.readExtraLines(arg.lines)

	# draw lines from percentiles
	if arg.percentileLines:
		svgStruct.readExtraPercLines(arg.percentileLines)
  
	#handle lines colored by 3D distances
  	if arg.distanceLines:
  		svgStruct.readExtraDistLines(arg.distanceLines)
  	
  	if arg.rawCountLines:
  		svgStruct.readExtraCountLines(arg.rawCountLines)
  		
  	if arg.categoryLines:
  		svgStruct.readExtraCatLines(arg.categoryLines)
  		
	# read enzymatic cleavage data
	if arg.enzyme:
		enzyme = Enzyme(arg.enzyme)

	# draw nucleotide outline
	outline = drawOutline(svgStruct)

	# draw circles behind nucleotides
	#circles = drawCircles(svgStruct)
	circles = drawCircles(svgStruct)

	# extra lines
	if arg.percentileLines:
		extraLines = drawExtraLinesPerc(svgStruct)
	elif arg.distanceLines:
		extraLines = drawExtraLinesDist(svgStruct,colorGradient=True,centralPoint=arg.centralPoint)
		#extraLines = drawExtraLinesDist(svgStruct)
	elif arg.rawCountLines:
		extraLines = drawExtraLinesCount(svgStruct)
	elif arg.categoryLines:
		extraLines = drawExtraLinesCat(svgStruct)
	else:
		extraLines = drawExtraLines(svgStruct)
	#print extraLines
	# draw lines and compare structures if needed
	if arg.ct:
		correct = RNA(arg.ct)
		drawnLines = hasCT(correct,svgStruct,arg.switch)
	else:
		drawnLines = noCT(svgStruct)

	#construct header and letters
	header = '<svg width="{0} px" height="{1} px" version="1.1" xmlns="http://www.w3.org/2000/svg">'.format(svgStruct.maxX+100, svgStruct.maxY+100)
	#background = '<rect x="0" y="0" width="100%" height="100%" style="fill:rgb(255,255,255)"/>'
	#letters = drawBases(svgStruct) + drawNums(svgStruct,arg.offset,arg.offsetStart)
	
	#don't draw nums
	if arg.hideNtNums:
		letters = drawBases(svgStruct)
	else:
		letters = drawBases(svgStruct) + drawNums(svgStruct,arg.offset,arg.offsetStart)
	
	if arg.enzyme:
		arrows = enzyme.drawArrows(svgStruct)
	else:
		arrows = ""

	#write file
	#out = header + background + outline + circles + letters + drawnLines + extraLines + arrows + '</svg>'
	out = header + outline + circles + letters + drawnLines + nonCanon + extraLines + arrows + '</svg>'

	w = open(arg.output,"w")
	w.write(out)
	w.close()

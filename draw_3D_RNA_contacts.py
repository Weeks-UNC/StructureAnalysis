#!/usr/bin/env python
#draw_3D_RNA_contacts.py

# This program will plot distances on a PDB structure from
# a list of pairs read from a file, in the correlation format

#Create by Tom Christy: 9/20/2016

#The script is designed to take in arguments and run via pymol on the command line.
#Here's how to type the command into the terminal:
#/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL -qcr draw_3D_RNA_contacts.py -- contactFile.txt structure.pdb 
#or if you have the alias working
#pymol -qcr draw_3D_RNA_contacts.py -- contactFile.txt structure.pdb 


from pymol import cmd
import argparse

def parseArgs():
	prs = argparse.ArgumentParser()
	prs.add_argument("ContactFile",type=str,help="3 Column Contact File with header")
	prs.add_argument("pdbFile",type=str,help="Input pdb of RNA structure")
	o = prs.parse_args()
	return o

def getPymolResnums():
    from pymol import cmd,stored
    stored.list=[]
    cmd.iterate("(name O2')","stored.list.append(resi)")
    return stored.list

#parse the arguments
args = parseArgs()

cmd.reinitialize()
cmd.load(args.pdbFile)
#cmd.show('cartoon', 'name p')
#cmd.hide('lines')
cmd.color('slate')
cmd.color('tv_red', 'elem O')

pairs = open(args.ContactFile, 'r')
pairlines = pairs.readlines()
pairs.close()
#remove header
pairlines.pop(0)
pairs = [[int(i.split()[0]), int(i.split()[1])] for i in pairlines]

prefix = "DistHelix_"
cmd.select(prefix, "none")
selected = []
residues = [int(i) for i in getPymolResnums()]
for i in pairs:
	j = i[0]
	k = i[1]
	if j in residues and k in residues:
		for at1 in cmd.index("name O2' in resi "+str(j)):
			for at2 in cmd.index("name O2' in resi "+str(k)):
				cmd.distance(prefix+str(j), "%s`%d"%at1, "%s`%d"%at2)
				cmd.color('red', prefix+str(j))
			# create selections containing relevant nts
		if j not in selected:
			selected.append(j)
			cmd.select(prefix, prefix+" + resi "+str(j))
		if k not in selected:
			selected.append(k)
			cmd.select(prefix, prefix+" + resi "+str(k))
            
cmd.show("spheres", "name O2\' in Link_*")
cmd.set('opaque_background', 'off')
cmd.center()
#create a file name with percentiles included
fileName = args.ContactFile[:-4]
fileName += ".pse"
print "Stuff Drawn!"
print fileName
cmd.save(fileName)
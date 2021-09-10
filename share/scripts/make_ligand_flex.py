#!/usr/bin/env python
import sys,os
from argparse import ArgumentParser 

# parse inputs
parser = ArgumentParser(description="Combine protein flexible parts and covligand")
parser.add_argument("-l", action="store",type=str, help="Covligand.pdb (contain single modified CYS)")
parser.add_argument("-r", action="store",type=str, help="Receptor file")
parser.add_argument("-x", action="store",type=str, help="Flexible residues list")
parser.add_argument("-o", action="store",type=str, help="Output file")
parser.add_argument("-fp",action="store",type=str, help="Flexstr file")

args = parser.parse_args()
ligandfile   = args.l
receptorfile = args.r
flexlistfile = args.x
outputfile   = args.o
fp           = args.fp

# load flex index
flexindex = [ x.strip() for x in open(flexlistfile) ] 

# load covligand lines
covligand_lines  = [ x.strip() for x in open(ligandfile) if "HETATM" in x ]
tmp = covligand_lines[0][21:26]
chain = tmp[0]
index = int(tmp[1:])
ligkey = "%s:%d"%(chain,index)
#print(ligkey)

# load all resi lines
allresidue_lines = [ x.strip() for x in open( receptorfile ) ]

# output
fpofp = open(fp,'w')
fps = set()
with open(outputfile,'w') as ofp:
    for key in flexindex:
        if key == ligkey:
            for line in covligand_lines:
                ofp.write(line+"\n")
        else:
            for line in allresidue_lines:
                if len(line)>=6 and line[:6] in ("ATOM  ","HETATM"):
                    tmp = line[21:26]
                    chain = tmp[0]
                    index = int(tmp[1:])
                    name = line[17:20]
                    if "%s:%d"%(chain,index) == key:
                        fps.add( "%s%d"%(name,index) )
                        ofp.write(line+"\n")
if len(fps) == 0 :
    pass
else:
    fpofp.write(","+",".join(fps))
fpofp.close()

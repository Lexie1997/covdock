#!/usr/bin/env python
import sys,os
import argparse

obabelbin = 'obabel'
chimerabin = 'chimera'
cd_base = os.environ["CONDA_PREFIX"]+"/share/covdock"

# parse input 
parser = argparse.ArgumentParser()
parser.add_argument("-i",type = str, help = "Input sdf", required = True)
parser.add_argument("-r",type = str, help = "Receptor file", required = True)
parser.add_argument("-s",type = str, help = "Receptor CYS specifer", required = True)
parser.add_argument("-l",type = str, help = "Flex residues list file" ,required = True)
parser.add_argument("-w",type = str, help = "Working dir", default = "default_covdock_results" )
args = parser.parse_args()
infile = args.i
receptor = args.r
flexlist = args.l
specifer = args.s
wd = args.w

assert infile[-4:] == ".sdf"
name = infile[:-4]

# cp files
os.system(f"mkdir {wd}")
os.system(f"cp {receptor} {wd}")
os.system(f"cp {flexlist} {wd}")
os.system(f"cp {infile}   {wd}")
os.chdir(f"{wd}")

if True:
    # make mol2
    os.system( f"{obabelbin} -isdf {name}.sdf -omol2 -O tmp.mol2" )
    os.system(f"cat tmp.mol2 |sed 's/XX/S/' |sed 's/Du/SH/' > {name}.mol2")
    with open("chimera.com",'w') as ofp:
        ofp.write(f"open {name}.mol2 \n")
        ofp.write("sel @/serialNumber=1 \n")
        ofp.write("addh spec sel\n")
        ofp.write("write format mol2 #0 tmp2.mol2\n")
        ofp.write("stop\n")
    os.system(f"{chimerabin} --nogui chimera.com")
    for line in open("tmp2.mol2"):
        if "HS" in line:
            h_index = int( line.split()[0] ) 
            s_index = 1
    # Done.  S serial number is 1. HS atom name is HS 
    os.system(f"mv tmp2.mol2 {name}.mol2")

    # make sh
    with open("run.sh","w") as ofp:
        ofp.write("source ~/.bashrc\n")
        ofp.write("conda activate covdock\n")
        ofp.write( "export RECEPTOR='fixed0'\n")
        ofp.write( f"export CYSSPECIFER='{specifer}'\n")
        ofp.write( "export FLEXINDEXFILE='flex.list'\n")
        ofp.write(f"export LIGAND='{name}'\n")
        ofp.write(f"export LIGANDINDICES='{h_index},{s_index}'\n")
        ofp.write(f"bash {cd_base}/share/scripts/run_covalent_dock.sh \n")
        ofp.write("\n")

    # submit 
    os.system("bash run.sh")

    # clean
    os.system("rm *.map")
    os.system(f"rm chimera.com empty {name}.mol2 fixed0_flex.pdbqt fixed0.glg fixed0.gpf fixed0.pdbqt fixed0_rigid.maps.fld fixed0_rigid.maps.xyz fixed0_rigid.pdbqt flex.fp ligcovalent_fixed0.dlg ligcovalent_fixed0.dpf ligcovalent_flex.pdbqt ligcovalent.pdbqt ligcovalent_rigid.pdbqt ligcovalent_single.pdb tmp.mol2")
    os.system("rm residues.log ")

    # grep log
    os.system('grep -E  "Run =|Estimated Free Energy of Binding" my_docking.pdbqt > score.log')

    print("Done")

os.chdir(f"..")


#!/usr/bin/env python
import sys, os
import argparse
import math

obabelbin = "obabel"
chimerabin = "chimera"
cd_base = os.environ["CONDA_PREFIX"] + "/share/covdock"

predefinedbondlength = {"S": 2.0, "C": 1.77}
SH_LENGTH = 1.34
angle = 108.5 * math.pi / 180.0

# parse input
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="Input sdf", required=True)
parser.add_argument("-r", type=str, help="Receptor file", required=True)
parser.add_argument("-s", type=str, help="Receptor CYS specifer", required=True)
parser.add_argument("-l", type=str, help="Flex residues list file", required=True)
parser.add_argument(
    "-w", type=str, help="Working dir", default="default_covdock_results"
)
args = parser.parse_args()
infile = args.i
receptor = args.r
assert receptor[-4:] == ".pdb"
receptorname = os.path.basename(receptor)[:-4]
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
    os.system(f"{obabelbin} -isdf {name}.sdf -omol2 -O tmp.mol2")

    if True:  # Use geometry methods to add attached hydrogen atoms
        # Load atoms: XX, C-XX
        if True:
            xxatomline = None
            xxindex = None
            cxxatomline = None
            cxxindex = None
            cxxelement = None
            bondlength = None

            # find header lines
            headerlines = list()
            flag = False
            for line in open("tmp.mol2"):
                if "@<TRIPOS>MOLECULE" in line:
                    flag = True
                    continue
                if "@<TRIPOS>ATOM" in line:
                    flag = False
                    continue
                if flag:
                    headerlines.append(line)

            # find atom lines
            atomlines = list()
            flag = False
            for line in open("tmp.mol2"):
                if "@<TRIPOS>ATOM" in line:
                    flag = True
                    continue
                if "@<TRIPOS>BOND" in line:
                    flag = False
                    continue
                if flag:
                    atomlines.append(line)

            # find bond lines
            bondlines = list()
            flag = False
            for line in open("tmp.mol2"):
                if "@<TRIPOS>BOND" in line:
                    flag = True
                    continue
                if flag and "@<TRIPOS>" in line:
                    flag = False
                    continue
                if flag:
                    bondlines.append(line)
            # find XX
            for line in atomlines:
                if "XX" in line and line.split()[1] == "XX":
                    xxatomline = line
                    xxindex = int(line.split()[0])
                    s_index = xxindex
            # find bond-XX atom index
            for line in bondlines:
                bonda1 = int(line.split()[1])
                bonda2 = int(line.split()[2])
                if xxindex == bonda1:
                    cxxindex = bonda2
                    break
                elif xxindex == bonda2:
                    cxxindex = bonda1
                    break
            # find bond-xx atom line
            for line in atomlines:
                index = int(line.split()[0])
                if index == cxxindex:
                    cxxatomline = line
                    break
            # determine cxx atom element ( S or C )
            cxxtype = cxxatomline.split()[5]
            if "S" in cxxtype:
                cxxelement = "S"
            else:
                cxxelement = "C"

            # determine bond length
            bondlength = predefinedbondlength[cxxelement]

            # extract coordinate
            x = float(xxatomline.split()[2])
            y = float(xxatomline.split()[3])
            z = float(xxatomline.split()[4])
            xxcoord = (x, y, z)
            x = float(cxxatomline.split()[2])
            y = float(cxxatomline.split()[3])
            z = float(cxxatomline.split()[4])
            cxxcoord = (x, y, z)

        # Compute new S coordinate
        cs = (
            xxcoord[0] - cxxcoord[0],
            xxcoord[1] - cxxcoord[1],
            xxcoord[2] - cxxcoord[2],
        )
        mod = (cs[0] ** 2 + cs[1] ** 2 + cs[2] ** 2) ** 0.5
        i = (cs[0] / mod, cs[1] / mod, cs[2] / mod)
        di = (i[0] * bondlength, i[1] * bondlength, i[2] * bondlength)
        scoord = (cxxcoord[0] + di[0], cxxcoord[1] + di[1], cxxcoord[2] + di[2])

        # Solve HS coordinate
        l_sc = bondlength
        l_sh = SH_LENGTH

        sc = ( -di[0], -di[1], -di[2])  # sh = (x,y,0) 
        xsc = sc[0]
        ysc = sc[1]
        const = l_sc * l_sh * math.cos(angle)

        a = 1 + (xsc / ysc) ** 2
        b = ( -2 * const * xsc ) / (ysc ** 2)
        c = (const / ysc) ** 2 - l_sh ** 2

        print("const:",const)
        print("l_sc:",l_sc,"l_sh:",l_sh)
        print("xsc:",xsc,"ysc:",ysc)
        print("A:",a,"B:",b,"C:",c)

        x = None
        y = None
        if ysc == 0.0:
            y = 0
            x = bondlength
        else:
            x = ( 0 - b + math.sqrt(b * b - 4 * a * c) ) / (2 * a)
            # y = math.sqrt( l_sh**2 - x**2 )
            y = (const-x*xsc)/ysc

        print("x:",x,"y:",y)
        print("cs:",cs)

        hcoord = (scoord[0] + x, scoord[1] + y, scoord[2] + 0)

        # save output mol2 file
        h_index = len(atomlines) + 1
        h_line = f"{h_index} HS       {hcoord[0]}   {hcoord[1]}  {hcoord[2]}  H      1  UNL1        0.0000\n"
        items = xxatomline.split()
        items[1] = "SH"
        items[2] = str(scoord[0])
        items[3] = str(scoord[1])
        items[4] = str(scoord[2])
        items[5] = "S.2"
        s_line = " ".join(items) + "\n"
        atomlines[xxindex - 1] = s_line

        bond_index = len(bondlines) + 1
        h_bond_line = f"{bond_index}   {h_index}   {xxindex}   1"
        count_line = headerlines[1]
        atom_count = int(count_line.split()[0])
        bond_count = int(count_line.split()[1])
        assert atom_count == len(atomlines)
        assert bond_count == len(bondlines)
        new_count_line = f"{atom_count+1} {bond_count+1} 0 0 0 \n"
        headerlines[1] = new_count_line
        with open(f"{name}.mol2", "w") as ofp:
            ofp.write("@<TRIPOS>MOLECULE\n")
            for line in headerlines:
                ofp.write(line)
            ofp.write("@<TRIPOS>ATOM\n")
            for line in atomlines:
                ofp.write(line)
            ofp.write(h_line)
            ofp.write("@<TRIPOS>BOND\n")
            for line in bondlines:
                ofp.write(line)
            ofp.write(h_bond_line)

    # make sh
    with open("run.sh", "w") as ofp:
        ofp.write("source ~/.bashrc\n")
        ofp.write("conda activate covdock\n")
        ofp.write(f"export RECEPTOR='{receptorname}'\n")
        ofp.write(f"export CYSSPECIFER='{specifer}'\n")
        ofp.write("export FLEXINDEXFILE='flex.list'\n")
        ofp.write(f"export LIGAND='{name}'\n")
        ofp.write(f"export LIGANDINDICES='{h_index},{s_index}'\n")
        ofp.write(f"bash {cd_base}/share/scripts/run_covalent_dock.sh \n")
        ofp.write("\n")

    # submit
    os.system("bash run.sh")

    # clean
    os.system("rm *.map")
    os.system(
        f"rm chimera.com empty {name}.mol2 {receptorname}_flex.pdbqt {receptorname}.glg {receptorname}.gpf {receptorname}.pdbqt {receptorname}_rigid.maps.fld {receptorname}_rigid.maps.xyz {receptorname}_rigid.pdbqt flex.fp ligcovalent_{receptorname}.dlg ligcovalent_{receptorname}.dpf ligcovalent_flex.pdbqt ligcovalent.pdbqt ligcovalent_rigid.pdbqt ligcovalent_single.pdb tmp.mol2"
    )
    os.system("rm residues.log ")

    # grep log
    os.system(
        'grep -E  "Run =|Estimated Free Energy of Binding" my_docking.pdbqt > score.log'
    )

    print("Done")

os.chdir(f"..")

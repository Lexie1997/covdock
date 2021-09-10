#!/bin/bash


export CD_BASE="$CONDA_PREFIX/share/covalent_dock"

# RECEPTOR='fixed0'
# CYSSPECIFER='A:CYS150'
# FLEXINDEXFILE='flex.list'
# LIGAND='313'
# LIGANDINDICES='32,1'

#change path to your own
python $CD_BASE/share/adCovalentDockResidue/adcovalent/prepareCovalent.py --ligand $LIGAND.mol2 --ligindices $LIGANDINDICES --receptor $RECEPTOR.pdb \
    --residue $CYSSPECIFER --outputfile ligcovalent_single.pdb
echo '1---ligcovalent.pdb prepared...'

# Combine
python $CD_BASE/share/scripts/make_ligand_flex.py -l ligcovalent_single.pdb -r $RECEPTOR.pdb -x $FLEXINDEXFILE -o ligcovalent.pdb -fp flex.fp
FLEXPART=$(cat flex.fp)
echo '1.5--done'

#protein.pdqt
$CD_BASE/share/mgltools_x86_64Linux2_1.5.6/bin/pythonsh $CD_BASE/share/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
	-r $RECEPTOR.pdb
echo '2---done'

#ligcovalent.pdbqt
$CD_BASE/share/mgltools_x86_64Linux2_1.5.6/bin/pythonsh $CD_BASE/share/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
	-r ligcovalent.pdb
echo '3---done'

#protein_flex.pdbqt protein_rigid.pdbqt	
$CD_BASE/share/mgltools_x86_64Linux2_1.5.6/bin/pythonsh $CD_BASE/share/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_flexreceptor4.py \
	-r $RECEPTOR.pdbqt -s $RECEPTOR:$CYSSPECIFER$FLEXPART

#ligcovalent_rigid.pdbqt ligcovalent_flex.pdbqt
$CD_BASE/share/mgltools_x86_64Linux2_1.5.6/bin/pythonsh $CD_BASE/share/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_flexreceptor4.py \
	-r ligcovalent.pdbqt -s ligcovalent:$CYSSPECIFER$FLEXPART
echo '4---done'

#protein.gpf	
$CD_BASE/share/mgltools_x86_64Linux2_1.5.6/bin/pythonsh $CD_BASE/share/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py \
	-r $RECEPTOR\_rigid.pdbqt -x ligcovalent_flex.pdbqt -l ligcovalent_flex.pdbqt -y -I 20 -o $RECEPTOR.gpf

#lig	
$CD_BASE/share/mgltools_x86_64Linux2_1.5.6/bin/pythonsh $CD_BASE/share/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py \
	-r $RECEPTOR\_rigid.pdbqt -x ligcovalent_flex.pdbqt -l ligcovalent_flex.pdbqt -o ligcovalent_$RECEPTOR.dpf -p move='empty'
echo '5---done'

#touch empty
echo "rmsdatoms all" > empty
sed -i 's/unbound_model extended/unbound_energy 0.0/g' ligcovalent_$RECEPTOR.dpf
#sed -i 's/ga_run 10/ga_run 30/g' ligcovalent_$RECEPTOR.dpf

autogrid4 -p $RECEPTOR.gpf -l $RECEPTOR.glg
sleep 2

autodock4 -p ligcovalent_$RECEPTOR.dpf -l ligcovalent_$RECEPTOR.dlg
sleep 2

grep '^DOCKED' ligcovalent_$RECEPTOR.dlg | cut -c9- > my_docking.pdbqt

echo Job finished, good luck!



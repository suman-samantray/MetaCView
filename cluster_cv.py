#!/opt/local/bin/python

import sys
import numpy, numpy.linalg
import MDAnalysis,MDAnalysis.analysis.hbonds

import CollectiveVariables


## read in cluster pdb file
u=MDAnalysis.Universe(sys.argv[1])

pro_select="protein"
back_select="backbone"

## check for unusual residues
if 'S1P' in u.residues.resnames:
    pro_select="protein or resname S1P L5A DML TPL"
    back_select="backbone or (resname S1P L5A DML TPL and name N CA C O)"


try:
    ff=sys.argv[2]
except:
    ff='charmm'

## make selections
protein=u.select_atoms(pro_select)
backbone=u.select_atoms(back_select)

## find alpha-helix hydrogen bonds
alpha_hbonds=CollectiveVariables.get_alpha_hbonds(protein,ff=ff)

## find phi and psi angles
phi_angles,psi_angles=CollectiveVariables.get_phi_psi(protein)

## get salt bridges
zeta_sel,coo_sel=CollectiveVariables.find_salt_bridges(protein)

## get number of protein and backbone hydrogen bonds (using built into MDAnalysis routines)
h=MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,pro_select,pro_select)
h.run()
nh=h.count_by_time()

hb=MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,back_select,back_select)
hb.run()
nhb=hb.count_by_time()

## loop over cluster structures
for i,ts in enumerate(u.trajectory):

    ## get radius of gyration
    rgyr=protein.radius_of_gyration()

    ## calculate gyration tensor
    gyrT,gyrVal,gyrVec=CollectiveVariables.calc_gyration_tensor(protein)
    val_max=max(gyrVal)
    val_min=min(gyrVal)
    val_mid=sum(gyrVal)-val_max-val_min

    val_max=numpy.sqrt(val_max)
    val_mid=numpy.sqrt(val_mid)
    val_min=numpy.sqrt(val_min)

    ## get number of alpha-helix H-bonds
    nh_alpha=CollectiveVariables.sum_alpha_hbonds(protein,alpha_hbonds)

    ## get dihedral angle correlation function
    dih_offset=CollectiveVariables.beta_dih_offset(phi_angles,psi_angles)

    ## get Cgamma contacts
    cg_contacts=CollectiveVariables.sum_cgamma_contacts(protein)

    sb_contacts=CollectiveVariables.sum_contacts(zeta_sel,coo_sel,4.5)

    i+=1
    print i,rgyr,val_max,val_mid,val_min,nh_alpha,dih_offset,cg_contacts,sb_contacts,nh[i-1][1],nhb[i-1][1]

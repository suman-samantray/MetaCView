#!/opt/local/bin/python
from __future__ import print_function
import sys
import numpy, numpy.linalg
import MDAnalysis,MDAnalysis.analysis.hbonds
sys.path.append('/Users/sumansamantray/Desktop/static-bias')
import CollectiveVariables

## read in simulation trajectory
u=MDAnalysis.Universe(sys.argv[1],sys.argv[2])

dt=float(sys.argv[4])

## check for force field
try:
    ff=sys.argv[5]
except:
    ff='charmm'

## make selections
pro_select="protein"
back_select="backbone"

## check for unusual residues
if 'S1P' in u.residues.resnames:
    pro_select="protein or resname S1P L5A DML TPL"
    back_select="backbone or (resname S1P L5A DML TPL and name N CA C O)"

## make selections
protein=u.select_atoms(pro_select)
backbone=u.select_atoms(back_select)

## find alpha-helix hydrogen bonds
alpha_hbonds=CollectiveVariables.get_alpha_hbonds(protein,ff=ff)

## find three-ten-helix hydrogen bonds
three_ten_hbonds=CollectiveVariables.get_three_ten_hbonds(protein,ff=ff)

## find phi and psi angles
phi_angles,psi_angles=CollectiveVariables.get_phi_psi(protein)

## get cgamma_atoms
cgamma_atoms=CollectiveVariables.get_cgamma_atoms(protein)

## get salt bridges
zeta_sel,coo_sel=CollectiveVariables.find_salt_bridges(protein,termini=True)

## loop over saved timesteps
ouf=open(sys.argv[3],'w')

for i,ts in enumerate(u.trajectory):

    ## get radius of gyration
    rgyr=protein.radius_of_gyration()

    ## calculate gyration tensor
    gyrT,gyrVal,gyrVec=CollectiveVariables.calc_gyration_tensor(protein,pbc=True,boxx=ts.dimensions)
    val_max=max(gyrVal)
    val_min=min(gyrVal)
    val_mid=sum(gyrVal)-val_max-val_min

    val_max=numpy.sqrt(val_max)
    val_mid=numpy.sqrt(val_mid)
    val_min=numpy.sqrt(val_min)

    ## get number of alpha-helix H-bonds
    nh_alpha=CollectiveVariables.sum_alpha_hbonds(protein,alpha_hbonds,pbc=True,boxx=ts.dimensions)

    ## get number of 3/10-helix H-bonds
    nh_three_ten=CollectiveVariables.sum_alpha_hbonds(protein,three_ten_hbonds,pbc=True,boxx=ts.dimensions)

    ## get dihedral angle correlation function
    dih_offset=CollectiveVariables.beta_dih_offset(phi_angles,psi_angles)

    ## get Cgamma contacts
    cg_contacts=CollectiveVariables.sum_cgamma_contacts(protein,cgamma_atoms,pbc=True,boxx=ts.dimensions)

    ## get salt bridges
    sb_contacts=CollectiveVariables.sum_contacts(zeta_sel,coo_sel,4.5,pbc=True,boxx=ts.dimensions)

    print ("%12.6f "*10 % (i*dt,rgyr,val_max,val_mid,val_min,nh_alpha,nh_three_ten,dih_offset,cg_contacts,sb_contacts),end=' ')
    print ("%12.6f "*10 % (i*dt,rgyr,val_max,val_mid,val_min,nh_alpha,nh_three_ten,dih_offset,cg_contacts,sb_contacts),end=' ',file=ouf)
    for phi,psi in zip(phi_angles,psi_angles):
        print ("%12.6f "*2 % (phi.dihedral(),psi.dihedral()),end='')
        print ("%12.6f "*2 % (phi.dihedral(),psi.dihedral()),end='',file=ouf)
    print ()
    print ('',file=ouf)

#!/opt/local/bin/python

import numpy, sys
import MDAnalysis
import CollectiveVariables

## chack command line arguments
if len(sys.argv)<4:
    print 'Usage: generate_plumed_cv.py <pdb/gro file> <out file> <ff> <CV1> <CV2> ...'
    print 'Where'
    print 'ff is one of'
    print 'amber/charmm/gromos'
    print 'CV1/2/.. is one of: '
    print 'alpha_hbonds'
    print 'beta_dih'
    print 'cgamma_contacts'


## open input file and create universe
u=MDAnalysis.Universe(sys.argv[1])

## open output file (plumed input file, missing some of the header and footer information)
ouf=open(sys.argv[2],'w')

ff=sys.argv[3]

## make list of collective variables
CVs=sys.argv[4:]

## make selections
protein=u.select_atoms('protein')

## alpha-helix CV
if 'alpha_hbonds' in CVs:

    ## get alpha-helix hydrogen bonds
    alpha_hbonds=CollectiveVariables.get_alpha_hbonds(protein,ff=ff)

    ## write into file
    print >> ouf, "CONTACTMAP ..."
    for i,hb in enumerate(alpha_hbonds):
        strn="ATOMS%d=%d,%d" % (i+1,hb[0]+1,hb[1]+1)
        print >> ouf, strn
    print >> ouf, "SWITCH={RATIONAL R_0=0.25 NN=8 MM=12}"
    print >> ouf, "SUM"
    print >> ouf, "LABEL=alpha_hbonds"
    print >> ouf, "... CONTACTMAP"
    print >> ouf

## beta-strand CV
if 'beta_dih' in CVs:

    ## get list of phi and psi angles
    phi_angles,psi_angles=CollectiveVariables.get_phi_psi(protein)

    print >>ouf, "ALPHABETA ..."
    i=1
    for phi_a,psi_a in zip(phi_angles,psi_angles):

        i1=phi_a.atoms[0].index+1
        i2=phi_a.atoms[1].index+1
        i3=phi_a.atoms[2].index+1
        i4=phi_a.atoms[3].index+1
        strn="ATOMS%d=%d,%d,%d,%d REFERENCE%d=-2.36" % (i,i1,i2,i3,i4,i)
        print >> ouf, strn
        i1=psi_a.atoms[0].index+1
        i2=psi_a.atoms[1].index+1
        i3=psi_a.atoms[2].index+1
        i4=psi_a.atoms[3].index+1
        strn="ATOMS%d=%d,%d,%d,%d REFERENCE%d=2.36" % (i+1,i1,i2,i3,i4,i+1)
        print >> ouf, strn
        i+=2
    print >> ouf, "LABEL=dih_offset"
    print >> ouf, "... ALPHABETA"
    print >> ouf

## cgamma contacts
if 'cgamma_contacts' in CVs:

    ## get number of residues
    nres=len(protein.residues)

    ## get list of cgamma atoms
    cgamma_list=[]
    for ires in range(nres):
        
        sel = u.select_atoms("resname VAL and name CG1 ").residues
        list(sel) == list(protein.residues)

        if False:
            icg=protein.residues[ires].atoms.CG.index

        else:
            icg=protein.residues[ires].atoms.CG1.index
   
        cgamma_list.append(icg)   

        nam = u.select.atoms("resname VAL and name CG2").residues 
        list(nam) == list(protein.residues)

        if True:
            icg=protein.residues[ires].atoms.CG2.index
        else:
            print "error"
        cgamma_list.append(icg)
        cgamma_list.sort(reverse = False) 
    ncg=len(cgamma_list)
    print >> ouf, "CONTACTMAP ..."
    ii=0

    ## double loop over Cgamma atoms
    for icg in range(ncg-1):
        atom1=cgamma_list[icg]
        for jcg in range(icg+1,ncg):
            atom2=cgamma_list[jcg]
            ii+=1
            strn="ATOMS%d=%d,%d" % (ii,atom1+1,atom2+1)
            print >> ouf, strn

    print >> ouf, "SWITCH={RATIONAL R_0=0.5 NN=8 MM=12}"
    print >> ouf, "SUM"
    print >> ouf, "LABEL=cgamma_contacts"
    print >> ouf, "... CONTACTMAP"
    print >> ouf

#print "METAD ARG=cm,ab SIGMA=0.4,0.1 HEIGHT=4.0 PACE=500 TEMP=300 BIASFACTOR=20 GRID_SPARSE LABEL=metad"
#print "PRINT ARG=cm,ab,metad.bias,uwall.bias STRIDE=100 FILE=COLVAR.dat"

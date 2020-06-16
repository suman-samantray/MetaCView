#!/opt/local/bin/python

"""
File containing functions for setting up and calculating
protein collective variables
"""

import sys
import numpy, numpy.linalg
import MDAnalysis,MDAnalysis.analysis.hbonds

## ----------------------------------------------------------------------------
## some global parameters

## indices in rational switch function
nn=8.0
mm=12.0

## phi and psi angles that define beta strand
beta_phi_ref=-(135.0/180.0)*numpy.pi
beta_psi_ref=(135.0/180.0)*numpy.pi

## cut-offs
rcgamma_cut=5.0

## ----------------------------------------------------------------------------

def anint(x):

    """
    anint

    Function to round to nearest integer

    Input
    - x: real number

    Output
    - anint: nearest integet to x
    """

    x1=numpy.floor(x)
    x2=x1+1
    if (x2-x)<(x-x1):
        return x2
    else:
        return x1

## ----------------------------------------------------------------------------

def switch(r,r0):

    """
    switch

    Rational switch function used for determing contacts
    S(r)=(1-(r/r0)^n)/(1-(r/r0)^m)
    Assuming n<m goes to 1 at small r and 0 at large r
    At r=r0 goes to n/m

    Input
    - r current distance
    - r0 reference distance

    Output
    - switch (value between 0 and 1)

    """

    ## check that r is positive
    if r>0:

        x=r/r0

        ## make sure r/=r0
        if abs(x-1.0)>0.01:
            x1=x**nn
            x2=x**mm
            return (1-x1)/(1-x2)
        else:

            ## r=r0 special case
            return nn/mm
    else:

        ## r<0 special case
        return 1

## ----------------------------------------------------------------------------

def get_alpha_hbonds(protein,ff='charmm'):

    """
    get_alpha_hbonds

    function to find atoms that would form hydrogen
    bonds in alpha-helix
    O(i)->H(i+4)
    This finds and omits proline residues

    Input
    - protein: MDAnalysis atom selection corresponding to protein of interest

    Output
    - alpha_hbonds: list consisting of pairs (O(i),H(i+4))
    """

    ## get number of residues
    nres=len(protein.residues)

    ## loop over residues
    alpha_hbonds=[]
    for ires in range(nres-4):
        iO=protein.residues[ires].atoms.O.index
        if protein.residues[ires+4].resname=='PRO':
            continue
        if ff=='charmm':
            iH=protein.residues[ires+4].atoms.HN.index
        elif ff=='amber':
            iH=protein.residues[ires+4].atoms.H.index
        elif ff=='gromos':
            iH=protein.residues[ires+4].atoms.H.index
        alpha_hbonds.append((iO,iH))
    return alpha_hbonds

## ----------------------------------------------------------------------------

def get_three_ten_hbonds(protein,ff='charmm'):

    """
    get_three_ten_hbonds

    function to find atoms that would form hydrogen
    bonds in 3/10-helix
    O(i)->H(i+3)
    This finds and omits proline residues

    Input
    - protein: MDAnalysis atom selection corresponding to protein of interest

    Output
    - alpha_hbonds: list consisting of pairs (O(i),H(i+4))
    """

    ## get number of residues
    nres=len(protein.residues)

    ## loop over residues
    three_ten_hbonds=[]
    for ires in range(nres-3):
        iO=protein.residues[ires].atoms.O.index
        if protein.residues[ires+3].resname=='PRO':
            continue
        if ff=='charmm':
            iH=protein.residues[ires+3].atoms.HN.index
        elif ff=='amber':
            iH=protein.residues[ires+3].atoms.H.index
        elif ff=='gromos':
            iH=protein.residues[ires+3].atoms.H.index
        three_ten_hbonds.append((iO,iH))
    return three_ten_hbonds

## ----------------------------------------------------------------------------

def get_cgamma_atoms(protein):

    cgamma_atoms=[]
    for residue in protein.residues:
        cgamma_atoms.append([])
        if 'CG' in residue.atoms.names:
            cgamma_atoms[-1].append(residue.atoms.CG.id)
        if 'CG1' in residue.atoms.names:
            cgamma_atoms[-1].append(residue.atoms.CG1.id)
        if 'CG2' in residue.atoms.names:
            cgamma_atoms[-1].append(residue.atoms.CG2.id)

    return cgamma_atoms

## ----------------------------------------------------------------------------

def sum_alpha_hbonds(protein,alpha_hbonds,pbc=False,boxx=None):

    """
    sum_alpha_hbonds

    Function to count number of alpha-helix hydrogen bonds
    are present in a protein structures

    Count H-bonds using rational switch function

    Input
    - protein: MDanalysis atom selection consisting of protein
    - alpha_hbonds: list of (O,H) atoms that form alpha-helix hydrogen bonds

    Output
    - nh_alpha: number of alpha-helix hydrogen bonds
    """

    ## initialise
    nh_alpha=0.0

    ## loop over H-bonds
    for iO,iH in alpha_hbonds:

        ## get separation (need to account for PBCs at some point)
        dr=protein[iO].position-protein[iH].position
        if pbc:
            dr[0]-=boxx[0]*anint(dr[0]/boxx[0])
            dr[1]-=boxx[1]*anint(dr[1]/boxx[1])
            dr[2]-=boxx[2]*anint(dr[2]/boxx[2])
        rr=numpy.sqrt(numpy.dot(dr,dr))
        nh_alpha+=switch(rr,2.5)

    return nh_alpha

## ----------------------------------------------------------------------------

def sum_three_ten_hbonds(protein,three_ten_hbonds,pbc=False,boxx=None):

    """
    sum_three_ten_hbonds

    Function to count number of 3/10-helix hydrogen bonds
    are present in a protein structures

    Count H-bonds using rational switch function

    Input
    - protein: MDanalysis atom selection consisting of protein
    - three_ten_hbonds: list of (O,H) atoms that form alpha-helix hydrogen bonds

    Output
    - nh_three_ten: number of alpha-helix hydrogen bonds
    """

    ## initialise
    nh_three_ten=0.0

    ## loop over H-bonds
    for iO,iH in three_ten_hbonds:

        ## get separation (need to account for PBCs at some point)
        dr=protein[iO].position-protein[iH].position
        if pbc:
            dr[0]-=boxx[0]*anint(dr[0]/boxx[0])
            dr[1]-=boxx[1]*anint(dr[1]/boxx[1])
            dr[2]-=boxx[2]*anint(dr[2]/boxx[2])
        rr=numpy.sqrt(numpy.dot(dr,dr))
        nh_three_ten+=switch(rr,2.5)

    return nh_three_ten


## ----------------------------------------------------------------------------

def get_phi_psi(protein):

    """
    get_phi_psi

    Function to get list of phi and psi angles

    Input
    - protein: MDAnalysis atom selection consisting of protein

    Output
    - phi_angles: list containing phi dihedrals
    - psi_angles: list containing psi dihedrals
    """

    ## initialise lists
    phi_angles=[]
    psi_angles=[]

    nres=len(protein.residues)

    ## loop over residues (omitting first and last)
    for ires in range(1,nres-1):

        ## use MDAnalysis intrinsic to get phi and psi selections
        phi_sel=protein.residues[ires].phi_selection()
        psi_sel=protein.residues[ires].psi_selection()
        i1=phi_sel[0].index
        i2=phi_sel[1].index
        i3=phi_sel[2].index
        i4=phi_sel[3].index
        phi_angles.append(MDAnalysis.core.topologyobjects.Dihedral(numpy.array([i1,i2,i3,i4]),protein))
        i1=psi_sel[0].index
        i2=psi_sel[1].index
        i3=psi_sel[2].index
        i4=psi_sel[3].index
        psi_angles.append(MDAnalysis.core.topologyobjects.Dihedral(numpy.array([i1,i2,i3,i4]),protein))

    return phi_angles,psi_angles

## ----------------------------------------------------------------------------

def beta_dih_offset(phi_angles,psi_angles):

    """
    beta_dih_offset

    Function to calculate dihedral offset compared to beta-strand
    dih_offset=0.5*(sum_phi (1+cos(phi-phi_ref))+ sum_phi (1+cos(psi-psi_ref)))

    Input
    - protein: MDAnalysis atom selection consisting of protein
    - phi_angles: list of phi-dihedral angles
    - psi_angles: list of psi-dihedral angles

    Output
    - dih_offset: value of dihedral offset

    phi_ref=-140 degrees, psi_ref=130 degrees
    """

    dh=0.0
    for p_angle in phi_angles:
        phi=(numpy.pi/180.0)*p_angle.dihedral()
        dh+=(1.0+numpy.cos(phi-beta_phi_ref))

    for p_angle in psi_angles:
        psi=(numpy.pi/180.0)*p_angle.dihedral()
        dh+=(1.0+numpy.cos(psi-beta_psi_ref))

    dh*=0.5
    return dh

## ----------------------------------------------------------------------------

def sum_cgamma_contacts(protein,cgamma_atoms,pbc=False,boxx=None):

    """
    sum_cgamma_contacts

    Function to sum contacts between Cgamma atoms in protein structure
    Uses rational switch function for counting
    Need to add PBCs

    Input
    - protein: MDAnalysis atom selection consisting of protein
    - cgamma_atoms: list of cgamma atoms for each residue

    Output
    - sum_cg: number of Cgamma contacts
    """

    nres=len(cgamma_atoms)

    ## loop over residues
    sum_cg=0.0
    for ires in range(nres):

        ## skip residue if it lacks a CG atom
        if len(cgamma_atoms[ires])==0:
            continue

        ## loop over CG atoms
        for iCG in cgamma_atoms[ires]:

            r1=protein.atoms[iCG-1].position

            ## loop over residues
            for jres in range(ires+1,nres):

                if len(cgamma_atoms[jres])==0:
                    continue

                for jCG in cgamma_atoms[jres]:

                    r2=protein.atoms[jCG-1].position

                    ## get separation
                    dr=r2-r1
                    if pbc:
                        dr[0]-=boxx[0]*anint(dr[0]/boxx[0])
                        dr[1]-=boxx[1]*anint(dr[1]/boxx[1])
                        dr[2]-=boxx[2]*anint(dr[2]/boxx[2])
                    rr=numpy.sqrt(numpy.dot(dr,dr))
                    sum_cg+=switch(rr,rcgamma_cut)

    # ## intialise
    # nres=len(protein.residues)
    # sum_cg=0.0
    #
    # ## double loop over residues
    # for ires in range(nres-1):
    #
    #     ## check residue has a Cgamma atom
    #
    #     try:
    #         icg=protein.residues[ires].atoms.CG.index
    #     except:
    #         continue
    #     for jres in range(ires+1,nres):
    #
    #         ## check residue has Cgamma atom
    #         try:
    #             jcg=protein.residues[jres].atoms.CG.index
    #         except:
    #             continue
    #
    #         ## get separation
    #         dr=protein[icg].position-protein[jcg].position
    #         if pbc:
    #             dr[0]-=boxx[0]*anint(dr[0]/boxx[0])
    #             dr[1]-=boxx[1]*anint(dr[1]/boxx[1])
    #             dr[2]-=boxx[2]*anint(dr[2]/boxx[2])
    #         rr=numpy.sqrt(numpy.dot(dr,dr))
    #         sum_cg+=switch(rr,5.0)
    #
    return sum_cg

## ----------------------------------------------------------------------------

def find_salt_bridges(protein,termini=False):

    """
    find_salt_bridges

    Function to find salt-bridges in proteins
    Salt bridge formed between Lys NZ/Arg CZ and Asp CG/Glu CD

    Input
    - protein: MDAnalysis atom selection consisting of protein
    - termini: include N/C-termini in calculation

    Output
    - zeta_sel: selection consisting of Lys NZ and Arg CZ atoms
    - coo_sel: selection consisting of carboxylic acid carbons
    """

    ## find NZ/CZ atoms
    zeta_sel=protein.select_atoms('resname LYS and name NZ')+protein.select_atoms('resname ARG and name CZ')
    if termini:
        n_sel='resid %4d and name N' % (protein.residues[0].resid)
        zeta_sel+=protein.select_atoms(n_sel)

    ## find carboxylic acid carbons
    coo_sel=protein.select_atoms('resname GLU and name CD')+protein.select_atoms('resname ASP and name CC')
    if termini:
        c_sel='resid %4d and name C' % (protein.residues[-1].resid)

    return zeta_sel,coo_sel

## ----------------------------------------------------------------------------

def calc_gyration_tensor(protein,pbc=False,boxx=None):

    """
    calc_gyration_tensor

    Function to calculate gyration tensor

    Input
    - protein: MDanalysis selection consisting of protein

    Output
    - gyrT: gyration tensor
    - gyrVal: gyration tensor eigenvalues
    - gyrVec: gyration tensor eigenvectors
    """

    ## get centre of mass
    rcom=protein.centroid()

    ## loop over atoms
    gyrT=numpy.zeros((3,3))
    for a in protein.atoms:
        dr=a.position-rcom
        if pbc:
            dr[0]-=boxx[0]*anint(dr[0]/boxx[0])
            dr[1]-=boxx[1]*anint(dr[1]/boxx[1])
            dr[2]-=boxx[2]*anint(dr[2]/boxx[2])
        gyrT+=numpy.outer(dr,dr)
    gyrT/=len(protein.atoms)

    ## diagonalise
    gyrVal,gyrVec=numpy.linalg.eig(gyrT)
    return gyrT,gyrVal,gyrVec

## ----------------------------------------------------------------------------

def sum_contacts(sel_1,sel_2,r0,pbc=False,boxx=None):

    """
    sum_contacts

    Function to get number of contacts between two selections

    Uses rational switch function
    Input
    - protein: MDAnalysis selection consisting of protein
    - sel_1: first selection
    - sel_2: second selection
    - r0: cutoff distance for contact

    Output
    - sum_c: number of contacts
    """

    ## intialise
    sum_c=0.0

    ## loop over selections
    for s_1 in sel_1:

        for s_2 in sel_2:

            dr=s_1.position-s_2.position
            if pbc:
                dr[0]-=boxx[0]*anint(dr[0]/boxx[0])
                dr[1]-=boxx[1]*anint(dr[1]/boxx[1])
                dr[2]-=boxx[2]*anint(dr[2]/boxx[2])

            rr=numpy.sqrt(numpy.dot(dr,dr))
            sum_c+=switch(rr,r0)

    return sum_c

## ----------------------------------------------------------------------------

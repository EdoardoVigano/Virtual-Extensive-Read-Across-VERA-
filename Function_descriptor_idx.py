import pandas as pd
from rdkit import Chem
import numpy as np

# NO rdkit

def f_steroid(molecule: Chem.rdchem.Mol):
    ''' steroid scaffold
        '''
    substructure = '[#6]12~[#6]~[#6]~[#6]3~[#6](~[#6]1~[#6]~[#6]~[#6]2)~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~[#6]34'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)
    
def f_five_ring_hetero(molecule: Chem.rdchem.Mol):
    ''' five_ring_hetero
        '''
    substructure = '[#6]1~[#6]~[#6]~[#6]~[!C;!c]1'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)
    
def f_six_ring_hetero(molecule: Chem.rdchem.Mol):
    ''' six_ring_hetero
        '''
    substructure = '[#6]1~[#6]~[#6]~[#6]~[#6]~[!C;!c]1'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)
    
def f_benzCH2(molecule: Chem.rdchem.Mol):
    ''' benzCH2 present in molecule
        '''
    substructure = 'c1c([C&H2;!R])cccc1' # old smile before 5/09 c1c([C&H2])cccc1
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)
    
def f_heteroatoms_heterocycles6(molecule: Chem.rdchem.Mol):
    ''' heteroatoms (O or S) bonded to a 6 membered
        heterocyclic compound
        '''  
    substructures = ['[O,S]-[#6]1~[#6]~[#6]~[#6]~[#6]~[!C;!c]1',
                    '[#6]1~[#6]([O,S])~[#6]~[#6]~[#6]~[!C;!c]1',
                    '[#6]1~[#6]~[#6]([O,S])~[#6]~[#6]~[!C;!c]1']
    
    substructures = map(Chem.MolFromSmarts, substructures)
    for j, fragment in enumerate(substructures):   
        hit_atss = list(molecule.GetSubstructMatches(fragment))
        hit_bondss = []
        for hit_ats in hit_atss:
            hit_bonds = []
            for bond in fragment.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            hit_bondss.append(hit_bonds)
        list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
        highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
        if j == 0:
            out_atm = list_hit_atss
            out_bonds = highlightBonds
        else:
            out_atm = out_atm + list_hit_atss
            out_bonds = out_bonds + highlightBonds
    return (out_bonds, out_atm)


def f_heteroatoms_heterocycles5(molecule: Chem.rdchem.Mol):

    ''' heteroatoms (O or S) bonded to a 5 membered
        heterocyclic compound
        '''   
    substructures = ['[#6]1~[#6]([O,S])~[#6]~[#6]~[!C;!c]1',
                    '[#6]1~[#6]~[#6]~[#6]([O,S])~[!C;!c]1']
    
    substructures = map(Chem.MolFromSmarts, substructures)
    for j, fragment in enumerate(substructures):   
        hit_atss = list(molecule.GetSubstructMatches(fragment))
        hit_bondss = []
        for hit_ats in hit_atss:
            hit_bonds = []
            for bond in fragment.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            hit_bondss.append(hit_bonds)
        list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
        highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
        if j == 0:
            out_atm = list_hit_atss
            out_bonds = highlightBonds
        else:
            out_atm = out_atm + list_hit_atss
            out_bonds = out_bonds + highlightBonds
    return (out_bonds, out_atm)

def f_biphenol(molecule: Chem.rdchem.Mol):
    
    ''' identify catechol, resorcinol,
        and hydroquinone
    '''   
    substructures = ['[OH]c1c([OH])cccc1',
                    '[OH]c1cc([OH])ccc1',
                    '[OH]c1ccc([OH])cc1']  

    substructures = map(Chem.MolFromSmarts, substructures)
    for j, fragment in enumerate(substructures):   
        hit_atss = list(molecule.GetSubstructMatches(fragment))
        hit_bondss = []
        for hit_ats in hit_atss:
            hit_bonds = []
            for bond in fragment.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            hit_bondss.append(hit_bonds)
        list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
        highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
        if j == 0:
            out_atm = list_hit_atss
            out_bonds = highlightBonds
        else:
            out_atm = out_atm + list_hit_atss
            out_bonds = out_bonds + highlightBonds
    return (out_bonds, out_atm)

def f_benzaldehyde(molecule: Chem.rdchem.Mol):
    ''' benzaldehyde present in molecule
        ''' 
    substructure = 'c1ccccc1-[C;H1](=O)'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

#16/12 add: modify the calculator function
#16/12 add: modify the calculator function
#16/12 add: modify the calculator function

def f_Ar_OR(molecule: Chem.rdchem.Mol):
    ''' aromatic ring with OR sostituent'''
    
    substructure = 'a-[OX2]-[C]-[!O;!N;!S]'
    substructure1 = 'a-[OX2]-[CD3]-[CX4H3]'
    substructure2 = 'a-[OX2]-[CD4C]-[CX4H3]' #carbonio quaternario legato all'O-Ar
    substructure3 = 'a-[OX2]-[CX4H3]'
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    substructure2 = Chem.MolFromSmarts(substructure2) 
    substructure3 = Chem.MolFromSmarts(substructure3) 
    
    sub = [substructure, substructure1, substructure2, 
           substructure3]
           
    for j, i in enumerate(sub):
        hit_atss = list(molecule.GetSubstructMatches(i))
        hit_bondss = []
        for hit_ats in hit_atss:
            hit_bonds = []
            for bond in i.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            hit_bondss.append(hit_bonds)
        list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
        highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
        if j == 0:
            out_atm = list_hit_atss
            out_bonds = highlightBonds
        else:
            out_atm = out_atm + list_hit_atss
            out_bonds = out_bonds + highlightBonds

    return (out_bonds, out_atm)
    
def f_Ar_R(molecule: Chem.rdchem.Mol):
    ''' aromatic ring with R sostituent'''
    
    substructure = 'a-[CX4H2]-[C]-[!O;!N;!S]' 
    substructure1 = 'a-[CD4]-[CX4H3]' #C(CH3)3
    substructure2 = 'a-[CD3]-[CX4H3]' #CH(CH3)2
    substructure3 = 'a-[CX4H3]' #CH(CH3)2
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    
    sub = [substructure, substructure1, substructure2, 
           substructure3]
           
    for j, i in enumerate(sub):
        hit_atss = list(molecule.GetSubstructMatches(i))
        hit_bondss = []
        for hit_ats in hit_atss:
            hit_bonds = []
            for bond in i.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            hit_bondss.append(hit_bonds)
        list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
        highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
        if j == 0:
            out_atm = list_hit_atss
            out_bonds = highlightBonds
        else:
            out_atm = out_atm + list_hit_atss
            out_bonds = out_bonds + highlightBonds

    return (out_bonds, out_atm)
    

def f_op_diphenolo_OR(molecule: Chem.rdchem.Mol):
    
    '''difenoli orto e para e sostituenti OR orto e para'''
    
    substructure = 'c1(-[OX2H])c(-[OX2H])cccc1'
    substructure1 = 'c1(-[OX2H])ccc(-[OX2H])cc1'
    substructure2 = 'c1(-[OX2]-[C;H])c(-[OX2]-[C;H])cccc1'
    substructure3 = 'c1(-[OX2]-[C;H])ccc(-[OX2]-[C;H])cc1'
    substructure4 = 'c1(-[OX2]-[C;H])c(-[OX2H])cccc1'
    substructure5 = 'c1(-[OX2H])ccc(-[OX2]-[C;H])cc1'
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    substructure4 = Chem.MolFromSmarts(substructure4)
    substructure5 = Chem.MolFromSmarts(substructure5)
    
    sub = [substructure, substructure1, substructure2, 
           substructure3, substructure4, substructure5]
           
    for j, i in enumerate(sub):
        hit_atss = list(molecule.GetSubstructMatches(i))
        hit_bondss = []
        for hit_ats in hit_atss:
            hit_bonds = []
            for bond in i.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            hit_bondss.append(hit_bonds)
        list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
        highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
        if j == 0:
            out_atm = list_hit_atss
            out_bonds = highlightBonds
        else:
            out_atm = out_atm + list_hit_atss
            out_bonds = out_bonds + highlightBonds

    return (out_bonds, out_atm)

def f_Ar_COR(molecule: Chem.rdchem.Mol):
    '''Chetone as ring sostituen'''    
    substructure = 'a-[CX3](=O)-[#6]~[!O,!N,!S]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_Ar_COO_R_H(molecule: Chem.rdchem.Mol):
    
    '''sostituente all'anello: estere 
       oppure acido benzoico'''
    
    substructure = 'c-C(=O)[O;H1,-]'
    substructure1 = 'c-C(=O)(-OC-[!O;!N;!S])'
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    
    sub = [substructure, substructure1]
           
    for j, i in enumerate(sub):
        hit_atss = list(molecule.GetSubstructMatches(i))
        hit_bondss = []
        for hit_ats in hit_atss:
            hit_bonds = []
            for bond in i.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            hit_bondss.append(hit_bonds)
        list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
        highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
        if j == 0:
            out_atm = list_hit_atss
            out_bonds = highlightBonds
        else:
            out_atm = out_atm + list_hit_atss
            out_bonds = out_bonds + highlightBonds

    return (out_bonds, out_atm)

def f_C_3phenols(molecule: Chem.rdchem.Mol):
    
    '''3 phenols bound to C'''
    
    substructure = 'C(-a)(-a)(-a)'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)
    
    return indices

def f_Ar_Cl_Br(molecule: Chem.rdchem.Mol):
    
    '''Cl, Br bound to aromatic ring '''
    substructure = 'a[Cl,Br]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)
    
    return indices

def f_Ring_3OH_3OR(molecule: Chem.rdchem.Mol):
    
    '''At least 3 OH or OR in the same ring'''
    
    substructure = 'c1(-[OX2H])cc(-[OX2H])cc(-[OX2H])c1'
    substructure1 = 'c1(-[OX2H])c(-[OX2H])ccc(-[OX2H])c1'
    substructure2 = 'c1(-[OX2H])c(-[OX2H])c(-[OX2H])ccc1'
    substructure3 = 'c1(-[OX2H])c(-[OX2H])c(-[OX2H])c(-[OX2H])cc1'
    substructure4 = 'c1(-[OX2H])c(-[OX2H])c(-[OX2H])cc(-[OX2H])c1'

    substructure5 = 'C1(-[OX2H])~C~C(-[OX2H])~C~C(-[OX2H])~C1'
    substructure6 = 'C1(-[OX2H])~C(-[OX2H])~C~C~C(-[OX2H])~C1'
    substructure7 = 'C1(-[OX2H])~C(-[OX2H])~C(-[OX2H])~C~C~C1'
    substructure8 = 'C1(-[OX2H])~C(-[OX2H])~C(-[OX2H])~C(-[OX2H])~C~C1'
    substructure9 = 'C1(-[OX2H])~C(-[OX2H])~C(-[OX2H])~C~C(-[OX2H])~C1'

    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    substructure4 = Chem.MolFromSmarts(substructure4)
    substructure5 = Chem.MolFromSmarts(substructure5)
    substructure6 = Chem.MolFromSmarts(substructure6)
    substructure7 = Chem.MolFromSmarts(substructure7)
    substructure8 = Chem.MolFromSmarts(substructure8)
    substructure9 = Chem.MolFromSmarts(substructure9)
    
    sub = [substructure, substructure1, substructure2, 
           substructure3, substructure4, substructure5, 
           substructure6, substructure7, substructure8, substructure9]
           
    for j, i in enumerate(sub):
        hit_atss = list(molecule.GetSubstructMatches(i))
        hit_bondss = []
        for hit_ats in hit_atss:
            hit_bonds = []
            for bond in i.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            hit_bondss.append(hit_bonds)
        list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
        highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
        if j == 0:
            out_atm = list_hit_atss
            out_bonds = highlightBonds
        else:
            out_atm = out_atm + list_hit_atss
            out_bonds = out_bonds + highlightBonds

    return (out_bonds, out_atm)

def f_Sulfoxide(molecule: Chem.rdchem.Mol):
    ''' gruppo di sulfoxide '''
    
    substructure = '[#6]-S(=O)-[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_CH2_Terminal(molecule: Chem.rdchem.Mol):
    ''' CH2 terminal '''
    
    substructure = '*~C=[CX3H2]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_Alcool_1(molecule: Chem.rdchem.Mol):
    ''' primary OH'''
    
    substructure = '*-[CH2]-[OX2H]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_triple(molecule: Chem.rdchem.Mol):
    '''triple bond'''
    
    substructure = '*#*'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)


#03/01 sistemare ether rdkit

def f_ethere2(molecule: Chem.rdchem.Mol):
    '''ethere'''
    
    substructure = '[OD2](-[#6]-[#6])-[#6]-[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices0 = molecule.GetSubstructMatches(substructure)

    substructure1= '[OD2](-[CH3])-[#6]-[#6]'
    substructure1= Chem.MolFromSmarts(substructure1)
    sub = [substructure, substructure1]
           
    for j, i in enumerate(sub):
        hit_atss = list(molecule.GetSubstructMatches(i))
        hit_bondss = []
        for hit_ats in hit_atss:
            hit_bonds = []
            for bond in i.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            hit_bondss.append(hit_bonds)
        list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
        highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
        if j == 0:
            out_atm = list_hit_atss
            out_bonds = highlightBonds
        else:
            out_atm = out_atm + list_hit_atss
            out_bonds = out_bonds + highlightBonds

    return (out_bonds, out_atm)

# 10/01

def f_Ar_OSO2(molecule: Chem.rdchem.Mol):
    '''gruppo Ar-OS(=O)2'''
    
    substructure = 'a-O-S(~O)(~O)(*)'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_CC4(molecule: Chem.rdchem.Mol):
    '''gruppo CC4'''
    
    substructure = '[#6D4#6]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

#13/01

def f_aniline_term(molecule: Chem.rdchem.Mol):
    ''' aniline groups
        '''
    substructure = 'c-[NX3H2]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

#FROM RDKIT

def f_benzene(molecule: Chem.rdchem.Mol):
    
    '''benzene '''
    substructure = 'c1ccccc1'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_phenol(molecule: Chem.rdchem.Mol):
    
    '''phenol '''
    substructure = 'c1ccccc1O'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

#Gruppi di RDKit da usare per i duplicati
'''
def f_Ar_OH(molecule: Chem.rdchem.Mol):
    # ArOH
    
    substructure = 'c[OH1]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_benzodiazepine(molecule: Chem.rdchem.Mol):
    #benzodiazepine
    
    substructure = '[c&R2]12[c&R1][c&R1][c&R1][c&R1][c&R2]1[N&R1][C&R1][C&R1][N&R1]=[C&R1]2'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_aryl_methyl(molecule: Chem.rdchem.Mol):
    #aryl_methyl
    
    substructure = '[$(a-[CH3]),$(a-[CH2]-[CH3]),$(a-[CH2]-[CH2]~[!N;!O]);!$(a(:a!:*):a!:*)]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_piperdine(molecule: Chem.rdchem.Mol):
    #piperdine
    substructure = 'N1CCCCC1'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_urea(molecule: Chem.rdchem.Mol):
    # urea
    substructure = 'C(=O)(-N)-N'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_amide(molecule: Chem.rdchem.Mol):
    #amide
    substructure = 'C(=O)-N'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_aniline(molecule: Chem.rdchem.Mol):
    # aniline groups
    substructure = 'c-[NX3;!$(N=*)]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)


def f_NH2(molecule: Chem.rdchem.Mol):
    # "Number of Primary amines"
    substructure = '[NH2,nH2]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_NH1(molecule: Chem.rdchem.Mol):
    # "Number of Secondary amines"
    substructure = '[NH1,nH1]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_NH0(molecule: Chem.rdchem.Mol):
    # "Number of Tertiary amines"
    substructure = '[NH0,nH0]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)


def f_nitro(molecule: Chem.rdchem.Mol):
    # "nitro"

    substructure = 'N(=O)(O)[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_nitro_arom_nonortho(molecule: Chem.rdchem.Mol):
    # nitro_arom_nonortho

    substructure = '[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1);!$(cc-!:*)]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_pyridine(molecule: Chem.rdchem.Mol):
    #"Number of pyridine rings"
    substructure = 'n1ccccc1'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)


def f_priamide(molecule: Chem.rdchem.Mol):
    # Number of primary amides
    substructure = 'C(=O)-[NH2]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_ketone_Topliss(molecule: Chem.rdchem.Mol):
    # "Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha"
    substructure = '[$([CX3](=[OX1])(C)([c,C]));!$([CX3](=[OX1])([CH1]=C)[c,C])]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_ketone(molecule: Chem.rdchem.Mol): 
    # "Number of ketones"
    substructure = '[#6][CX3](=O)[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_priamide(molecule: Chem.rdchem.Mol):
    # Number of primary amides
    substructure = 'C(=O)-[NH2]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_hdrzine(molecule: Chem.rdchem.Mol):
    # "Number of hydrazine groups"
    substructure = '[NX3]-[NX3]' #io farei [NX3H]-[NX3H2]
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_ester(molecule: Chem.rdchem.Mol):
    # "Number of esters"
    substructure = '[#6][CX3](=O)[OX2H0][#6]' 
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_ether(molecule: Chem.rdchem.Mol):
    # "Number of ether oxygens (including phenoxy)"
    substructure = '[OD2]([#6])[#6]' 
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)


def f_imidazole(molecule: Chem.rdchem.Mol):
    
    "Number of imidazole rings"
    substructure = 'n1cncc1'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)'''
# aggiunte 22/07

def f_Imidothioesters(molecule: Chem.rdchem.Mol):
    "Number of Imidothioesters"
    substructure = '[#1,#6][CX3](=[NX2][#1,#6])[SX2][C&!$([CX3]=[OX1,SX1,NX2])]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_anhydrides(molecule: Chem.rdchem.Mol):
    '''anhydrides'''
    
    substructure = '[#6][#6X3](=[OX1])[#8X2][#6X3](=[OX1])[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_alkyl_halide(molecule: Chem.rdchem.Mol):
    '''Number of aniline'''
    
    substructure = '[C]-[Cl,Br,I,F]' 
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_nitro_aliphatic(molecule: Chem.rdchem.Mol):
    '''nitro aliphatic'''
    substructure = 'N(=O)(O)-[C]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_ketone_aliphatic(molecule: Chem.rdchem.Mol):
    '''ketone aliphatic'''
    substructure = '[C]-[CX3](=O)-[C]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_Ar_ketone(molecule: Chem.rdchem.Mol):
    '''ketone aromatic'''
    substructure = '[c]-[CX3](=O)-[C]'
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

def f_carbamate(molecule: Chem.rdchem.Mol):
    #carbamate
    substructure = '[$(N([!H])C(=O)[OH0])]C(=O)[OH0]C' 
    substructure = Chem.MolFromSmarts(substructure)
    hit_atss = list(molecule.GetSubstructMatches(substructure))
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in substructure.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]    
    return (highlightBonds, list_hit_atss)

# rdkit groups
def f_C_O(molecule: Chem.rdchem.Mol):
	'''Number of carbonyl O'''
	substructure = '[CX3]=[OX1]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_C_O_noCOO(molecule: Chem.rdchem.Mol):
	'''Number of carbonyl O, excluding COOH'''
	substructure = '[C!$(C-[OH])]=O'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Al_OH(molecule: Chem.rdchem.Mol):
	'''Number of aliphatic hydroxyl groups'''
	substructure = '[C!$(C=O)]-[OH]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Ar_OH(molecule: Chem.rdchem.Mol):
	'''Number of aromatic hydroxyl groups'''
	substructure = 'c[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_methoxy(molecule: Chem.rdchem.Mol):
	'''Number of methoxy groups -OCH3'''
	substructure = '[OX2](-[#6])-[CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_oxime(molecule: Chem.rdchem.Mol):
	'''Number of oxime groups'''
	substructure = '[CX3]=[NX2]-[OX2]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_ester(molecule: Chem.rdchem.Mol):
	'''Number of esters'''
	substructure = '[#6][CX3](=O)[OX2H0][#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Al_COO(molecule: Chem.rdchem.Mol):
	'''Number of aliphatic carboxylic acids'''
	substructure = 'C-C(=O)[O;H1,-]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Ar_COO(molecule: Chem.rdchem.Mol):
	'''Number of Aromatic carboxylic acide'''
	substructure = 'c-C(=O)[O;H1,-]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_COO(molecule: Chem.rdchem.Mol):
	'''Number of carboxylic acids'''
	substructure = '[#6]C(=O)[O;H,-1]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_COO2(molecule: Chem.rdchem.Mol):
	'''Number of carboxylic acids'''
	substructure = '[CX3](=O)[OX1H0-,OX2H1]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_ketone(molecule: Chem.rdchem.Mol):
	'''Number of ketones'''
	substructure = '[#6][CX3](=O)[#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_ether(molecule: Chem.rdchem.Mol):
	'''Number of ether oxygens (including phenoxy)'''
	substructure = '[OD2]([#6])[#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_phenol(molecule: Chem.rdchem.Mol):
	'''Number of phenols'''
	substructure = '[OX2H]-c1ccccc1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_aldehyde(molecule: Chem.rdchem.Mol):
	'''Number of aldehydes'''
	substructure = '[CX3H1](=O)[#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_quatN(molecule: Chem.rdchem.Mol):
	'''Number of quarternary nitrogens'''
	substructure = '[$([NX4+]),$([NX4]=*)]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_NH2(molecule: Chem.rdchem.Mol):
	'''Number of Primary amines'''
	substructure = '[NH2,nH2]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_NH1(molecule: Chem.rdchem.Mol):
	'''Number of Secondary amines'''
	substructure = '[NH1,nH1]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_NH0(molecule: Chem.rdchem.Mol):
	'''Number of Tertiary amines'''
	substructure = '[NH0,nH0]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Ar_N(molecule: Chem.rdchem.Mol):
	'''Number of aromatic nitrogens'''
	substructure = 'n'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Ar_NH(molecule: Chem.rdchem.Mol):
	'''Number of aromatic amines'''
	substructure = '[nH]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_aniline(molecule: Chem.rdchem.Mol):
	'''Number of anilines'''
	substructure = 'c-[NX3;!$(N=*)]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Imine(molecule: Chem.rdchem.Mol):
	'''Number of Imines'''
	substructure = '[Nv3](=C)-[#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_nitrile(molecule: Chem.rdchem.Mol):
	'''Number of nitriles'''
	substructure = '[NX1]#[CX2]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_hdrzine(molecule: Chem.rdchem.Mol):
	'''Number of hydrazine groups'''
	substructure = '[NX3]-[NX3]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_hdrzone(molecule: Chem.rdchem.Mol):
	'''Number of hydrazone groups'''
	substructure = 'C=N-[NX3]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_nitroso(molecule: Chem.rdchem.Mol):
	'''Number of nitroso groups, excluding NO2'''
	substructure = '[N!$(N-O)]=O'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_N_O(molecule: Chem.rdchem.Mol):
	'''Number of hydroxylamine groups'''
	substructure = '[N!$(N=O)](-O)-C'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_nitro(molecule: Chem.rdchem.Mol):
	'''Number of nitro groups'''
	substructure = '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_azo(molecule: Chem.rdchem.Mol):
	'''Number of azo groups'''
	substructure = '[#6]-N=N-[#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_diazo(molecule: Chem.rdchem.Mol):
	'''Number of diazo groups'''
	substructure = '[N+]#N'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_azide(molecule: Chem.rdchem.Mol):
	'''Number of azide groups'''
	substructure = '[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_amide(molecule: Chem.rdchem.Mol):
	'''Number of amides'''
	substructure = 'C(=O)-N'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_priamide(molecule: Chem.rdchem.Mol):
	'''Number of primary amides'''
	substructure = 'C(=O)-[NH2]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_amidine(molecule: Chem.rdchem.Mol):
	'''Number of amidine groups'''
	substructure = 'C(=N)(-N)-[!#7]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_guanido(molecule: Chem.rdchem.Mol):
	'''Number of guanidine groups'''
	substructure = 'C(=N)(N)N'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Nhpyrrole(molecule: Chem.rdchem.Mol):
	'''Number of H-pyrrole nitrogens'''
	substructure = '[nH]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_imide(molecule: Chem.rdchem.Mol):
	'''Number of imide groups'''
	substructure = 'N(-C(=O))-C=O'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_isocyan(molecule: Chem.rdchem.Mol):
	'''Number of isocyanates'''
	substructure = 'N=C=O'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_isothiocyan(molecule: Chem.rdchem.Mol):
	'''Number of isothiocyanates'''
	substructure = 'N=C=S'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_thiocyan(molecule: Chem.rdchem.Mol):
	'''Number of thiocyanates'''
	substructure = 'S-C#N'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_halogen(molecule: Chem.rdchem.Mol):
	'''Number of halogens'''
	substructure = '[#9,#17,#35,#53]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_alkyl_halide(molecule: Chem.rdchem.Mol):
	'''Number of alkyl halides'''
	substructure = '[CX4]-[Cl,Br,I,F]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_sulfide(molecule: Chem.rdchem.Mol):
	'''Number of thioether'''
	substructure = '[SX2](-[#6])-C'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_SH(molecule: Chem.rdchem.Mol):
	'''Number of thiol groups'''
	substructure = '[SH]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_C_S(molecule: Chem.rdchem.Mol):
	'''Number of thiocarbonyl'''
	substructure = 'C=[SX1]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_sulfone(molecule: Chem.rdchem.Mol):
	'''Number of sulfone groups'''
	substructure = 'S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_sulfone(molecule: Chem.rdchem.Mol):
	'''Number of sulfone groups'''
	substructure = 'S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_sulfonamd(molecule: Chem.rdchem.Mol):
	'''Number of sulfonamides'''
	substructure = 'N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_prisulfonamd(molecule: Chem.rdchem.Mol):
	'''Number of primary sulfonamides'''
	substructure = '[NH2]-S(=,-[OX1;+0;-1])(=,-[OX1;+0;-1])-[#6]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_barbitur(molecule: Chem.rdchem.Mol):
	'''Number of barbiturate groups'''
	substructure = 'C1C(=O)NC(=O)NC1=O'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_urea(molecule: Chem.rdchem.Mol):
	'''Number of urea groups'''
	substructure = 'C(=O)(-N)-N'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_term_acetylene(molecule: Chem.rdchem.Mol):
	'''Number of terminal acetylenes'''
	substructure = 'C#[CH]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_imidazole(molecule: Chem.rdchem.Mol):
	'''Number of imidazole rings'''
	substructure = 'n1cncc1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_furan(molecule: Chem.rdchem.Mol):
	'''Number of furan rings'''
	substructure = 'o1cccc1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_thiophene(molecule: Chem.rdchem.Mol):
	'''Number of thiophene rings'''
	substructure = 's1cccc1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_thiazole(molecule: Chem.rdchem.Mol):
	'''Number of thiazole rings'''
	substructure = 'c1scnc1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_oxazole(molecule: Chem.rdchem.Mol):
	'''Number of oxazole rings'''
	substructure = 'c1ocnc1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_pyridine(molecule: Chem.rdchem.Mol):
	'''Number of pyridine rings'''
	substructure = 'n1ccccc1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_piperdine(molecule: Chem.rdchem.Mol):
	'''Number of piperdine rings'''
	substructure = 'N1CCCCC1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_piperzine(molecule: Chem.rdchem.Mol):
	'''Number of piperzine rings'''
	substructure = 'N1CCNCC1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_morpholine(molecule: Chem.rdchem.Mol):
	'''Number of morpholine rings'''
	substructure = 'O1CCNCC1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_lactam(molecule: Chem.rdchem.Mol):
	'''Number of beta lactams'''
	substructure = 'N1C(=O)CC1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_lactone(molecule: Chem.rdchem.Mol):
	'''Number of cyclic esters (lactones)'''
	substructure = '[C&R1](=O)[O&R1][C&R1]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_tetrazole(molecule: Chem.rdchem.Mol):
	'''Number of tetrazole rings'''
	substructure = 'c1nnnn1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_epoxide(molecule: Chem.rdchem.Mol):
	'''Number of epoxide rings'''
	substructure = 'O1CC1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_unbrch_alkane(molecule: Chem.rdchem.Mol):
	'''Number of unbranched alkanes  of at least 4 members (excludes halogenated alkanes)'''
	substructure = '[R0;D2][R0;D2][R0;D2][R0;D2]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_bicyclic(molecule: Chem.rdchem.Mol):
	'''Bicyclic'''
	substructure = '[R2][R2]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_benzene(molecule: Chem.rdchem.Mol):
	'''Number of benzene rings'''
	substructure = 'c1ccccc1'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_phos_acid(molecule: Chem.rdchem.Mol):
	'''Number of phosphoric acid groups'''
	substructure = '[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_phos_ester(molecule: Chem.rdchem.Mol):
	'''Number of phosphoric ester groups'''
	substructure = '[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_nitro_arom(molecule: Chem.rdchem.Mol):
	'''Number of nitro benzene ring substituents'''
	substructure = '[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1)]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_nitro_arom_nonortho(molecule: Chem.rdchem.Mol):
	'''Number of non-ortho nitro benzene ring substituents'''
	substructure = '[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1);!$(cc-!:*)]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_dihydropyridine(molecule: Chem.rdchem.Mol):
	'''Number of dihydropyridines'''
	substructure = '[$([NX3H1]1-C=C-C-C=C1),$([Nv3]1=C-C-C=C-C1),$([Nv3]1=C-C=C-C-C1),$([NX3H1]1-C-C=C-C=C1)]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_phenol_noOrthoHbond(molecule: Chem.rdchem.Mol):
	'''Number of phenolic OH excluding ortho intramolecular Hbond substituents'''
	substructure = '[$(c1(-[OX2H])ccccc1);!$(cc-!:[CH2]-[OX2H]);!$(cc-!:C(=O)[O;H1,-]);!$(cc-!:C(=O)-[NH2])]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Al_OH_noTert(molecule: Chem.rdchem.Mol):
	'''Number of aliphatic hydroxyl groups excluding tert-OH'''
	substructure = '[$(C-[OX2H]);!$([CX3](-[OX2H])=[OX1]);!$([CD4]-[OX2H])]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_benzodiazepine(molecule: Chem.rdchem.Mol):
	'''Number of benzodiazepines with no additional fused rings'''
	substructure = '[c&R2]12[c&R1][c&R1][c&R1][c&R1][c&R2]1[N&R1][C&R1][C&R1][N&R1]=[C&R1]2'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_para_hydroxylation(molecule: Chem.rdchem.Mol):
	'''Number of para-hydroxylation sites'''
	substructure = '[$([cH]1[cH]cc(c[cH]1)~[$([#8,$([#8]~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)~[$([#7X3,$([#7](~[H,c,C])~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)-!:[$([NX3H,$(NC(=O)[H,c,C])])])]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_allylic_oxid(molecule: Chem.rdchem.Mol):
	'''Number of allylic oxidation sites excluding steroid dienone'''
	substructure = '[$(C=C-C);!$(C=C-C-[N,O,S]);!$(C=C-C-C-[N,O]);!$(C12=CC(=O)CCC1C3C(C4C(CCC4)CC3)CC2)]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_aryl_methyl(molecule: Chem.rdchem.Mol):
	'''Number of aryl methyl sites for hydroxylation'''
	substructure = '[$(a-[CH3]),$(a-[CH2]-[CH3]),$(a-[CH2]-[CH2]~[!N;!O]);!$(a(:a!:*):a!:*)]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Ndealkylation1(molecule: Chem.rdchem.Mol):
	'''Number of XCCNR groups'''
	substructure = '[$(N(-[CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)]),$(N(-[CH2][CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)])]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Ndealkylation2(molecule: Chem.rdchem.Mol):
	'''Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)'''
	substructure = '[$([N&R1]1(-C)CCC1),$([N&R1]1(-C)CCCC1),$([N&R1]1(-C)CCCCC1),$([N&R1]1(-C)CCCCCC1),$([N&R1]1(-C)CCCCCCC1)]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_alkyl_carbamate(molecule: Chem.rdchem.Mol):
	'''Number of alkyl carbamates (subject to hydrolysis)'''
	substructure = 'C[NH1]C(=O)OC'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_ketone_Topliss(molecule: Chem.rdchem.Mol):
	'''Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha'''
	substructure = '[$([CX3](=[OX1])(C)([c,C]));!$([CX3](=[OX1])([CH1]=C)[c,C])]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_ArN(molecule: Chem.rdchem.Mol):
	'''Number of N functional groups attached to aromatics'''
	substructure = '[$(a-[NX3H2]),$(a-[NH1][NH2]),$(a-C(=[OX1])[NH1][NH2]),$(a-C(=[NH])[NH2])]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_HOCCN(molecule: Chem.rdchem.Mol):
	'''Number of C(OH)CCN-Ctert-alkyl or  C(OH)CCNcyclic'''
	substructure = '[$([OX2H1][CX4][CX4H2][NX3&R1]),$([OH1][CX4][CX4H2][NX3][CX4](C)(C)C)]'
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)



def f_ketone_dehydro(molecule: Chem.rdchem.Mol):
	'''ketone that can generate deihydro product with OH in beta e H in alpha'''	
	substructure = ('[C]-[Cv2](=O)-[Cv3h]-[Cv3h](-[Oh])[*]')
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Ar_bicycle(molecule: Chem.rdchem.Mol):
	'''Aromatic bicyles'''	
	substructure = ('[R2;a][R2;a]')
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)

def f_Al_bicycle(molecule: Chem.rdchem.Mol):
	'''Alifatic bicyles'''	
	substructure = ('[R2;!a][R2;!a]')
	substructure = Chem.MolFromSmarts(substructure)
	hit_atss = list(molecule.GetSubstructMatches(substructure))
	hit_bondss = []
	for hit_ats in hit_atss:
		hit_bonds = []
		for bond in substructure.GetBonds():
			aid1 = hit_ats[bond.GetBeginAtomIdx()]
			aid2 = hit_ats[bond.GetEndAtomIdx()]
			hit_bonds.append(molecule.GetBondBetweenAtoms(aid1,aid2).GetIdx())
			hit_bondss.append(hit_bonds)
	list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
	highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
	return (highlightBonds, list_hit_atss)





def inx_groups(mol: Chem.rdchem.Mol):

    Al_COO=list(); Al_OH=list(); Al_OH_noTert=list(); ArN=list(); Ar_N=list(); Ar_NH=list()
    Ar_OH=list(); COO2=list(); C_O=list(); C_O_noCOO=list(); C_S=list(); HOCCN=list(); Imine=list(); NH0=list(); NH1=list(); NH2=list()
    N_O=list(); Ndealkylation1=list(); Ndealkylation2=list(); Nhpyrrole=list(); SH=list(); aldehyde=list(); alkyl_carbamate=list(); 
    alkyl_halide=list(); allylic_oxid=list(); amide=list(); amidine=list(); aniline=list(); aryl_methyl=list(); azide=list(); 
    azo=list(); barbitur=list(); benzene=list(); benzodiazepine=list(); bicyclic=list(); diazo=list(); 
    dihydropyridine=list(); epoxide=list(); ester=list(); ether=list(); furan=list(); guanido=list(); halogen=list(); hdrzine=list(); 
    hdrzone=list(); imidazole=list(); imide=list(); isocyan=list(); isothiocyan=list(); 
    ketone=list(); ketone_Topliss=list(); lactam=list(); lactone=list(); methoxy=list(); morpholine=list(); nitrile=list(); 
    nitro=list(); nitro_arom=list(); nitro_arom_nonortho=list(); nitroso=list(); oxazole=list()
    oxime=list(); para_hydroxylation=list(); phenol=list(); phenol_noOrthoHbond=list(); phos_acid=list(); phos_ester=list()
    piperdine=list(); piperzine=list(); priamide=list(); prisulfonamd=list(); pyridine=list(); 
    quatN=list(); sulfide=list(); sulfonamd=list(); sulfone=list(); term_acetylene=list()
    tetrazole=list(); thiazole=list(); thiocyan=list(); thiophene=list(); unbrch_alkane=list(); urea=list(); 
    
    #ADD:
    heteroatoms_heterocycles6_list = []
    heteroatoms_heterocycles5_list = []
    steroid_list = []
    hetero5_list = []
    hetero6_list = []
    benzCH2_list = []
    benzaldehyde_list = []
    biphenol_list = []
    # 17/12
    Ar_OR_list = []
    Ar_R_list = []
    op_diphenolo_OR_list = []
    Ar_COR_list = []
    Ar_COO_R_H_list = []
    C_3phenols_list = []
    Ar_Cl_Br_list = []
    Ring_3OH_3OR_list = []
    Sulfoxide_list = []
    CH2_Terminal_list = []
    Alcool_1_list = []
    
    #22/12
    triple_list = []
    
    # 10/01 dimenticato?
    Ar_COO = []
    
    # 10/01
    Ar_OSO2 = []
    CC4 = []
    
    #erika
    
    Imidothioesters = []
    anhydrides = []
    carbamate = []
    
    #13/01
    aniline_term = []
    
    #25/01 modify fr_alkyl_halide
    
    #26/01
    charge2 = []
    charge1 = []
    nitro_aliphatic = []
    ketone_aliphatic = []
    Ar_ketone = []
	# 09/09
    ketone_dehydro = []
    
    #23 Nov 23 - Ar and Al bicycle
    Ar_bicycle = []
    Al_bicycle = []


    lst_chemical_groups=[Al_COO,Al_OH,Al_OH_noTert,ArN,Ar_N,Ar_NH,Ar_OH,COO2,C_O,C_O_noCOO,C_S,HOCCN,Imine,NH0,NH1,NH2,N_O,
                         Ndealkylation1,Ndealkylation2,Nhpyrrole,SH,aldehyde,alkyl_carbamate,alkyl_halide,allylic_oxid,amide,
                         amidine,aniline,aryl_methyl,azide,azo,barbitur,benzene,benzodiazepine,bicyclic,diazo,dihydropyridine,
                         epoxide,ester,ether,furan,guanido,halogen,hdrzine,hdrzone,imidazole,imide,isocyan,isothiocyan,ketone,
                         ketone_Topliss,lactam,lactone,methoxy,morpholine,nitrile,nitro,
                         nitro_arom,nitro_arom_nonortho,nitroso,oxazole,oxime,para_hydroxylation,phenol,
                         phenol_noOrthoHbond,phos_acid,phos_ester,piperdine,piperzine,priamide,
                         prisulfonamd,pyridine,quatN,sulfide,sulfonamd,sulfone,term_acetylene,tetrazole,thiazole,thiocyan, thiophene, unbrch_alkane,urea,
						 heteroatoms_heterocycles5_list,heteroatoms_heterocycles6_list,steroid_list,hetero5_list,hetero6_list,benzCH2_list,benzaldehyde_list, 
						 biphenol_list, Ar_OR_list, Ar_R_list, op_diphenolo_OR_list, Ar_COR_list, Ar_COO_R_H_list, C_3phenols_list, Ar_Cl_Br_list, 
						 Ring_3OH_3OR_list, Sulfoxide_list, CH2_Terminal_list, Alcool_1_list, triple_list, Ar_COO, Ar_OSO2, CC4, Imidothioesters, 
						 anhydrides, carbamate, aniline_term, charge2, charge1, nitro_aliphatic, ketone_aliphatic, Ar_ketone, ketone_dehydro, Ar_bicycle, Al_bicycle]
    
    if mol != 'none':    
        Al_COO.append(f_Al_COO(mol))
        Al_OH.append(f_Al_OH(mol))
        Al_OH_noTert.append(f_Al_OH_noTert(mol))
        ArN.append(f_ArN(mol))
        Ar_N.append(f_Ar_N(mol))
        Ar_NH.append(f_Ar_NH(mol))
        Ar_OH.append(f_Ar_OH(mol))
        COO2.append(f_COO2(mol))
        C_O.append(f_C_O(mol))
        C_O_noCOO.append(f_C_O_noCOO(mol))
        C_S.append(f_C_S(mol))
        HOCCN.append(f_HOCCN(mol))
        Imine.append(f_Imine(mol))
        NH0.append(f_NH0(mol))
        NH1.append(f_NH1(mol))
        NH2.append(f_NH2(mol))
        N_O.append(f_N_O(mol))
        Ndealkylation1.append(f_Ndealkylation1(mol))
        Ndealkylation2.append(f_Ndealkylation2(mol))
        Nhpyrrole.append(f_Nhpyrrole(mol))
        SH.append(f_SH(mol))
        aldehyde.append(f_aldehyde(mol))
        alkyl_carbamate.append(f_alkyl_carbamate(mol))
        #alkyl_halide.append(f_alkyl_halide(mol))
        allylic_oxid.append(f_allylic_oxid(mol))
        amide.append(f_amide(mol))
        amidine.append(f_amidine(mol))
        aniline.append(f_aniline(mol))
        aryl_methyl.append(f_aryl_methyl(mol))
        azide.append(f_azide(mol))
        azo.append(f_azo(mol))
        barbitur.append(f_barbitur(mol))
        benzene.append(f_benzene(mol))
        benzodiazepine.append(f_benzodiazepine(mol))
        bicyclic.append(f_bicyclic(mol))
        diazo.append(f_diazo(mol))
        dihydropyridine.append(f_dihydropyridine(mol))
        epoxide.append(f_epoxide(mol))
        ester.append(f_ester(mol))
     
        furan.append(f_furan(mol))
        guanido.append(f_guanido(mol))
        halogen.append(f_halogen(mol))
        hdrzine.append(f_hdrzine(mol))
        hdrzone.append(f_hdrzone(mol))
        imidazole.append(f_imidazole(mol))
        imide.append(f_imide(mol))
        isocyan.append(f_isocyan(mol))
        isothiocyan.append(f_isothiocyan(mol))
        ketone.append(f_ketone(mol))
        ketone_Topliss.append(f_ketone_Topliss(mol))
        lactam.append(f_lactam(mol))
        lactone.append(f_lactone(mol))
        methoxy.append(f_methoxy(mol))
        morpholine.append(f_morpholine(mol))
        nitrile.append(f_nitrile(mol))
        nitro.append(f_nitro(mol))
        nitro_arom.append(f_nitro_arom(mol))
        nitro_arom_nonortho.append(f_nitro_arom_nonortho(mol))
        nitroso.append(f_nitroso(mol))
        oxazole.append(f_oxazole(mol))
        oxime.append(f_oxime(mol))
        para_hydroxylation.append(f_para_hydroxylation(mol))
        phenol.append(f_phenol(mol))
        phenol_noOrthoHbond.append(f_phenol_noOrthoHbond(mol))
        phos_acid.append(f_phos_acid(mol))
        phos_ester.append(f_phos_ester(mol))
        piperdine.append(f_piperdine(mol))
        piperzine.append(f_piperzine(mol))
        priamide.append(f_priamide(mol))
        prisulfonamd.append(f_prisulfonamd(mol))
        pyridine.append(f_pyridine(mol))
        quatN.append(f_quatN(mol))
        sulfide.append(f_sulfide(mol))
        sulfonamd.append(f_sulfonamd(mol))
        sulfone.append(f_sulfone(mol))
        term_acetylene.append(f_term_acetylene(mol))
        tetrazole.append(f_tetrazole(mol))
        thiazole.append(f_thiazole(mol))
        thiocyan.append(f_thiocyan(mol))
        thiophene.append(f_thiophene(mol))
        unbrch_alkane.append(f_unbrch_alkane(mol))
        urea.append(f_urea(mol))
        
        heteroatoms_heterocycles6_list.append(f_heteroatoms_heterocycles6(mol))
        heteroatoms_heterocycles5_list.append(f_heteroatoms_heterocycles5(mol))
        steroid_list.append(f_steroid(mol))
        hetero5_list.append(f_five_ring_hetero(mol))
        hetero6_list.append(f_six_ring_hetero(mol))
        benzCH2_list.append(f_benzCH2(mol))
        benzaldehyde_list.append(f_benzaldehyde(mol))
        biphenol_list.append(f_biphenol(mol))
        
        Ar_OR_list.append(f_Ar_OR(mol))
        Ar_R_list.append(f_Ar_R(mol))
        op_diphenolo_OR_list.append(f_op_diphenolo_OR(mol))
        Ar_COR_list.append(f_Ar_COR(mol))
        Ar_COO_R_H_list.append(f_Ar_COO_R_H(mol))
        C_3phenols_list.append(f_C_3phenols(mol))
        Ar_Cl_Br_list.append(f_Ar_Cl_Br(mol))
        Ring_3OH_3OR_list.append(f_Ring_3OH_3OR(mol))
        Sulfoxide_list.append(f_Sulfoxide(mol))
        CH2_Terminal_list.append(f_CH2_Terminal(mol))
        Alcool_1_list.append(f_Alcool_1(mol))
        
        triple_list.append(f_triple(mol))
        
        ether.append(f_ethere2(mol)) #sostituisce ethere di rdkit
        
        #ti eri scordato Ar_COO?
        Ar_COO.append(f_Ar_COO(mol))
        # nuovo 10/1
        Ar_OSO2.append(f_Ar_OSO2(mol))
        CC4.append(f_CC4(mol))
        
        #erika
        Imidothioesters.append(f_Imidothioesters(mol))
        anhydrides.append(f_anhydrides(mol))
        carbamate.append(f_carbamate(mol))
        
        #13/01
        aniline_term.append(f_aniline_term(mol))
        
        #25/01 rdkit alkyl_halide have to be replaced
        alkyl_halide.append(f_alkyl_halide(mol))
        
        #26/01
        nitro_aliphatic.append(f_nitro_aliphatic(mol))
        ketone_aliphatic.append(f_ketone_aliphatic(mol))
        Ar_ketone.append(f_Ar_ketone(mol))
        # 09/09
        ketone_dehydro.append(f_ketone_dehydro(mol))
        #23 Nov 23
        Ar_bicycle.append(f_Ar_bicycle(mol))
        Al_bicycle.append(f_Al_bicycle(mol))
    else:
        for chem_group in lst_chemical_groups:
            chem_group.append('none')
    
    
    lista1 = [Al_COO, Al_OH, Al_OH_noTert, ArN, Ar_N, Ar_NH, Ar_OH,COO2, C_O, C_O_noCOO, C_S,
            HOCCN, Imine, NH0, NH1, NH2, N_O, Ndealkylation1, Ndealkylation2, Nhpyrrole, SH,
            aldehyde, alkyl_carbamate, alkyl_halide, allylic_oxid, amide,amidine, aniline, aryl_methyl,
            azide, azo, barbitur, benzene, benzodiazepine, bicyclic, diazo, dihydropyridine, epoxide,
            ester,ether, furan, guanido, halogen, hdrzine, hdrzone, imidazole, imide, isocyan, isothiocyan,
            ketone, ketone_Topliss, lactam, lactone, methoxy, morpholine, nitrile, nitro, 
            nitro_arom, nitro_arom_nonortho, nitroso, oxazole, oxime, para_hydroxylation, phenol, phenol_noOrthoHbond, 
            phos_acid, phos_ester, piperdine, piperzine, priamide, prisulfonamd, pyridine, quatN, 
            sulfide, sulfonamd, sulfone, term_acetylene, tetrazole, thiazole, thiocyan, thiophene,  unbrch_alkane, urea, 
            heteroatoms_heterocycles5_list,  heteroatoms_heterocycles6_list, steroid_list, hetero5_list, 
            hetero6_list, benzCH2_list, benzaldehyde_list, biphenol_list, Ar_OR_list, Ar_R_list, op_diphenolo_OR_list,
            Ar_COR_list, Ar_COO_R_H_list, C_3phenols_list, Ar_Cl_Br_list, Ring_3OH_3OR_list, Sulfoxide_list, 
            CH2_Terminal_list, Alcool_1_list, triple_list,Ar_COO, Ar_OSO2, CC4, Imidothioesters, 
            anhydrides,carbamate,aniline_term,charge2,charge1,nitro_aliphatic,ketone_aliphatic,Ar_ketone,ketone_dehydro,
			Ar_bicycle,Al_bicycle]
    
    lista2 = ["Al_COO","Al_OH","Al_OH_noTert","ArN","Ar_N","Ar_NH","Ar_OH","COO2","C_O",
            "C_O_noCOO","C_S","HOCCN","Imine","NH0","NH1","NH2","N_O","Ndealkylation1",
            "Ndealkylation2","Nhpyrrole","SH","aldehyde","alkyl_carbamate","alkyl_halide","allylic_oxid",
            "amide","amidine","aniline","aryl_methyl","azide","azo","barbitur","benzene","benzodiazepine",
            "bicyclic","diazo","dihydropyridine","epoxide","ester","ether","furan","guanido","halogen","hdrzine",
            "hdrzone","imidazole","imide","isocyan","isothiocyan","ketone","ketone_Topliss","lactam",
            "lactone","methoxy","morpholine","nitrile","nitro","nitro_arom","nitro_arom_nonortho","nitroso",
            "oxazole","oxime","para_hydroxylation","phenol","phenol_noOrthoHbond","phos_acid",
            "phos_ester","piperdine","piperzine","priamide","prisulfonamd","pyridine","quatN",
            "sulfide","sulfonamd","sulfone","term_acetylene","tetrazole","thiazole","thiocyan",
            "thiophene","unbrch_alkane","urea","heteroatoms_heterocycles5",
            "heteroatoms_heterocycles6","steroid","hetero5","hetero6","benzCH2",
            "benzaldehyde","biphenol", "Ar_OR", "Ar_R", "op_diphenolo_OR", "Ar_COR", 
            "Ar_COO_R_H", "C_3phenols", "Ar_Cl_Br", "Ring_3OH_3OR", "Sulfoxide", "CH2_Terminal", 
            "Alcool_1", "triple_bond", "Ar_COO", "Ar_OSO2", "CC4", "Imidothioesters","anhydrides",
            "carbamate","aniline_term","charge+","charge-","nitro_aliphatic",'ketone_aliphatic',
            'Ar_ketone','ketone_dehydro','Ar_bicycle','Al_bicycle']

    out = {lista2[i]:lista1[i] for i in range(len(lista1))}
    return out


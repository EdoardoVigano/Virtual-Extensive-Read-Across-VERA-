import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments
import numpy as np

'''Note:
Every time you add something here you have to add the functions in file.py alerts_index, calculator and alerts
'''

'''
Note for add Alerts:


    template for alerts with only one smart:

    def f_SA1(molecule: Chem.rdchem.Mol):
        substructure = 'SMARTS'
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

    template for alerts with more than one smart:

    def f_SA4(molecule: Chem.rdchem.Mol):

        substructure = 'SMART'
        substructure1 = 'SMART1
        
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

After you have to add this count to functions:
    def all_index_atoms(molecule: Chem.rdchem.Mol):
    def all_index_bonds(molecule: Chem.rdchem.Mol):

    at the end of the code!
'''

def fish1(molecule: Chem.rdchem.Mol):
    ''' FISH 1
        '''
    substructure = 'O(c1ccccc1)c1ccccc1'
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
    
def fish2(molecule: Chem.rdchem.Mol):
    ''' FISH 2
        '''
    substructure = '[O-][N+](=O)c1ccccc1'
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

def fish3(molecule: Chem.rdchem.Mol):
    ''' FISH 3
        '''
    substructure = 'C[S]c1ccccc1'
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

def fish4(molecule: Chem.rdchem.Mol):
    ''' FISH 4
        '''
    substructure = '*O[CX4H2]C#N'
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

def fish5(molecule: Chem.rdchem.Mol):
    ''' FISH 5
        '''
    substructure = '[*;A][N+]([*;A])([*;A])[*;A]'
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

def fish6(molecule: Chem.rdchem.Mol):
    ''' FISH 6
        '''
    substructure = '*SS*'
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

def fish7(molecule: Chem.rdchem.Mol):
    ''' FISH 7
        '''
    substructure = 'Cl[C]=[C;A]'
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

def fish8(molecule: Chem.rdchem.Mol):
    ''' FISH 8
        '''
    substructure = 'Cn1cc(C(N)=O)c(C)n1'
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

def fish9(molecule: Chem.rdchem.Mol):
    ''' FISH 9
        '''
    substructure = '*P(*)(*)=*'
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

def fish10(molecule: Chem.rdchem.Mol): #Modifiche dell'8/09/22 fino a fish 13
    ''' FISH 10: Nucleophilic substitution: allylic activation 
        '''
    substructure = 'C=C[CX4H2][Br,I,Cl]'
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

def fish11(molecule: Chem.rdchem.Mol): #Modifiche dell'8/09/22
    ''' FISH 11: Nucleophilic substitution: propargylic activation
        '''
    substructure = 'C#C[CX4H2][Br,I,Cl]'
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

def fish12(molecule: Chem.rdchem.Mol): #Modifiche dell'8/09/22
    ''' FISH 12 Nucleophilic substitution: benzylic activation
        '''
    substructure = '[S]([S]c1ccccc1)c1ccccc1'
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = 'c1ccccc1[CX4H2][Cl,Br,I]'
    substructure2 = Chem.MolFromSmarts(substructure2)
    sub = [substructure, substructure2]
    
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

def fish13(molecule: Chem.rdchem.Mol):
    ''' FISH 13: Nucleophilic substitution: alpha-alo-(C=X,C#X)
        '''   
    substructure = 'Br[CX4H2]C(=O)O*'
    substructure2 = "N#C[CX4H2][Br,I,Cl]"
    #[Br,Cl,I][CX4H2]C(=O)c1ccc([$([NX3](=O)(O));$([NX4H3])])cc1 è giusta ma RDKit non la legge
    substructure3 = "[Br,Cl,I][CX4H2]C(=O)c1ccc([$([NX3](=O)(O))])cc1" 
    substructure4 = "[Br,Cl,I][CX4H2]C(=O)c1ccc([$([NX4H3])])cc1"
    substructure5 = "[Br,Cl,I][CX4H2]C(=O)c1ccc([$([NX4H1])])cc1"
    substructure6 = "[Br,Cl,I][CX4H2]C(=O)c1ccc([$([NX4H2])])cc1"
    substructure7 = "[Br,Cl,I][CX4H2]C(=O)c1ccc([$([S](=O)(=O)(O))])cc1"
    substructure8 = "[Br,Cl,I][CX4H2]C(=O)c1ccc([$(C#N)])cc1"
    substructure9 = "[Br,Cl,I][CX4H2]c1ccc([$([NX3](=O)(O))])cc1"
    substructure10 = "[Br,Cl,I][CX4H2]c1ccc([$([NX4H3])])cc1"
    substructure11 = "[Br,Cl,I][CX4H2]c1ccc([$([NX4H1])])cc1"
    substructure12 = "[Br,Cl,I][CX4H2]c1ccc([$([NX4H2])])cc1"
    substructure13 = "[Br,Cl,I][CX4H2]c1ccc([$([S](=O)(=O)(O))])cc1"
    substructure14 = "[Br,Cl,I][CX4H2]c1ccc([$(C#N)])cc1"
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    substructure4 = Chem.MolFromSmarts(substructure4)
    substructure5 = Chem.MolFromSmarts(substructure5)
    substructure6 = Chem.MolFromSmarts(substructure6)
    substructure7 = Chem.MolFromSmarts(substructure7)
    substructure8 = Chem.MolFromSmarts(substructure8)
    substructure9 = Chem.MolFromSmarts(substructure9)
    substructure10 = Chem.MolFromSmarts(substructure10)
    substructure11 = Chem.MolFromSmarts(substructure11)
    substructure12 = Chem.MolFromSmarts(substructure12)
    substructure13 = Chem.MolFromSmarts(substructure13)
    substructure14 = Chem.MolFromSmarts(substructure14)
    
    sub = [substructure, substructure2,substructure3, substructure4, substructure5, substructure6, substructure7, substructure8, substructure9, substructure10, substructure11, substructure12, substructure13, substructure14]
    
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

def fish14(molecule: Chem.rdchem.Mol):
    ''' FISH 14: Acid anhydrides (acylation)
        '''
    substructure = 'CC(=O)OC(=O)C'
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

def fish15(molecule: Chem.rdchem.Mol):
    ''' FISH 15: Schiff Base
        '''
    #substructure = '[*]-[Ch2]-[Ch]=O'
    substructure = '[c,C][CX3H1]=O'
    substructure2 = "c1(F)c(F)c(F)c(F)c(F)c1C(=O)"
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    sub = [substructure, substructure2]
    
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

def fish16(molecule: Chem.rdchem.Mol):
    ''' FISH 16
        '''
    substructure = '[Ch1]=[Ch1]C(=O)C'#chetone alpha-beta insaturo
    substructure2 = "[Ch1]=[Ch1]C(=O)[NX3H2]" #amide alpha-beta instatura
    substructure3 = "[Ch1]=[Ch1][CX3H1](=O)" #aldeide alpha-beta insatura
    substructure4 = "[Ch1]=[Ch1]C(=O)C#N" #alpha beta insaturo con C#N
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    substructure4 = Chem.MolFromSmarts(substructure4)
    
    sub = [substructure, substructure2, substructure3, substructure4]
    
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

''' 
def fish17(molecule: Chem.rdchem.Mol):
    #FISH 17
     
    substructure = 'C=CC#N'
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
    '''

def fish18(molecule: Chem.rdchem.Mol):
    ''' FISH 18
        '''
    substructure = 'C1[O,N]C1'
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

def fish19(molecule: Chem.rdchem.Mol):
    ''' FISH 19: Proelectrophiles - alcohol dehydrogenase activators
        '''
    substructure = '[Ch]=[Ch][CX4H2][OX2H1]' #con alcool primario
    substructure2 = "[Ch]=[Ch][CX4H1][OX2H1]" #con alcool secondario - la reazione non avviene con alcoli terziari
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    sub = [substructure, substructure2]
    
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

def fish20(molecule: Chem.rdchem.Mol):
    ''' FISH 20
        '''
    substructure = 'BrCCCBr'
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

def fish21(molecule: Chem.rdchem.Mol):
    ''' FISH 21: Proelectrophiles - alcohol dehydrogenase activators
        '''
    substructure = 'C#C[CX4H2][OX2H1]' #come per fish 19
    substructure2 = 'C#C[CX4H1][OX2H1]'
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)

    sub = [substructure, substructure2]
    
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

def fish22(molecule: Chem.rdchem.Mol):
    ''' FISH 22
        '''
    substructure = 'C1(C(C1C)(C)C)C(O)=O'
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

def fish23(molecule: Chem.rdchem.Mol):
    ''' FISH 23
        '''
    substructure = 'O=C(O)C(C)C'
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

def fish24(molecule: Chem.rdchem.Mol):
    ''' FISH 24
        '''
    substructure = 'C1(=CC=C(C=C1)O)[Cl]'
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
'''
def fish25(molecule: Chem.rdchem.Mol):
    #FISH 25 alert di non tossicità
        
    substructure = 'O-C'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def fish26(molecule: Chem.rdchem.Mol):
    #FISH 26- Erika: ho aggiunto * prima della N altrimenti era una molecola a sè
        
        #alert di non tossicità!
    substructure = '*N-C=O-N-C=O'
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
    '''

def fish27(molecule: Chem.rdchem.Mol):
    ''' FISH 27: Cyanogenic Mechanisms - cyanohydrin (rilasciano formaldeide e/o gruppo ciano)
        '''
    substructure = '*C[CX4H1]([OX2H1])C#N'
    substructure2 = "[Cl]c1ccc([NX3H1][CX4H2]C#N)cc1"
    substructure3 = "CC(=O)O[CX4H2]C#N"
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)

    sub = [substructure, substructure2, substructure3]
    
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

def fish28(molecule: Chem.rdchem.Mol):
    ''' FISH 28: Cyanogenic Mechanisms - cyanide
        '''
    substructure = '[Ch1]=[Ch1]C#N'
    substructure2 = "N#C[CX4H2]C#N"
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)

    sub = [substructure, substructure2]
    
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

def fish29(molecule: Chem.rdchem.Mol):
    ''' FISH 29: Cyanogenic Mechanisms - Multistep or Multiple mechanisms - thiocyanates
        '''
    substructure = '[C,c]SC#N'
    substructure2 = "[NX3](=O)(O)aaaa(SC#N)"
    substructure3 = "[NX3H2]aaaa(SC#N)"
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)

    sub = [substructure, substructure2, substructure3]
    
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

def fish30(molecule: Chem.rdchem.Mol):
    ''' FISH 30: Multistep or Multiple mechanisms - alpha-chloro benzyldene malononitrile (addizione di Michael su doppio legame o addizione di acqua e formazione di base di schiff)
        '''
    substructure = '[Cl]c1c([Ch1]=C(C#N)(C#N))cccc1'
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

def fish31(molecule: Chem.rdchem.Mol):
    ''' FISH 31: Multistep or Multiple mechanisms - Tautomerization
        '''
    substructure = 'O=Nc1ccc(N(C)(C))cc1'
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

def ToxRead0(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'O=C(OCc1cccc(Oc2ccccc2)c1)C3CC3'
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

def ToxRead1(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'C1=CC2CC1C3CCCC23'
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

def ToxRead2(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'O=C(OCc1ccccc1)C2CC2'
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

def ToxRead3(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'O=C(OC)C1C(C=C)C1(C)C'
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

def ToxRead4(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = '[#7]-1-[#16]-c2ccccc2-[#6]-1=O'
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

def ToxRead5(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'N#Cc1cn(c(c1))'
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

def ToxRead6(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Cc2ccccc2)cc1'
	substructure2  = 'Cl-C(-Cl)(-Cl)-C(-c1:c:c:c:c:c:1)-c2:c:c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead8(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Cc2ccccc2)cc1'
	substructure2  = 'O-c1:c(:c:c(:c:c:1)-Cl)-C-c2:c:c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead10(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N#C-C(-O-C(-C)=O)-c1:c:c(:c:c:c:1)-O-c2:c:c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead12(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N#C-C(-O-C(-C1-C(-C-1-C)(-C)-C)=O)-c2:c:c(:c:c:c:2)-O-c3:c:c:c:c:c:3'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead14(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N#C-C(-O-C(-C1-C(-C-1-C=C)(-C)-C)=O)-c2:c:c(:c:c:c:2)-O-c3:c:c:c:c:c:3'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead16(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'O=C(Nc1ccccc1)c2ccccc2'
	substructure2  = 'Cl-c1:c:c(:c:c:c:1)-N'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead18(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'C1=CCCCC1'
	substructure2  = 'Cl-C12-C3-C(-C(-C-2(-Cl)-Cl)(-C(=C-1-Cl)-Cl)-Cl)-C-C-C-3'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead20(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'C1=CC2CCC1C2'
	substructure2  = 'Cl-C12-C3-C(-C(-C-2(-Cl)-Cl)(-C(=C-1-Cl)-Cl)-Cl)-C-C-C-3'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead22(molecule: Chem.rdchem.Mol):
	# Alerts NoTOx  ( > 100 mg/l)
	substructure = 'O=C1NC=CC=C1'
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

def ToxRead23(molecule: Chem.rdchem.Mol):
	# coppie NoTOx alerts 
	substructure = 'c1ccncc1'
	substructure2  = 'O-c1:c:n:c:c:c:1'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead25(molecule: Chem.rdchem.Mol):
	# coppie NoTOx alerts 
	substructure = 'c1ccncc1'
	substructure2  = 'O=C1-N-C=C-C=C-1'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead27(molecule: Chem.rdchem.Mol):
	# coppie NoTOx alerts 
	substructure = 'c1ccncc1'
	substructure2  = 'n1:c:c:c(:c:c:1)-C-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead29(molecule: Chem.rdchem.Mol):
	# coppie NoTOx alerts 
	substructure = 'c1cncnc1'
	substructure2  = 'O=C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead31(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc2ccccc2c1'
	substructure2  = 'O'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead33(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc2ccccc2c1'
	substructure2  = 'O-c1:c2:c(:c:c:c:1):c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead35(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc2ccccc2c1'
	substructure2  = 'O-c1:c:c:c:c:c:1'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead37(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead39(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N-c1:c:c:c:c:c:1'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead41(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'O(-c1:c:c:c(:c:c:1)-O)-c2:c:c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead43(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'O=C1NC(=O)c2ccccc21'
	substructure2  = 'N1(-C(-c2:c(-C-1=O):c:c:c:c:2)=O)-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead45(molecule: Chem.rdchem.Mol):
	# coppie Alerts 10-100 mg/l Toxic
	substructure = 'C1CCCCC1'
	substructure2  = 'N(-C1-C-C-C-C-C-1)-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead47(molecule: Chem.rdchem.Mol):
	# coppie Alerts 10-100 mg/l Toxic
	substructure = 'C1CCCCC1'
	substructure2  = 'N-C1-C-C-C(-C-C-1)-C-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead49(molecule: Chem.rdchem.Mol):
	# coppie Alerts 10-100 mg/l Toxic
	substructure = 'C1CCCCC1'
	substructure2  = 'O=C1-C(-C-C-C-C-1)-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

def ToxRead51(molecule: Chem.rdchem.Mol):
	# coppie Alerts 10-100 mg/l Toxic
	substructure = 'O=C1CCCCC1'
	substructure2  = 'O=C1-C(-C-C-C-C-1)-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	sub = [substructure, substructure2]
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
		if j==0:
			out_atm = list_hit_atss
			out_bonds = highlightBonds
		else:
			out_atm = out_atm + list_hit_atss
			out_bonds = out_bonds + highlightBonds
	return (out_bonds, out_atm)

# indices:

def all_index_bonds(molecule: Chem.rdchem.Mol):
    idx = [ fish1(molecule)[0]+
            fish2(molecule)[0]+
            fish3(molecule)[0]+
            fish4(molecule)[0]+
            fish5(molecule)[0]+
            fish6(molecule)[0]+
            fish7(molecule)[0]+
            fish8(molecule)[0]+
            fish9(molecule)[0]+
            fish10(molecule)[0]+
            fish11(molecule)[0]+
            fish12(molecule)[0]+
            fish13(molecule)[0]+
            fish14(molecule)[0]+
            fish15(molecule)[0]+
            fish16(molecule)[0]+
            #fish17(molecule)[0]+
            fish18(molecule)[0]+
            fish19(molecule)[0]+
            fish20(molecule)[0]+
            fish21(molecule)[0]+
            fish22(molecule)[0]+
            fish23(molecule)[0]+
            fish24(molecule)[0]+
            fish27(molecule)[0]+
            fish28(molecule)[0]+
            fish29(molecule)[0]+
            fish30(molecule)[0]+
            fish31(molecule)[0]+
            #fish26(molecule)[0]+
            ToxRead0(molecule)[0]+
            ToxRead1(molecule)[0]+
            ToxRead2(molecule)[0]+
            ToxRead3(molecule)[0]+
            ToxRead4(molecule)[0]+
            ToxRead5(molecule)[0]+
            ToxRead6(molecule)[0]+            
            ToxRead8(molecule)[0]+           
            ToxRead10(molecule)[0]+
            ToxRead12(molecule)[0]+
            ToxRead14(molecule)[0]+
            ToxRead16(molecule)[0]+
            ToxRead18(molecule)[0]+
            ToxRead20(molecule)[0]+
            ToxRead22(molecule)[0]+
            ToxRead23(molecule)[0]+
            ToxRead25(molecule)[0]+
            ToxRead27(molecule)[0]+
            ToxRead29(molecule)[0]+
            ToxRead31(molecule)[0]+
            ToxRead33(molecule)[0]+
            ToxRead35(molecule)[0]+
            ToxRead37(molecule)[0]+
            ToxRead39(molecule)[0]+
            ToxRead41(molecule)[0]+
            ToxRead43(molecule)[0]+
            ToxRead45(molecule)[0]+
            ToxRead47(molecule)[0]+
            ToxRead49(molecule)[0]+
            ToxRead51(molecule)[0]]
                
    return idx[0]
    
def all_index_atoms(molecule: Chem.rdchem.Mol):

    idx = [ np.array(fish1(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish2(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish3(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish4(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish5(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish6(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish7(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish8(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish9(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish10(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish11(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish12(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish13(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish14(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish15(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish16(molecule)[1]).reshape(1,-1).tolist()[0]+
            #np.array(fish17(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish18(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish19(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish20(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish21(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish22(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish23(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish24(molecule)[1]).reshape(1,-1).tolist()[0]+
            #np.array(fish26(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish27(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish28(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish29(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish30(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(fish31(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead0(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead1(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead2(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead3(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead4(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead5(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead6(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead8(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead10(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead12(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead14(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead16(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead18(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead20(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead22(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead23(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead25(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead27(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead29(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead31(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead33(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead35(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead37(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead39(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead41(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead43(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead45(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead47(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead49(molecule)[1]).reshape(1,-1).tolist()[0]+
            np.array(ToxRead51(molecule)[1]).reshape(1,-1).tolist()[0]]

    return idx[0]



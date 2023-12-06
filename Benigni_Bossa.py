'''Missing some count about alogens and others'''

#import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import rdqueries

def f_SA1(molecule: Chem.rdchem.Mol):
    
    ''' SA1_gen. Acyl halide RC(=O)[Br,Cl,F,I], where R is not OH or SH
        '''
    substructure = '[!$([OH1,SH1])]C(=O)[Br,Cl,F,I]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA2(molecule: Chem.rdchem.Mol):
    
    ''' SA2_gen. Methyl, ethyl, propyl, butyl or benzyl esters of sulphonic or phosphonic acid. P(=O)(O)(O)R or S(=O)(O)(O)R where R is not S or O The alkyl chains can have halogen substituents
        '''
    substructure = 'S([!$([OH1,SH1])])(=O)(=O)O([$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$([CH2]c1ccccc1)])' #error
    substructure2 =    'P(=O)([!$([OH1,SH1])])(O([$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$([CH2]c1ccccc1)]))O([$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$(C([#1,Cl,Br,I,F])(C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F]))C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])C([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])([#1,Cl,Br,I,F])),$([CH2]c1ccccc1)])'#errore
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices = molecule.GetSubstructMatches(substructure)
    indices2 = molecule.GetSubstructMatches(substructure2)      
    return len(indices)+len(indices2)

def f_SA3(molecule: Chem.rdchem.Mol):
    ''' SA3_gen. N-methylol derivativesH
        '''
    substructure = "[CX4H2](N)([OX2H1])"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA4(molecule: Chem.rdchem.Mol):
    ''' SA4_gen. This alert contains halogenated olefins where at least one hydrogen or alkyl group is attached to each carbon atom'''
     
    substructure = '[C;D1]=[C;D2]([Cl,F,Br,I])'
    substructure1 = '[C;D1]=[C;D3]([Cl,F,Br,I])([C;!$(C=*);!$(C#*)])'
    substructure2 = "[C;D2]([!Cl;!Br;!F;!I;!$(C=O)])=[C;D2]([Cl,F,Br,I])"
    substructure3 = '[C;D2]([!Cl;!Br;!F;!I;!$(C=O)])=[C;D3]([Cl,F,Br,I])([C;!$(C=*);!$(C#*)])'
    substructure4 = '[C;D3]([!Cl;!Br;!F;!I;!$(C=O)])([C;!$(C=*);!$(C#*)])=[C;D2]([Cl,F,Br,I])'
    substructure5 = '[C;D3]([!Cl;!Br;!F;!I;!$(C=O)])([C;!$(C=*);!$(C#*)])=[C;D3]([Cl,F,Br,I])([C;!$(C=*);!$(C#*)])'
    
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    substructure4 = Chem.MolFromSmarts(substructure4)
    substructure5 = Chem.MolFromSmarts(substructure5)
    
    indices = molecule.GetSubstructMatches(substructure)
    indices1 = molecule.GetSubstructMatches(substructure1)  
    indices2 = molecule.GetSubstructMatches(substructure2)  
    indices3 = molecule.GetSubstructMatches(substructure3)  
    indices4 = molecule.GetSubstructMatches(substructure4)  
    indices5 = molecule.GetSubstructMatches(substructure5)   
    tot = len(indices)+len(indices1)+len(indices2)+len(indices3)+len(indices4)+len(indices5)
    return tot

def f_SA5(molecule: Chem.rdchem.Mol):
    ''' SA5_gen. S or N mustard
        '''
    substructure = "[F,Cl,Br,I][CX4H2][CX4H2][N,S][CX4H2][CX4H2][F,Cl,Br,I]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA6(molecule: Chem.rdchem.Mol):
    ''' SA6_gen. propiolactones and propiosultones
        '''
    substructure = "[O,S]=C1[O,S]CC1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    substructure2 = "O=S1(=O)(CCCO1)"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)     
    return len(indices) + len(indices2)

def f_SA7(molecule: Chem.rdchem.Mol):
    ''' SA7_gen. epoxides and aziridines
        '''
    substructure = "C1[O,N]C1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA8(molecule: Chem.rdchem.Mol):
    ''' SA8_gen. This alert contains non tertiary aliphatic halogens
        '''
    substructure = "[$([CX4!H0;R0]);!$(C([#1,C])=[O,C]);!$([CX4H2]([F,Cl,Br,I])[CX4H2][N,S][CX4H2][CX4H2][F,Cl,Br,I]);!$(COP(O)(=O));!$(COS(=O)(=O));!$(CCOP(O)(=O));!$(CCOS(=O)(=O));!$(CCCOP(O)(=O));!$(CCCOS(=O)(=O));!$(CCCCOP(O)(=O));!$(CCCCOS(=O)(=O));!$(CCCCCOP(O)(=O));!$(CCCCCOS(=O)(=O));!$(CCCCCCOP(O)(=O));!$(CCCCCCOS(=O)(=O))][Cl,Br,I]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA9(molecule: Chem.rdchem.Mol):
    ''' SA9_gen.
        '''
    substructure = "O=[NX2]OC"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
def f_SA11(molecule: Chem.rdchem.Mol):
    ''' SA11_gen. Aliphatic and aromatic aldehydes
        '''
    substructure = "[#6][$([CX3H1]);!$(CC=C)](=O)"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA12(molecule: Chem.rdchem.Mol):
    ''' SA12_gen. Quinones
        '''
    substructure = "O=[#6]1[#6]=,:[#6][#6](=O)[#6]=,:[#6]1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    
    substructure2 = "O=[#6]1[#6]=,:[#6][#6]=,:[#6][#6]1(=O)"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)       
    return len(indices) + len(indices2)

def f_SA13(molecule: Chem.rdchem.Mol):
    ''' SA13_gen. This applies to molecules that contain a NN group not in a ring, and not NN=O.
        '''
    substructure = "[N+0]!@;-[N+0](=[!O;!N])"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    
    # substructure2 = "[N+0]([#1,*])!@;-[N+0]([#1,*])"
    substructure2 = "[$([N+0;D1]),$([N+0;D2](-*)(-N)),$([N+0;D3](-*)(-*)(-N))]!@;-[$([N+0;D1]),$([N+0;D2](-*)(-N)),$([N+0;D3](-*)(-*)(-N))]"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)  
    '''
    5/09 modifiche alert 13
    #substructure3 = "[N+0]([#1,*])!@;-[N+0]" #modifche 12/01 maybe better this for hydrazine '[NX3H2]-[NX3H]'
    substructure3 = "[NX3H]-[NX3H]"
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices3 = molecule.GetSubstructMatches(substructure3) 
    
    substructure4 = "[NX3H]-[NX3H2]"
    substructure4 = Chem.MolFromSmarts(substructure4)
    indices4 = molecule.GetSubstructMatches(substructure4) 
    '''
    return len(indices) + len(indices2) # + len(indices3) + len(indices4) 
    
def f_SA14(molecule: Chem.rdchem.Mol):
    ''' SA14_gen. Aliphatic azo and azoxy.
        '''
    substructure = "[C,#1]N=[NX2][C,#1]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    
    substructure2 = "[$(C=[N+]=[N-]);!$(C=[N+]=[N-]=N);!$(C=[N+]=[N-]N)]"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)
    
    substructure3 = "C=[$(N=N);!$(N=N=N);!$(N=NN)]"
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices3 = molecule.GetSubstructMatches(substructure3)
    
    substructure4 = "CN=NO"
    substructure4 = Chem.MolFromSmarts(substructure4)
    indices4 = molecule.GetSubstructMatches(substructure4)
          
    return len(indices) + len(indices2) + len(indices3) + len(indices4)
    
def f_SA15(molecule: Chem.rdchem.Mol):
    ''' SA15_gen. Isocyanate and isothiocyanate groups
        '''
    substructure = "[NX2]=C=[O,S]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA16(molecule: Chem.rdchem.Mol):
    ''' SA16_gen. Alkyl carbamate and thiocarbamate
        ''' #alberto
    substructure = "[NH2;D1]C(=[O,S])[O,S][CX4]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    substructure2 = "[NH;D2]([C;!$(C=*);!$(C#*)])C(=[O,S])[O,S][CX4]"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2) 
    
    substructure3 = "[N;D3]([C;!$(C=*);!$(C#*)])C(=[O,S])[O,S][CX4]" #Modificato da Erika
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices3 = molecule.GetSubstructMatches(substructure3) 
    return len(indices) + len(indices2) + len(indices3)

def f_SA18(mol, includeSpiro=False):
    
    '''SA18 Polycyclic Aromatic Hydrocarbons:
    Polycyclic Aromatic Hydrocarbons, with three or more fused rings. Does not include heterocyclic compounds'''
    q = rdqueries.IsAromaticQueryAtom()
    test_arom = [x.GetIdx() for x in mol.GetAtomsMatchingQuery(q)]
    out = 0
    if len(test_arom)>10:
        ri = mol.GetRingInfo()
        systems = []
        for ring in ri.AtomRings():
            ringAts = set(ring)
            nSystems = []
            for system in systems:
                nInCommon = len(ringAts.intersection(system))
                if nInCommon and (includeSpiro or nInCommon>1):
                    ringAts = ringAts.union(system)
                else:
                    nSystems.append(system)
            nSystems.append(ringAts)
            systems = nSystems
            if len(systems)>0:
                for ringsist in systems:
                    result =  all(elem in test_arom  for elem in ringsist)
                    if len(ringsist)>10 and result:
                        carbon_atoms = []
                        for atom in mol.GetAtoms():
                            if atom.GetIdx() in list(ringsist):
                                carbon_atoms.append(atom.GetAtomicNum())
                        if set(carbon_atoms) == {6}:
                            # case = [3,4] 
                            # for n_case in case:
                                # if len(ringsist)-6*n_case + ((n_case-1)*2) == 0:
                            out = out + 1 
    return out


def f_SA19(mol, includeSpiro=False):
    '''  
    SA19 Heterocyclic Polycyclic Aromatic Hydrocarbons: 
        Heterocyclic Polycyclic Aromatic Hydrocarbons (3 or more fused rings).
    '''
    q = rdqueries.IsAromaticQueryAtom()
    test_arom = [x.GetIdx() for x in mol.GetAtomsMatchingQuery(q)]
    out = 0
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon>1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
        if len(systems)>0:
            for ringsist in systems:
                result =  all(elem in test_arom  for elem in ringsist)
                if len(ringsist)>10 and result:
                    carbon_atoms = []
                    for atom in mol.GetAtoms():
                        if atom.GetIdx() in list(ringsist):
                            carbon_atoms.append(atom.GetAtomicNum())
                    if len(set(carbon_atoms))>1:
                        # case = [3,4,5] 
                        # for n_case in case:
                            # if len(ringsist)-6*n_case + ((n_case-1)*2) == 0:
                        out = out + 1 
    return out


'''
def f_SA18(molecule: Chem.rdchem.Mol):
    # SA18_gen. Polycyclic Aromatic Hydrocarbons, with three or more fused rings
        
    substructure = "c:c(:c):c"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    if len(indices) >3:
        a = 1
    else:
        a = 0
    return a

def f_SA19(molecule: Chem.rdchem.Mol):
    # SA19_gen. Heterocyclic Polycyclic Aromatic Hydrocarbons (3 or more fused rings).
        
    substructure = "c:c(:c):c" #ma se sono condensati tipo antracene?
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    if len(indices) >3:
        a = 1
    else:
        a = 0
    return a
'''
def f_SA20(mol,includeSpiro=False):
    '''SA20 (Poly) Halogenated Cycloalkanes (Nongenotoxic carcinogens)'''
    # smart alogeni
    substructure = Chem.MolFromSmarts('[Br,Cl,F,I]')
    # smart cicli non aromatici
    substructure1 = Chem.MolFromSmarts('[C;R]')
    
    indices = mol.GetSubstructMatches(substructure)
    ind_ring = mol.GetSubstructMatches(substructure1)
    
    out = 0
    # sono presenti degli alogeni e degli anelli
    if len(indices)>2 and len(ind_ring):
        ri = mol.GetRingInfo()
        systems = []
        for ring in ri.AtomRings():
            ringAts = set(ring)
            nSystems = []
            for system in systems:
                nInCommon = len(ringAts.intersection(system))
                if nInCommon and (includeSpiro or nInCommon>1):
                    ringAts = ringAts.union(system)
                else:
                    nSystems.append(system)
            nSystems.append(ringAts)
            systems = nSystems
        # ciclo sui ring trovati
        for ring in systems:
            n = 0
            # se i-esimo ring è alifatico.. continua
            test = [x[0] for x in ind_ring]
            result =  all(elem in test  for elem in ring)
            if result:
                # verifica che su uno stesso sistema ci siano almeno tre alogeni
                for bond in mol.GetBonds():
                    Begin = bond.GetBeginAtomIdx()
                    End = bond.GetEndAtomIdx()
                    # se almeno uno dei due ind del bond sono sul ciclo
                    if (End in ring or Begin in ring) and not (End in ring and Begin in ring):
                        test_alogen = [x[0] for x in indices]
                        # se l'altro è tra ind degli alogeni
                        if End in test_alogen or Begin in test_alogen:
                            n = n+1
                            if n>2:
                                out = out+1
                        
    return out

def f_SA21(molecule: Chem.rdchem.Mol):
    ''' SA21_gen. Alkyl and aryl N-nitroso groups
        '''
    substructure = "[C,c]N[NX2;v3]=O"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA22(molecule: Chem.rdchem.Mol):
    ''' SA22_gen. Azide and triazene groups
        '''
    substructure = "[N]=[N]-[N]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    
    substructure2 = "[N]=[N]=[N]"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)     
    return len(indices) + len(indices2)
    
def f_SA23(molecule: Chem.rdchem.Mol):
    ''' SA23_gen. Aliphatic N-nitro
        '''
    substructure = "[C!r][NH1]N(=O)O"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    
    substructure2 = "[C!r]N(A)N(=O)O"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)     
    return len(indices) + len(indices2)
    
def f_SA24(molecule: Chem.rdchem.Mol):
    ''' SA24_gen. unsaturated alkoxy
        '''
    substructure = "[!$([#6](=O)[!O]),#1][C!H0;!R]([!$([#6](=O)[!O]),#1])!@;=[C!H0;!R]O[#6]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA25(molecule: Chem.rdchem.Mol):
    ''' SA25_gen. Aromatic nitroso group
        '''
    substructure = "[a;!$(a(a[A;!#1])(a[A;!#1]));!$(aa[CX3](=O)[OX2H1]);!$(aa[SX4](=[OX1])(=[OX1])([O]));!$(aaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaaaa[SX4](=[OX1])(=[OX1])([O]))]!@[$([NX2]=O)]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA26(molecule: Chem.rdchem.Mol):
    ''' SA26_gen. Aromatic ring N-oxide
        '''
    substructure = "[n+]!@[O-]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA27(molecule: Chem.rdchem.Mol):
    #alberto
    ''' SA27_gen. Nitro aromatic,SA27 Nitro aromatic: Nitro aromatic. However: Aromatic nitro groups with ortho-disubstitution or with a carboxylic acid substituent in ortho position should be excluded. Please note that a molecule like this <b>CC1=CC=CC(=C1[N+](=O)[O-])[N+](=O)[O-]</b> should be included in the alert: one of the two nitro groups is ortho disubstituted, but the other one is ortho-monosubstituted. Also the following molecule <b>CC2=CC=CC(CCC1=CC=CC(=C1)[N+](=O)[O-])=C2[N+](=O)[O-]</b> Should fire the alert (one nitro group is ortho disubstituted, but the other is not). If a sulfonic acid group (-SO3H) is present on the ring that contains also the nitro group, the substance should be excluded. 
        '''
    substructure = "[a;!$(a(a[A;!#1;!H])(a[A;!#1;!H]));!$(aa[CX3](=O)[OX2H1]);!$(aa[SX4](=[OX1])(=[OX1])([OX2H1]));!$(aaa[SX4](=[OX1])(=[OX1])([OX2H1]));!$(aaaa[SX4](=[OX1])(=[OX1])([OX2H1]));!$(aaaaa[SX4](=[OX1])(=[OX1])([OX2H1]));!$(aaaaaa[SX4](=[OX1])(=[OX1])([OX2H1]))]([$([N+]([O-])=O)])" 
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA28(molecule: Chem.rdchem.Mol):
    ''' SA28_gen. Primary aromatic amine, hydroxyl amine and its derived esters (with restrictions). SA28 Primary aromatic amine, hydroxyl amine and its derived esters (with restrictions): Primary aromatic amine, hydroxyl amine and its derived esters (with restrictions). However: Aromatic amino groups with ortho disubstitutions or with a carboxylic acid substituent in ortho position are excluded. If a sulfonic acid group (-SO3H) is present on the ring that contains also the amino group, the substance should be excluded from the alert. The following structures should also be included: O=C=NC1=CC=CC=C1 and C([H])([H])=NC1=CC=CC=C1. The possibility that the Nitrogen atom of hydroxyl amine is part of a cycle, should be excluded.
        '''
    substructure = "[a;!$(a(a[A;!#1])(a[A;!#1]));!$(aa[CX3](=O)[OX2H1]);!$(aa[SX4](=[OX1])(=[OX1])([O]));!$(aaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaaaa[SX4](=[OX1])(=[OX1])([O]))]!@[$([N+0;H2;D1]),$([N+0;H1;D2][OH;D1]),$([N+0;H0;D3]([OH;D1])C),$([N+0;H1;D2]OC=O),$([N+0;H0;D3](C)OC=O)]" #error (sia alberto che non)
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    substructure2 = "aN=C=O"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)  
    
    substructure3 = "aN=[CH2]"
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices3 = molecule.GetSubstructMatches(substructure3)  
           
    return len(indices)+ len(indices2) + len(indices3)

def f_SA28bis(molecule: Chem.rdchem.Mol):
    ''' SA28bis_gen. Mono- or di- methyl or ethyl aromatic amines, are included. SA28bis Aromatic mono- and dialkylamine: Mono- or di- methyl or ethyl     aromatic amines, are included. However:Aromatic amino groups with ortho-disubstitution or with a carboxylic acid substituent in ortho position should be excluded. If a sulfonic acid group (-SO3H) is present on the ring that contains also the amino group, the substance should be excluded from the alert.
        '''
    substructure = "[a;!$(a(a[A;!#1])(a[A;!#1]));!$(aa[CX3](=O)[OX2H1]);!$(aa[SX4](=[OX1])(=[OX1])([O]));!$(aaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaaaa[SX4](=[OX1])(=[OX1])([O]))]!@[$([NX3;v3]([#1,CH3])([CH3])),$([NX3;v3]([#1,CH3])([CH2][CH3])),$([NX3;v3]([CH2][CH3])([CH2][CH3]))]"#error (sia alberto che non)
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA28ter(molecule: Chem.rdchem.Mol):
    '''SA28ter_gen. Aromatic N-acyl amine
        '''
    substructure = "[a;!$(a(a[A;!#1])(a[A;!#1]));!$(aa[CX3](=O)[OX2H1]);!$(aa[SX4](=[OX1])(=[OX1])([O]));!$(aaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaaa[SX4](=[OX1])(=[OX1])([O]));!$(aaaaaa[SX4](=[OX1])(=[OX1])([O]))]!@[$([NX3;v3]([#1,CH3])C(=O)([#1,CH3]))]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA29(molecule: Chem.rdchem.Mol):
    ''' SA29_gen. Aromatic diazo
        '''
    substructure = "a[N]=[N]a"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    substructure2 = "[$([N](a)=[N]a);!$([N](a:aS(=[OD1])(=[OD1])[OD1])=[N]a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1]);!$([N](a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])=[N]a:a:a:a:a:a:aS(=[OD1])(=[OD1])[OD1])]" #alberto
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)     
    return len(indices) + len(indices2)
 
def f_SA30(molecule: Chem.rdchem.Mol):
    ''' SA30_gen. Coumarins and Furocoumarins
        '''
    substructure = "O=c1ccc2ccccc2(o1)"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    substructure = "O=C1C=Cc2ccccc2O1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)        
    return len(indices)
    
def f_SA37(molecule: Chem.rdchem.Mol):
    ''' SA37_gen. Genotoxic mechanism. Pyrrolizidine Alkaloids
        '''
    substructure = "C12CCCN1CC=C2"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA38(molecule: Chem.rdchem.Mol):
    ''' SA38_gen. Genotoxic mechanism. Alkenylbenzenes: Alkenylbenzenes
        '''
    substructure = "c1ccccc1C[C;!R]=C"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA39(molecule: Chem.rdchem.Mol):
    ''' SA39_gen_and_nogen. Mixed genotoxic and non genotoxic mechanism. Steroidal estrogens: Steroidal estrogens
        '''
    substructure = "C[C@@]12CCC3c4c(CCC3C1CC[C@H]2O)cc(O)cc4"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return len(indices)
       
def f_SA10(molecule: Chem.rdchem.Mol): #28/03 modifica alert 10
    ''' SA10_gen. ?,? unsaturated carbonyls.
     '''    

    substructure = "[!a,#1;!$(C1(=O)C=CC(=O)C=C1)][#6]([!a,#1;!$(C1(=O)C=CC(=O)C=C1)])!:;=[#6][#6](=O)[!O;!$([#6]1:,=[#6][#6](=O)[#6]:,=[#6][#6](=O)1)]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    substructure1 = "[$([#6]);!$([#6]1C(=O)[#6]:,=[#6]C(=O)[#6]:,=1);!$([#6]C!=;-[C;R]!=;-[C;R]!=;-[C;R]!=;-[C;R]!=;-[C;R])]!:;=[$([#6]);!$(C=C[a])][$([#6]);!$([#6]-O);!$(C1(=O)[#6]:,=[#6]C(=O)[#6]:,=[#6]1)](=O)"
    substructure1 = Chem.MolFromSmarts(substructure1)
    indices1 = molecule.GetSubstructMatches(substructure1)
    
    
    return len(indices1) + len(indices)
    
def f_SA17(molecule: Chem.rdchem.Mol):
    ''' SA17_nogen. Thiocarbonyl (Nongenotoxic carcinogens)
        '''
    substructure = "[#7X3][#6](=[SX1])[!$([O,S][CX4])!$([OH,SH])!$([O-,S-])]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA31a(molecule: Chem.rdchem.Mol):
    ''' SA31a_nogen. Halogenated benzene (Nongenotoxic carcinogens)
        '''
    substructure = "[c;!$(c[N+](=O)[O]);!$(c[N]([#1,C])([#1,C]));!$(cN([OX2H])([#1,C]));!$(cN([#1,C])OC=O);!$(c[NX3v3]([#1,CH3])([#1,CH3]));!$(c[NX3v3]([#1,CH3])([CH2][CH3]));!$(c[NX3v3]([CH2][CH3])([CH2][CH3]));!$(cNC(=O)[#1,CH3]);!$(cN=[N]a);!$(c!@[cR1r6]1ccccc1);!$(c!@*!@c1ccccc1);!$([R2])]1[c;!$(c[N+](=O)[O]);!$(c[N]([#1,C])([#1,C]));!$(cN([OX2H])([#1,C]));!$(cN([#1,C])OC=O);!$(c[NX3v3]([#1,CH3])([#1,CH3]));!$(c[NX3v3]([#1,CH3])([CH2][CH3]));!$(c[NX3v3]([CH2][CH3])([CH2][CH3]));!$(cNC(=O)[#1,CH3]);!$(cN=[N]a);!$(c!@[cR1r6]1ccccc1);!$(c!@*!@c1ccccc1);!$([R2])][c;!$(c[N+](=O)[O]);!$(c[N]([#1,C])([#1,C]));!$(cN([OX2H])([#1,C]));!$(cN([#1,C])OC=O);!$(c[NX3v3]([#1,CH3])([#1,CH3]));!$(c[NX3v3]([#1,CH3])([CH2][CH3]));!$(c[NX3v3]([CH2][CH3])([CH2][CH3]));!$(cNC(=O)[#1,CH3]);!$(cN=[N]a);!$(c!@[cR1r6]1ccccc1);!$(c!@*!@c1ccccc1);!$([R2])][c;!$(c[N+](=O)[O]);!$(c[N]([#1,C])([#1,C]));!$(cN([OX2H])([#1,C]));!$(cN([#1,C])OC=O);!$(c[NX3v3]([#1,CH3])([#1,CH3]));!$(c[NX3v3]([#1,CH3])([CH2][CH3]));!$(c[NX3v3]([CH2][CH3])([CH2][CH3]));!$(cNC(=O)[#1,CH3]);!$(cN=[N]a);!$(c!@[cR1r6]1ccccc1);!$(c!@*!@c1ccccc1);!$([R2])][c;!$(c[N+](=O)[O]);!$(c[N]([#1,C])([#1,C]));!$(cN([OX2H])([#1,C]));!$(cN([#1,C])OC=O);!$(c[NX3v3]([#1,CH3])([#1,CH3]));!$(c[NX3v3]([#1,CH3])([CH2][CH3]));!$(c[NX3v3]([CH2][CH3])([CH2][CH3]));!$(cNC(=O)[#1,CH3]);!$(cN=[N]a);!$(c!@[cR1r6]1ccccc1);!$(c!@*!@c1ccccc1);!$([R2])]c1([Cl,Br,F,I;!$([Cl,Br,I,F]cc[Cl,Br,I,F]);!$([Cl,Br,I,F]ccc[Cl,Br,I,F]);!$([Cl,Br,I,F]c1c([OX2H])c([OX2H])c([OX2H])cc1);!$([Cl,Br,I,F]c1c([OX2H])c([OX2H])cc([OX2H])c1);!$([Cl,Br,I,F]c1c([OX2H])c([OX2H])ccc1([OX2H]));!$([Cl,Br,I,F]c1c([OX2H])cc([OX2H])c([OX2H])c1);!$([Cl,Br,I,F]c1c([OX2H])cc([OX2H])cc1([OX2H]));!$([Cl,Br,I,F]c1cc([OX2H])c([OX2H])c([OX2H])c1)])"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA31b(molecule: Chem.rdchem.Mol):
    ''' SA31b_nogen. Halogenated PAH (naphthalenes, biphenyls, diphenyls) (Nongenotoxic carcinogens)
        '''
    substructure = "[Cl,Br,F,I]c1ccc2ccccc2(c1)"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    
    substructure2 = "[Cl,Br,F,I]c1ccc(cc1)!@c2ccc(cc2)[Cl,Br,F,I]"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2) 
    
    substructure3 = "c1cc(ccc1[!R]c2ccc(cc2)[Cl,Br,F,I])[Cl,Br,F,I]"
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices3 = molecule.GetSubstructMatches(substructure3)      
    return len(indices) + len(indices2) + len(indices3) 
    
def f_SA31c(molecule: Chem.rdchem.Mol):
    ''' SA31c_nogen. Halogenated dibenzodioxins (Nongenotoxic carcinogens)
        '''
    substructure = "c1ccc2Oc3cc(ccc3(Oc2(c1)))[Cl,Br,F,I]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA40(molecule: Chem.rdchem.Mol):
    ''' SA40_nogen. Nongenotoxic mechanism. Substituted phenoxyacid: Substituted phenoxyacid
        '''
    substructure = "c1(OC(C)(C)C(=O)O)ccc([#6,#17])cc1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    substructure2 = "c1(OCC(=O)[O;H0])cc(Cl)c(Cl)cc1"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)
    
    substructure3 = "c1(OCC(=O)[O;H0])c(Cl)cc(Cl)cc1"
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices3 = molecule.GetSubstructMatches(substructure3)     
    return len(indices) + len(indices2) + len(indices3)
    
def f_SA41(molecule: Chem.rdchem.Mol):
    ''' SA41_nogen. Nongenotoxic mechanism. Substituted n-alkylcarboxylic acids: Substituted n-alkylcarboxylic acids
        '''
    substructure = "[C;!R&$(C([C;!R])([C;!R])[C;!R][C;!R])&!$(CCCCCCCCCCCC)][C;!R&$(C[OX2;!R]),$(C(=O)[OX2;!R])]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA42(molecule: Chem.rdchem.Mol):
    ''' SA42_nogen. Nongenotoxic mechanism. Phthalate diesters and monoesters: Phthalate diesters and monoesters
        '''
    substructure = "O=C(O)c1ccccc1C(=O)O"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    
    substructure2 = "O=C(O)[CX4;!R][CX4;!R][CX4;!R][CX4;!R]C(=O)O"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)             
    return len(indices) + len(indices2)
    
def f_SA43(molecule: Chem.rdchem.Mol):
    ''' SA43_nogen. Nongenotoxic mechanism. SA43 Perfluorooctanoic acid (PFOA): Perfluorooctanoic acid (PFOA)
        '''
    substructure = "CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA44(molecule: Chem.rdchem.Mol):
    ''' SA44_nogen. Nongenotoxic mechanism. SA44 Trichloro (or fluoro) ethylene and Tetrachloro (or fluoro) ethylene: Trichloro (or fluoro) ethylene and Tetrachloro (or fluoro) ethylene
        '''
    substructure = "[Cl,F][C;!$(Cc)]=C([Cl,F])[Cl,F]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)  
    
    substructure2 = "Cl[C;!$(Cc)]=C(Cl)Cl"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2) 
    
    substructure3 = "[Cl,F]C#C[Cl,F]"
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices3 = molecule.GetSubstructMatches(substructure3)
    return len(indices) + len(indices2) +len(indices3)
    
def f_SA45(molecule: Chem.rdchem.Mol):
    ''' SA45_nogen. Ngenotoxic mechanisms. Indole-3-carbinol: Indole-3-carbinol
        '''
    substructure = "OCc1c[nH]c2ccccc12"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA46(molecule: Chem.rdchem.Mol):
    ''' SA46_nogen. Nongenotoxic mechanism. Pentachlorophenol: Pentachlorophenol
        '''
    substructure = "Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1O"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_SA47(molecule: Chem.rdchem.Mol):
    ''' SA47_nogen. Nongenotoxic mechanisms. O-phenylphenol: O-phenylphenol
        '''
    substructure = "Oc2ccccc2c1ccccc1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    
    substructure2 = "Oc1c(c2ccccc2)cccc1"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)     
    return len(indices) + len(indices2)
    
def f_SA48(molecule: Chem.rdchem.Mol):
    ''' SA48_nogen. Nongenotoxic mechanisms. SA48 Quercetin-type flavonoids: Quercetin-type flavonoids
        '''
    substructure = "Oc1cc(O)c2C(=O)C(O)=C(Oc2c1)c3ccc(O)c(O)c3"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA49(molecule: Chem.rdchem.Mol):
    ''' SA49_nogen. Nongenotoxic mechanism. Imidazole and benzimidazole: Imidazole and benzimidazole
        '''
    substructure = "n1c[nH]cc1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    substructure2 = "n2c1ccccc1nc2"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)      
    return len(indices)

def f_SA50(molecule: Chem.rdchem.Mol):
    ''' SA50_nogen. Nongenotoxic mechanism. Dicarboximide: Dicarboximide
        '''
    substructure = "[#6]1[#6](=O)[#7][#6](=O)[#6]1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA51(molecule: Chem.rdchem.Mol):
    ''' SA51_nogen. Nongenotoxic mechanism. Dimethylpyridine: Dimethylpyridine
        '''
    substructure = "[CX4H3]c1cccc([CX4H3])n1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA52(molecule: Chem.rdchem.Mol):
    ''' SA52_nogen. Nongenotoxic mechanism. Metals, oxidative stress: Metals, oxidative stress
        '''
    substructure = "[As,Cu,Cr,Hg,Co]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA53(molecule: Chem.rdchem.Mol):
    ''' SA53_nogen. Nongenotoxic mechanism. SA53 Benzensulfonic ethers: Benzensulfonic ethers
        '''
    substructure = "c1cc(N)ccc1S(=O)(=O)N"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    
    substructure2 = "c1cc(S)ccc1S(=O)(=O)N"
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)
    
    substructure3 = "c1ccccc1S(=O)(=O)[N;-1]"
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices3 = molecule.GetSubstructMatches(substructure3)
    
    substructure4 = "c1cc(C)ccc1S(=O)(=O)O"
    substructure4 = Chem.MolFromSmarts(substructure4)
    indices4 = molecule.GetSubstructMatches(substructure4)
    return len(indices) + len(indices2) + len(indices3) +len(indices4)

def f_SA54(molecule: Chem.rdchem.Mol):
    ''' SA54_nogen. Nongenotoxic mechanism. 1,3-Benzodioxoles: 1,3-Benzodioxoles
        '''
    substructure = "[C;!$(C(C)(C))]1Oc2ccccc2O1"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)

    return len(indices)

def f_SA55(molecule: Chem.rdchem.Mol):
    ''' SA55_nogen. Nongenotoxic mechanisms. Phenoxy herbicides: Phenoxy herbicides
        '''
    substructure = "c1cc(O)ccc1OC(C)C(=O)[O;!$(OCCCC)]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA56(molecule: Chem.rdchem.Mol):
    ''' SA56_nogen. Nongenotoxic mechanism
        '''
    substructure = "[C;H1,!$(CCOP)&!$(CCP)&!$(CF)&!$(CN)&!$(C[CH3])](Cl)(Cl)(Cl)"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_SA57(molecule: Chem.rdchem.Mol):
    ''' SA57 DNA Intercalating Agents with a basic side chain: DNA intercalating agents are defined as those compounds that are able to insert partially or completely between adjacent DNA base pairs.
        '''
    substructure = "[N+0;!R;!$(N=*);!$(N#*)][C+0;!R;!$(C=*);!$(C#*)][C+0;!R;!$(C=*);!$(C#*)][N+0;!R;!$(N=*);!$(N#*)]Cc1c2aaaac2aaa1"
    substructure1 = "[N+0;!R;!$(N=*);!$(N#*)][C+0;!R;!$(C=*);!$(C#*)][C+0;!R;!$(C=*);!$(C#*)][N+0;!R;!$(N=*);!$(N#*)]Cc1c2[A]=[A][A]c2aaa1"
    substructure2 = "[N+0;!R;!$(N=*);!$(N#*)][C+0;!R;!$(C=*);!$(C#*)][C+0;!R;!$(C=*);!$(C#*)][N+0;!R;!$(N=*);!$(N#*)]Cc1c2[A][A]=[A][A]c2aaa1"
    substructure3 = "[N+0;!R;!$(N=*);!$(N#*)][C+0;!R;!$(C=*);!$(C#*)][C+0;!R;!$(C=*);!$(C#*)][N+0;!R;!$(N=*);!$(N#*)]Cc1c2[A][A]=[A]c2aaa1"
    substructure4 = "[N+0;!R;!$(N=*);!$(N#*)][C+0;!R;!$(C=*);!$(C#*)][C+0;!R;!$(C=*);!$(C#*)][N+0;!R;!$(N=*);!$(N#*)]Cc1c2Cc3aaaac3c2aaa1"
    substructure5 = "[N+0;!R;!$(N=*);!$(N#*)][C+0;!R;!$(C=*);!$(C#*)][C+0;!R;!$(C=*);!$(C#*)][N+0;!R;!$(N=*);!$(N#*)]Cc1c2c3aaaac3Cc2aaa1"
    
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    substructure4 = Chem.MolFromSmarts(substructure4)
    substructure5 = Chem.MolFromSmarts(substructure5)
    
    indices = molecule.GetSubstructMatches(substructure)
    indices1 = molecule.GetSubstructMatches(substructure1)
    indices2 = molecule.GetSubstructMatches(substructure2)
    indices3 = molecule.GetSubstructMatches(substructure3)
    indices4 = molecule.GetSubstructMatches(substructure4)
    indices5 = molecule.GetSubstructMatches(substructure5)
    tot = len(indices)+len(indices1)+len(indices2)+len(indices3)+len(indices4)+len(indices5)
    return tot
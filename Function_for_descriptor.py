
from rdkit import Chem

def f_steroid(molecule: Chem.rdchem.Mol):
    ''' steroid scaffold
        '''
    substructure = '[#6]12~[#6]~[#6]~[#6]3~[#6](~[#6]1~[#6]~[#6]~[#6]2)~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~[#6]34'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def f_five_ring_hetero(molecule: Chem.rdchem.Mol):
    ''' five_ring_hetero
        '''
    substructure = '[#6]1~[#6]~[#6]~[#6]~[!C;!c]1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)   
    return len(indices)
    
def f_six_ring_hetero(molecule: Chem.rdchem.Mol):
    ''' six_ring_hetero '''
    substructure = '[#6]1~[#6]~[#6]~[#6]~[#6]~[!C;!c]1'
    substructure = Chem.MolFromSmarts(substructure)
    # find substructures
    indices = molecule.GetSubstructMatches(substructure)    
    return len(indices)
    
def f_benzCH2(molecule: Chem.rdchem.Mol):
    ''' benzCH2 present in molecule'''
    substructure = 'c1c([C&H2])cccc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)    
    return len(indices)
    
def f_heteroatoms_heterocycles6(molecule: Chem.rdchem.Mol):
    ''' heteroatoms (O or S) bonded to a 6 membered heterocyclic compound'''  
    all_indices = []  
    substructures = ['[O,S]-[#6]1~[#6]~[#6]~[#6]~[#6]~[!C;!c]1',
                    '[#6]1~[#6]([O,S])~[#6]~[#6]~[#6]~[!C;!c]1',
                    '[#6]1~[#6]~[#6]([O,S])~[#6]~[#6]~[!C;!c]1']
    
    substructures = map(Chem.MolFromSmarts, substructures)  
    for fragment in substructures:
        indices = molecule.GetSubstructMatches(fragment)
        all_indices += list(indices)     
    return len(tuple(all_indices))

def f_heteroatoms_heterocycles5(molecule: Chem.rdchem.Mol):
    ''' heteroatoms (O or S) bonded to a 5 membered
        heterocyclic compound'''  
    all_indices = []  
    substructures = ['[#6]1~[#6]([O,S])~[#6]~[#6]~[!C;!c]1',
                    '[#6]1~[#6]~[#6]~[#6]([O,S])~[!C;!c]1']
    
    substructures = map(Chem.MolFromSmarts, substructures)  
    for fragment in substructures:
        indices = molecule.GetSubstructMatches(fragment)
        all_indices += list(indices)     
    return len(tuple(all_indices))
    
def f_biphenol(molecule: Chem.rdchem.Mol):
    ''' identify catechol, resorcinol,
        and hydroquinone
        ''' 
    all_indices = []   
    substructures = ['[OH]c1c([OH])cccc1',
                    '[OH]c1cc([OH])ccc1',
                    '[OH]c1ccc([OH])cc1']  
    substructures = map(Chem.MolFromSmarts, substructures)
    for fragment in substructures:
        indices = molecule.GetSubstructMatches(fragment)
        all_indices += list(indices)      
    return len(all_indices)
    
def f_benzaldehyde(molecule: Chem.rdchem.Mol):
    ''' benzaldehyde present in molecule
        ''' 
    substructure = 'c1ccccc1-[C;H1](=O)'
    substructure = Chem.MolFromSmarts(substructure)
    # find substructures
    indices = molecule.GetSubstructMatches(substructure)     
    return len(indices)

#16/12 add: modify the calculator function
#16/12 add: modify the calculator function
#16/12 add: modify the calculator function

def f_Ar_OR(molecule: Chem.rdchem.Mol):
    ''' aromatic ring with OR sostituent'''
    
    substructure = 'a-[OX2]-[C^3]'
    #substructure1 = 'a-[OX2]-[CD3]-[CX4H3]'
    #substructure2 = 'a-[OX2]-[CD4C]-[CX4H3]' #carbonio quaternario legato all'O-Ar
    #substructure3 = 'a-[OX2]-[CX4H3]'
    
    substructure = Chem.MolFromSmarts(substructure)
    #substructure1 = Chem.MolFromSmarts(substructure1)
    #substructure2 = Chem.MolFromSmarts(substructure2) 
    #substructure3 = Chem.MolFromSmarts(substructure3) 
    
    indices = molecule.GetSubstructMatches(substructure)
    #indices1 = molecule.GetSubstructMatches(substructure1)
    #indices2 = molecule.GetSubstructMatches(substructure2)
    #indices3 = molecule.GetSubstructMatches(substructure3)
    
    #return len(indices)+len(indices2)+len(indices1)+len(indices3)
    return len(indices)
    
def f_Ar_R(molecule: Chem.rdchem.Mol):
    ''' aromatic ring with R sostituent'''
    
    substructure = 'a-[CX4H2]-[C]-[!O;!N;!S]' 
    substructure1 = 'a-[CD4]-[CX4H3]' #C(CH3)3
    substructure2 = 'a-[CD3]-[CX4H3]' #CH(CH3)2
    substructure3 = 'a-[CX4H3]' #because aryl methyl exist
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    
    indices = molecule.GetSubstructMatches(substructure) 
    indices1 = molecule.GetSubstructMatches(substructure1) 
    indices2 = molecule.GetSubstructMatches(substructure2)
    indices3 = molecule.GetSubstructMatches(substructure3)
    
    #return (len(indices1)//3)+(len(indices2)//2)+len(indices)+len(indices3)
    return len(indices)+len(indices1)+len(indices2)+len(indices3)

def f_op_diphenolo_OR(molecule: Chem.rdchem.Mol):
    
    '''difenoli orto e para e sostituenti OR orto e para'''
    
    substructure = 'c1(-[OX2H])c(-[OX2H])cccc1'
    substructure1 = 'c1(-[OX2H])ccc(-[OX2H])cc1'
    substructure2 = 'c1(-[OX2]-[C^3])c(-[OX2]-[C^3])cccc1'
    substructure3 = 'c1(-[OX2]-[C^3])ccc(-[OX2]-[C^3])cc1'
    substructure4 = 'c1(-[OX2]-[C^3])c(-[OX2H])cccc1'
    substructure5 = 'c1(-[OX2H])ccc(-[OX2]-[C^3])cc1'
    
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
    
    a = len(indices)+len(indices1)+len(indices2)+len(indices3)+len(indices4)+len(indices5)
    return a

def f_Ar_COR(molecule: Chem.rdchem.Mol):
    '''Chetone as ring sostituen'''    
    substructure = 'a-[CX3](=O)-[#6]-[!O,!N,!S]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)     
    return len(indices)

def f_Ar_COO_R_H(molecule: Chem.rdchem.Mol):
    
    '''sostituente all'anello: estere 
       oppure acido benzoico'''
    
    substructure = 'c-C(=O)[O;H1,-]'
    substructure1 = 'c-C(=O)(-OC-[!O;!N;!S])'
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    
    indices = molecule.GetSubstructMatches(substructure)  
    indices1 = molecule.GetSubstructMatches(substructure1)
    return len(indices)+len(indices1)

def f_C_3phenols(molecule: Chem.rdchem.Mol):
    
    '''3 aromatic structures bound to C'''
    
    substructure = 'C(-a)(-a)(-a)'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return len(indices)

def f_Ar_Cl_Br(molecule: Chem.rdchem.Mol):
    
    '''Cl, Br bound to aromatic ring '''
    substructure = 'a[Cl,Br]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return len(indices)

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
    
    #2OH e 1 OR
    
    substructure10 = 'c1(-[OX2H])c(-[OX2H])c(-[OX2]-[C^3])ccc1'
    substructure11 = 'c1(-[OX2H])c(-[OX2H])cc(-[OX2]-[C^3])cc1'
    substructure12 = 'c1(-[OX2H])c(-[OX2]-[C^3])c(-[OX2H])ccc1'
    substructure13 = 'c1(-[OX2H])cc(-[OX2H])c(-[OX2]-[C^3])cc1'
    substructure14 = 'c1(-[OX2H])cc(-[OX2H])cc(-[OX2]-[C^3])c1'
    substructure15 = 'c1(-[OX2H])c(-[OX2]-[C^3])cc(-[OX2H])cc1'
    
    #2OR e 1 OH

    substructure16 = 'c1(-[OX2]-[C^3])c(-[OX2]-[C^3])c(-[OX2H])ccc1'
    substructure17 = 'c1(-[OX2]-[C^3])c(-[OX2]-[C^3])cc(-[OX2H])cc1'
    substructure18 = 'c1(-[OX2]-[C^3])c(-[OX2H])c(-[OX2]-[C^3])ccc1'
    substructure19 = 'c1(-[OX2]-[C^3])cc(-[OX2]-[C^3])c(-[OX2H])cc1'
    substructure20 = 'c1(-[OX2]-[C^3])cc(-[OX2]-[C^3])cc(-[OX2H])c1'
    substructure21 = 'c1(-[OX2]-[C^3])c(-[OX2H])cc(-[OX2]-[C^3])cc1'
    
     #3OR o 4OR
        
    substructure22 = 'c1(-[OX2]-[C^3])cc(-[OX2]-[C^3])cc(-[OX2]-[C^3])c1'
    substructure23 = 'c1(-[OX2]-[C^3])c(-[OX2]-[C^3])ccc(-[OX2]-[C^3])c1'
    substructure24 = 'c1(-[OX2]-[C^3])c(-[OX2]-[C^3])c(-[OX2]-[C^3])ccc1'
    substructure25 = 'c1(-[OX2]-[C^3])c(-[OX2]-[C^3])c(-[OX2]-[C^3])c(-[OX2]-[C^3])cc1'
    substructure26 = 'c1(-[OX2]-[C^3])c(-[OX2]-[C^3])c(-[OX2]-[C^3])cc(-[OX2]-[C^3])c1'   
    
    

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
    substructure10 = Chem.MolFromSmarts(substructure10)
    substructure11 = Chem.MolFromSmarts(substructure11)
    substructure12 = Chem.MolFromSmarts(substructure12)
    substructure13 = Chem.MolFromSmarts(substructure13)
    substructure14 = Chem.MolFromSmarts(substructure14)
    substructure15 = Chem.MolFromSmarts(substructure15)
    substructure16 = Chem.MolFromSmarts(substructure16)
    substructure17 = Chem.MolFromSmarts(substructure17)
    substructure18 = Chem.MolFromSmarts(substructure18)
    substructure19 = Chem.MolFromSmarts(substructure19)
    substructure20 = Chem.MolFromSmarts(substructure20)
    substructure21 = Chem.MolFromSmarts(substructure21)
    substructure22 = Chem.MolFromSmarts(substructure22)
    substructure23 = Chem.MolFromSmarts(substructure23)
    substructure24 = Chem.MolFromSmarts(substructure24)
    substructure25 = Chem.MolFromSmarts(substructure25)
    substructure26 = Chem.MolFromSmarts(substructure26)
    
    indices = molecule.GetSubstructMatches(substructure)
    indices1 = molecule.GetSubstructMatches(substructure1)
    indices2 = molecule.GetSubstructMatches(substructure2)
    indices3 = molecule.GetSubstructMatches(substructure3)
    indices4 = molecule.GetSubstructMatches(substructure4)
    indices5 = molecule.GetSubstructMatches(substructure5)
    indices6 = molecule.GetSubstructMatches(substructure6)
    indices7 = molecule.GetSubstructMatches(substructure7)
    indices8 = molecule.GetSubstructMatches(substructure8)
    indices9 = molecule.GetSubstructMatches(substructure9)
    indices10 = molecule.GetSubstructMatches(substructure10)
    indices11 = molecule.GetSubstructMatches(substructure11)
    indices12 = molecule.GetSubstructMatches(substructure12)
    indices13 = molecule.GetSubstructMatches(substructure13)
    indices14 = molecule.GetSubstructMatches(substructure14)
    indices15 = molecule.GetSubstructMatches(substructure15)
    indices16 = molecule.GetSubstructMatches(substructure16)
    indices17 = molecule.GetSubstructMatches(substructure17)
    indices18 = molecule.GetSubstructMatches(substructure18)
    indices19 = molecule.GetSubstructMatches(substructure19)
    indices20 = molecule.GetSubstructMatches(substructure20)
    indices21 = molecule.GetSubstructMatches(substructure21)
    indices22 = molecule.GetSubstructMatches(substructure22)
    indices23 = molecule.GetSubstructMatches(substructure23)
    indices24 = molecule.GetSubstructMatches(substructure24)
    indices25 = molecule.GetSubstructMatches(substructure25)
    indices26 = molecule.GetSubstructMatches(substructure26)
    
    a = len(indices)+len(indices1)+len(indices2)+len(indices3)+len(indices4)+len(indices5)+len(indices6)+len(indices7)+len(indices8)+len(indices9)+len(indices10)+len(indices11)+len(indices12)+len(indices13)+len(indices14)+len(indices15)+len(indices16)+len(indices17)+len(indices18)+len(indices19)+len(indices20)+len(indices21)+len(indices22)+len(indices23)+len(indices24)+len(indices25)+len(indices26)
    return a

def f_Sulfoxide(molecule: Chem.rdchem.Mol):
    ''' gruppo di sulfoxide '''
    
    substructure = '[#6]-S(=O)-[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return len(indices)

def f_CH2_Terminal(molecule: Chem.rdchem.Mol):
    ''' CH2 terminal '''
    
    substructure = '*~C=[CX3H2]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return len(indices)

def f_Alcool_1(molecule: Chem.rdchem.Mol):
    ''' primary OH'''
    
    substructure = '*-[CH2]-[OX2H]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return len(indices)

def f_triple(molecule: Chem.rdchem.Mol):
    '''triple bond'''
    
    substructure = '*#*'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return len(indices)


#03/01 sistemare ether rdkit

def f_ethere2(molecule: Chem.rdchem.Mol):
    '''ethere'''
    
    substructure = '[OD2](-[#6]-[#6])-[#6]-[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)

    substructure1= '[OD2](-[CH3])-[#6]-[#6]'
    substructure1= Chem.MolFromSmarts(substructure1)
    indices1= molecule.GetSubstructMatches(substructure1)
    
    return len(indices) + len(indices1)

# 10/01

def f_Ar_OSO2(molecule: Chem.rdchem.Mol):
    '''gruppo Ar-OS(=O)2'''
    
    substructure = 'a-O-S(~O)(~O)(*)'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return len(indices) 

def f_CC4(molecule: Chem.rdchem.Mol):
    '''gruppo CC4'''
    
    substructure = '[#6D4#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return len(indices)

#Erika aggiunge gruppi funzionali carini

def f_Imidothioesters(molecule: Chem.rdchem.Mol):
    '''Imidothioesters'''
    
    substructure = '[#1,#6][CX3](=[NX2][#1,#6])[SX2][C&!$([CX3]=[OX1,SX1,NX2])]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return len(indices)


def f_anhydrides(molecule: Chem.rdchem.Mol):
    '''anhydrides'''
    
    substructure = '[#6][#6X3](=[OX1])[#8X2][#6X3](=[OX1])[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return len(indices)

def f_carbamate(molecule: Chem.rdchem.Mol):
    '''Number of alkyl carbamates (subject to hydrolysis)'''
    
    substructure = '[$(N([!H])C(=O)[OH0])]C(=O)[OH0]C' 
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return len(indices)

#13/01

def f_aniline_term(molecule: Chem.rdchem.Mol):
    '''Number of aniline'''
    
    substructure = 'c1ccccc1-[NX3H2]' 
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return len(indices)

#25/01
def f_alkyl_halide(molecule: Chem.rdchem.Mol):
    '''Number of aniline'''
    
    substructure = '[C]-[Cl,Br,I,F]' 
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return len(indices)

def f_charge2(molecule: Chem.rdchem.Mol):
    '''charge +'''
    substructure = "[*+]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_charge1(molecule: Chem.rdchem.Mol):
    '''charge -'''
    substructure = '[*-]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_nitro_aliphatic(molecule: Chem.rdchem.Mol):
    '''nitro aliphatic'''
    substructure = 'N(=O)(O)-[C]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_ketone_aliphatic(molecule: Chem.rdchem.Mol):
    '''ketone aliphatic'''
    substructure = '[C]-[CX3](=O)-[C]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def f_Ar_ketone(molecule: Chem.rdchem.Mol):
    '''ketone aromatic'''
    substructure = '[c]-[CX3](=O)-[C]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
# 9/09
def f_ketone_dehydro(molecule: Chem.rdchem.Mol):
    '''ketone that can generate deihydro product with OH in beta e H in alpha'''
    substructure = '[C]-[C](=O)-[Ch1]-[C](-[Oh1])'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
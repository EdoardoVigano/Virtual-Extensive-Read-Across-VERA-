from rdkit import Chem
# NO rdkit

def f_steroid(molecule: Chem.rdchem.Mol):
    ''' steroid scaffold
        '''
    substructure = '[#6]12~[#6]~[#6]~[#6]3~[#6](~[#6]1~[#6]~[#6]~[#6]2)~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~[#6]34'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices
    
def f_five_ring_hetero(molecule: Chem.rdchem.Mol):
    ''' five_ring_hetero
        '''
    substructure = '[#6]1~[#6]~[#6]~[#6]~[!C;!c]1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)   
    return indices
    
def f_six_ring_hetero(molecule: Chem.rdchem.Mol):
    ''' six_ring_hetero
        '''
    substructure = '[#6]1~[#6]~[#6]~[#6]~[#6]~[!C;!c]1'
    substructure = Chem.MolFromSmarts(substructure)
    # find substructures
    indices = molecule.GetSubstructMatches(substructure)    
    return indices
    
def f_benzCH2(molecule: Chem.rdchem.Mol):
    ''' benzCH2 present in molecule
        '''
    substructure = 'c1c([C&H2])cccc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)    
    return indices
    
def f_heteroatoms_heterocycles6(molecule: Chem.rdchem.Mol):
    ''' heteroatoms (O or S) bonded to a 6 membered
        heterocyclic compound
        '''  
    all_indices = []  
    substructures = ['[O,S]-[#6]1~[#6]~[#6]~[#6]~[#6]~[!C;!c]1',
                    '[#6]1~[#6]([O,S])~[#6]~[#6]~[#6]~[!C;!c]1',
                    '[#6]1~[#6]~[#6]([O,S])~[#6]~[#6]~[!C;!c]1']
    
    substructures = map(Chem.MolFromSmarts, substructures)
    indices = []
    for fragment in substructures:
        a = molecule.GetSubstructMatches(fragment)
        if len(a)>0:
            indices.append(molecule.GetSubstructMatches(fragment))          
    return indices

def f_heteroatoms_heterocycles5(molecule: Chem.rdchem.Mol):
    ''' heteroatoms (O or S) bonded to a 5 membered
        heterocyclic compound
        '''  
    all_indices = []  
    substructures = ['[#6]1~[#6]([O,S])~[#6]~[#6]~[!C;!c]1',
                    '[#6]1~[#6]~[#6]~[#6]([O,S])~[!C;!c]1']
    
    substructures = map(Chem.MolFromSmarts, substructures)
    indices = []
    for fragment in substructures:
        a = molecule.GetSubstructMatches(fragment)
        if len(a)>0:
            indices.append(molecule.GetSubstructMatches(fragment))          
    return indices

def f_biphenol(molecule: Chem.rdchem.Mol):
    ''' identify catechol, resorcinol,
        and hydroquinone
        ''' 
    all_indices = []   
    substructures = ['[OH]c1c([OH])cccc1',
                    '[OH]c1cc([OH])ccc1',
                    '[OH]c1ccc([OH])cc1']  
    substructures = map(Chem.MolFromSmarts, substructures)
    indices = []
    for fragment in substructures:
        a = molecule.GetSubstructMatches(fragment)
        if len(a)>0:
            indices.append(molecule.GetSubstructMatches(fragment))          
    return indices

def f_benzaldehyde(molecule: Chem.rdchem.Mol):
    ''' benzaldehyde present in molecule
        ''' 
    substructure = 'c1ccccc1-[C;H1](=O)'
    substructure = Chem.MolFromSmarts(substructure)
    # find substructures
    indices = molecule.GetSubstructMatches(substructure)     
    return indices

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
    
    indices0 = molecule.GetSubstructMatches(substructure)
    indices1 = molecule.GetSubstructMatches(substructure1)
    indices2 = molecule.GetSubstructMatches(substructure2)
    indices3 = molecule.GetSubstructMatches(substructure3)
    
    indices = [indices0, indices1, indices2, indices3]
    
    return indices
    
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
    
    indices0 = molecule.GetSubstructMatches(substructure) 
    indices1 = molecule.GetSubstructMatches(substructure1) 
    indices2 = molecule.GetSubstructMatches(substructure2)
    indices3 = molecule.GetSubstructMatches(substructure3)
    indices = [indices0, indices1, indices2, indices3]
    
    tot=[]
    for ind in indices:
        if len(ind)>0:
            tot.append(ind)
    
    return tot

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
    
    indices0 = molecule.GetSubstructMatches(substructure)
    indices1 = molecule.GetSubstructMatches(substructure1)
    indices2 = molecule.GetSubstructMatches(substructure2)
    indices3 = molecule.GetSubstructMatches(substructure3)
    indices4 = molecule.GetSubstructMatches(substructure4)
    indices5 = molecule.GetSubstructMatches(substructure5)
    
    indices = [indices0,indices1,indices2,indices3,indices4,indices5]
    
    tot=[]
    for ind in indices:
        if len(ind)>0:
            tot.append(ind)
    
    
    return tot

def f_Ar_COR(molecule: Chem.rdchem.Mol):
    '''Chetone as ring sostituen'''    
    substructure = 'a-[CX3](=O)-[#6]-[!O,!N,!S]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)     
    return indices

def f_Ar_COO_R_H(molecule: Chem.rdchem.Mol):
    
    '''sostituente all'anello: estere 
       oppure acido benzoico'''
    
    substructure = 'c-C(=O)[O;H1,-]'
    substructure1 = 'c-C(=O)(-OC-[!O;!N;!S])'
    
    substructure = Chem.MolFromSmarts(substructure)
    substructure1 = Chem.MolFromSmarts(substructure1)
    
    indices0 = molecule.GetSubstructMatches(substructure)  
    indices1 = molecule.GetSubstructMatches(substructure1)
    
    indices = [indices0,indices1]
    tot=[]
    for ind in indices:
        if len(ind)>0:
            tot.append(ind)   
    return tot

def f_C_3phenols(molecule: Chem.rdchem.Mol):
    
    '''3 phenols bound to C'''
    
    substructure = 'C(-a)(-a)(-a)'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return indices

def f_Ar_Cl_Br(molecule: Chem.rdchem.Mol):
    
    '''Cl, Br bound to aromatic ring '''
    substructure = 'a[Cl,Br]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
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
    
    indices0 = molecule.GetSubstructMatches(substructure)
    indices1 = molecule.GetSubstructMatches(substructure1)
    indices2 = molecule.GetSubstructMatches(substructure2)
    indices3 = molecule.GetSubstructMatches(substructure3)
    indices4 = molecule.GetSubstructMatches(substructure4)
    indices5 = molecule.GetSubstructMatches(substructure5)
    indices6 = molecule.GetSubstructMatches(substructure6)
    indices7 = molecule.GetSubstructMatches(substructure7)
    indices8 = molecule.GetSubstructMatches(substructure8)
    indices9 = molecule.GetSubstructMatches(substructure9)
    
    indices = [indices0,indices1,indices2,indices3,indices4,indices6,indices7,indices8,indices9]
    tot = []
    for ind in indices:
        if len(ind)>0:
            tot.append(ind)   
    return tot

def f_Sulfoxide(molecule: Chem.rdchem.Mol):
    ''' gruppo di sulfoxide '''
    
    substructure = '[#6]-S(=O)-[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return indices

def f_CH2_Terminal(molecule: Chem.rdchem.Mol):
    ''' CH2 terminal '''
    
    substructure = '*~C=[CX3H2]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return indices

def f_Alcool_1(molecule: Chem.rdchem.Mol):
    ''' primary OH'''
    
    substructure = '*-[CH2]-[OX2H]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return indices

def f_triple(molecule: Chem.rdchem.Mol):
    '''triple bond'''
    
    substructure = '*#*'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return indices


#03/01 sistemare ether rdkit

def f_ethere2(molecule: Chem.rdchem.Mol):
    '''ethere'''
    
    substructure = '[OD2](-[#6]-[#6])-[#6]-[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices0 = molecule.GetSubstructMatches(substructure)

    substructure1= '[OD2](-[CH3])-[#6]-[#6]'
    substructure1= Chem.MolFromSmarts(substructure1)
    indices1= molecule.GetSubstructMatches(substructure1)
    
    indices = [indices0,indices1]
    tot = []
    for ind in indices:
        if len(ind)>0:
            tot.append(ind)   
    return tot

# 10/01

def f_Ar_OSO2(molecule: Chem.rdchem.Mol):
    '''gruppo Ar-OS(=O)2'''
    
    substructure = 'a-O-S(~O)(~O)(*)'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return indices 

def f_CC4(molecule: Chem.rdchem.Mol):
    '''gruppo CC4'''
    
    substructure = '[#6D4#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return indices

#13/01

def f_aniline_term(molecule: Chem.rdchem.Mol):
    ''' aniline groups
        '''
    substructure = 'c-[NX3H2]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

#FROM RDKIT

def f_benzene(molecule: Chem.rdchem.Mol):
    
    '''benzene '''
    substructure = 'c1ccccc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return indices

def f_phenol(molecule: Chem.rdchem.Mol):
    
    '''phenol '''
    substructure = 'c1ccccc1O'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    
    return indices
'''
def f_Imidothioesters(molecule: Chem.rdchem.Mol):
 
    substructure = '[#1,#6][CX3](=[NX2][#1,#6])[SX2][C&!$([CX3]=[OX1,SX1,NX2])]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return indices


def f_anhydrides(molecule: Chem.rdchem.Mol):
    
    substructure = '[#6][#6X3](=[OX1])[#8X2][#6X3](=[OX1])[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return indices'''

#Gruppi di RDKit da usare per i duplicati

def f_Ar_OH(molecule: Chem.rdchem.Mol):
    '''ArOH'''
    
    substructure = 'c[OH1]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return indices

def f_benzodiazepine(molecule: Chem.rdchem.Mol):
    '''benzodiazepine'''
    
    substructure = '[c&R2]12[c&R1][c&R1][c&R1][c&R1][c&R2]1[N&R1][C&R1][C&R1][N&R1]=[C&R1]2'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return indices

def f_aryl_methyl(molecule: Chem.rdchem.Mol):
    '''aryl_methyl'''
    
    substructure = '[$(a-[CH3]),$(a-[CH2]-[CH3]),$(a-[CH2]-[CH2]~[!N;!O]);!$(a(:a!:*):a!:*)]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)
    return indices

def f_piperdine(molecule: Chem.rdchem.Mol):
    ''' piperdine
        '''
    substructure = 'N1CCCCC1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_urea(molecule: Chem.rdchem.Mol):
    ''' urea
        '''
    substructure = 'C(=O)(-N)-N'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_amide(molecule: Chem.rdchem.Mol):
    ''' amide
        '''
    substructure = 'C(=O)-N'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_aniline(molecule: Chem.rdchem.Mol):
    ''' aniline groups
        '''
    substructure = 'c-[NX3;!$(N=*)]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices


def f_NH2(molecule: Chem.rdchem.Mol):
    ''' "Number of Primary amines"s
        '''
    substructure = '[NH2,nH2]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_NH1(molecule: Chem.rdchem.Mol):
    ''' "Number of Secondary amines"
        '''
    substructure = '[NH1,nH1]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_NH0(molecule: Chem.rdchem.Mol):
    ''' "Number of Tertiary amines"
        '''
    substructure = '[NH0,nH0]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices


def f_nitro(molecule: Chem.rdchem.Mol):
    ''' "nitro"
        '''
    substructure = 'N(=O)(O)[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_nitro_arom_nonortho(molecule: Chem.rdchem.Mol):
    ''' nitro_arom_nonortho
        '''
    substructure = '[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1);!$(cc-!:*)]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_pyridine(molecule: Chem.rdchem.Mol):
    ''' "Number of pyridine rings"
        '''
    substructure = 'n1ccccc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_alkyl_carbamate(molecule: Chem.rdchem.Mol):
    ''' "Number of alkyl carbamates (subject to hydrolysis)"
        '''
    substructure = 'C[NH1]C(=O)OC'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_priamide(molecule: Chem.rdchem.Mol):
    ''' Number of primary amides
        '''
    substructure = 'C(=O)-[NH2]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_ketone_Topliss(molecule: Chem.rdchem.Mol):
    ''' "Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha"
        '''
    substructure = '[$([CX3](=[OX1])(C)([c,C]));!$([CX3](=[OX1])([CH1]=C)[c,C])]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_ketone(molecule: Chem.rdchem.Mol): 
    ''' "Number of ketones"
        '''
    substructure = '[#6][CX3](=O)[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_priamide(molecule: Chem.rdchem.Mol):
    ''' Number of primary amides
        '''
    substructure = 'C(=O)-[NH2]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_hdrzine(molecule: Chem.rdchem.Mol):
    ''' "Number of hydrazine groups"
        '''
    substructure = '[NX3]-[NX3]' #io farei [NX3H]-[NX3H2]
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_ester(molecule: Chem.rdchem.Mol):
    ''' "Number of esters"
        '''
    substructure = '[#6][CX3](=O)[OX2H0][#6]' 
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_ether(molecule: Chem.rdchem.Mol):
    '''"Number of ether oxygens (including phenoxy)"
        '''
    substructure = '[OD2]([#6])[#6]' 
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_carbamate(molecule: Chem.rdchem.Mol):
    '''carbamate
        '''
    substructure = '[$(N([!H])C(=O)[OH0])]C(=O)[OH0]C' 
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_charge2(molecule: Chem.rdchem.Mol):
    '''charge +'''
    substructure = "[*+]"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_charge1(molecule: Chem.rdchem.Mol):
    '''charge -'''
    substructure = '[*-]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_imidazole(molecule: Chem.rdchem.Mol):
    "Number of imidazole rings"
    substructure = 'n1cncc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

# aggiunte 22/07

def f_Imidothioesters(molecule: Chem.rdchem.Mol):
    "Number of Imidothioesters"
    substructure = '[#1,#6][CX3](=[NX2][#1,#6])[SX2][C&!$([CX3]=[OX1,SX1,NX2])]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_anhydrides(molecule: Chem.rdchem.Mol):
    '''anhydrides'''
    
    substructure = '[#6][#6X3](=[OX1])[#8X2][#6X3](=[OX1])[#6]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_alkyl_halide(molecule: Chem.rdchem.Mol):
    '''Number of aniline'''
    
    substructure = '[C]-[Cl,Br,I,F]' 
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_nitro_aliphatic(molecule: Chem.rdchem.Mol):
    '''nitro aliphatic'''
    substructure = 'N(=O)(O)-[C]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_ketone_aliphatic(molecule: Chem.rdchem.Mol):
    '''ketone aliphatic'''
    substructure = '[C]-[CX3](=O)-[C]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

def f_Ar_ketone(molecule: Chem.rdchem.Mol):
    '''ketone aromatic'''
    substructure = '[c]-[CX3](=O)-[C]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices

# 9/09
def f_ketone_dehydro(molecule: Chem.rdchem.Mol):
    '''ketone that can generate deihydro product with OH in beta e H in alpha'''
    substructure = '[C]-[C](=O)-[Ch1]-[C](-[Oh1])'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return indices



import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments

'''Note:
Every time you add something here you have to add the same changes in thesefile.py: 

 1.alerts_index
 2.calculator 
 3.alerts
'''

def fish1(molecule: Chem.rdchem.Mol):
    ''' FISH 1
        '''
    substructure = 'O(c1ccccc1)c1ccccc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    
def fish2(molecule: Chem.rdchem.Mol):
    ''' FISH 2
        '''
    substructure = '[O-][N+](=O)c1ccccc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish3(molecule: Chem.rdchem.Mol):
    ''' FISH 3
        '''
    substructure = 'C[S]c1ccccc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish4(molecule: Chem.rdchem.Mol):
    ''' FISH 4
        '''
    substructure = '*O[CX4H2]C#N'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish5(molecule: Chem.rdchem.Mol):
    ''' FISH 5
        '''
    substructure = '[*;A][N+]([*;A])([*;A])[*;A]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish6(molecule: Chem.rdchem.Mol):
    ''' FISH 6
        '''
    substructure = '*SS*'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish7(molecule: Chem.rdchem.Mol):
    ''' FISH 7
        '''
    substructure = 'Cl[C]=[C;A]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish8(molecule: Chem.rdchem.Mol):
    ''' FISH 8
        '''
    substructure = 'Cn1cc(C(N)=O)c(C)n1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish9(molecule: Chem.rdchem.Mol):
    ''' FISH 9
        '''
    substructure = '*P(*)(*)=*'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

#modifiche del 9 settembre da 10 a 21

def fish10(molecule: Chem.rdchem.Mol):
    ''' FISH 10: Nucleophilic substitution: allylic activation
        '''
    substructure = 'C=C[CX4H2][Br,I,Cl]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish11(molecule: Chem.rdchem.Mol):
    ''' FISH 11: Nucleophilic substitution: propargylyc activation
        '''
    substructure = 'C#C[CX4H2][Br,I,Cl]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish12(molecule: Chem.rdchem.Mol):
    ''' FISH 12: Nucleophilic substitution: benzylic activation
        '''
    substructure = '[S]([S]c1ccccc1)c1ccccc1'
    substructure2 = 'c1ccccc1[CX4H2][Cl,Br,I]'
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices = molecule.GetSubstructMatches(substructure)   
    indices2 = molecule.GetSubstructMatches(substructure2)
    return len(indices)+len(indices2)

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
    
    indices = molecule.GetSubstructMatches(substructure)
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
    return len(indices)+len(indices2)+len(indices3)+len(indices4)+len(indices5)+len(indices6)+len(indices7)+len(indices8)+len(indices9)+len(indices10)+len(indices11)+len(indices12)+len(indices13)+len(indices14)

def fish14(molecule: Chem.rdchem.Mol):
    ''' FISH 14 : Acid anhydrides (acylation)
        '''
    substructure = 'CC(=O)OC(=O)C'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish15(molecule: Chem.rdchem.Mol):
    ''' FISH 15: Schiff base
        '''
	# [*]CC=O
    substructure = '[c,C][CX3H1]=O'
    substructure2 = "c1(F)c(F)c(F)c(F)c(F)c1C(=O)"
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure) 
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices2 = molecule.GetSubstructMatches(substructure2)
    return len(indices)+ len(indices2)

def fish16(molecule: Chem.rdchem.Mol):
    ''' FISH 16: Michael Addition
        '''
    substructure = '[Ch1]=[Ch1]C(=O)C'#chetone alpha-beta insaturo
    substructure2 = "[Ch1]=[Ch1]C(=O)[NX3H2]" #amide alpha-beta instatura
    substructure3 = "[Ch1]=[Ch1][CX3H1](=O)" #aldeide alpha-beta insatura
    substructure4 = "[Ch1]=[Ch1]C(=O)C#N" #alpha beta insaturo con C#N
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    substructure4 = Chem.MolFromSmarts(substructure4)
    indices = molecule.GetSubstructMatches(substructure)
    indices2 = molecule.GetSubstructMatches(substructure) 
    indices3 = molecule.GetSubstructMatches(substructure) 
    indices4 = molecule.GetSubstructMatches(substructure)
    return len(indices) + len(indices2) + len(indices3) + len(indices4)
'''
def fish17(molecule: Chem.rdchem.Mol):
     #FISH 17: Michael Addition - è compreso nel 16
        
    substructure = 'C=CC#N'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
'''

def fish18(molecule: Chem.rdchem.Mol):
    ''' FISH 18: Strained three-membered heterocyclic ring
        '''
    substructure = 'C1[O,N]C1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish19(molecule: Chem.rdchem.Mol):
    ''' FISH 19: Proelectrophiles - alcohol dehydrogenase activators
        '''
	# OCC=C
    substructure = '[Ch]=[Ch][CX4H2][OX2H1]' #con alcool primario
    substructure2 = "[Ch]=[Ch][CX4H1][OX2H1]" #con alcool secondario - la reazione non avviene con alcoli terziari
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices = molecule.GetSubstructMatches(substructure) 
    indices2 = molecule.GetSubstructMatches(substructure2) 
    return len(indices) + len(indices2)

def fish20(molecule: Chem.rdchem.Mol):
    ''' FISH 20: Proelectrophiles - Glutathione transferase activator
        '''
    substructure = 'BrCCCBr'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish21(molecule: Chem.rdchem.Mol):
    ''' FISH 21: Proelectrophiles - alcohol dehydrogenase activators
        '''
	# 'OCC#C'
    substructure = 'C#C[CX4H2][OX2H1]'
    substructure2 = 'C#C[CX4H1][OX2H1]'
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices = molecule.GetSubstructMatches(substructure)
    indices2 = molecule.GetSubstructMatches(substructure2)
    return len(indices) + len(indices2)

def fish22(molecule: Chem.rdchem.Mol):
    ''' FISH 22
        '''
    substructure = 'C1(C(C1C)(C)C)C(O)=O'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish23(molecule: Chem.rdchem.Mol):
    ''' FISH 23
        '''
    substructure = 'O=C(O)C(C)C'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish24(molecule: Chem.rdchem.Mol):
    ''' FISH 24
        '''
    substructure = 'C1(=CC=C(C=C1)O)[Cl]'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
'''
def fish25(molecule: Chem.rdchem.Mol):
    #FISH 25 alert di non tossicità
        
    substructure = 'O-C'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish26(molecule: Chem.rdchem.Mol):
     #FISH 26- Erika: ho aggiunto * prima della N altrimenti era una molecola a sè
        
        #alert di non tossicità
    substructure = 'NC=ONC=O'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    '''
    
    #Aggiunte dal paper

def fish27(molecule: Chem.rdchem.Mol):
    ''' FISH 27: Cyanogenic Mechanisms - cyanohydrin (rilasciano formaldeide e/o gruppo ciano)
        '''
    substructure = '*C[CX4H1]([OX2H1])C#N'
    substructure2 = "[Cl]c1ccc([NX3H1][CX4H2]C#N)cc1"
    substructure3 = "CC(=O)O[CX4H2]C#N"
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices = molecule.GetSubstructMatches(substructure)
    indices2 = molecule.GetSubstructMatches(substructure2) 
    indices3 = molecule.GetSubstructMatches(substructure3) 
    return len(indices)+ len(indices2) + len(indices3)

def fish28(molecule: Chem.rdchem.Mol):
    ''' FISH 28: Cyanogenic Mechanisms - cyanide
        '''
    substructure = '[Ch1]=[Ch1]C#N'
    substructure2 = "N#C[CX4H2]C#N"
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    indices = molecule.GetSubstructMatches(substructure)
    indices2 = molecule.GetSubstructMatches(substructure2) 
    return len(indices)+ len(indices2)

def fish29(molecule: Chem.rdchem.Mol):
    ''' FISH 29: Multistep or Multiple mechanisms - thiocyanates (può fare addizione di Michael sul C del ciano oppure attacco nucleofilo su S)
        '''
    substructure = '[C,c]SC#N'
    substructure2 = "[NX3](=O)(O)aaaa(SC#N)"
    substructure3 = "[NX3H2]aaaa(SC#N)"
    substructure = Chem.MolFromSmarts(substructure)
    substructure2 = Chem.MolFromSmarts(substructure2)
    substructure3 = Chem.MolFromSmarts(substructure3)
    indices = molecule.GetSubstructMatches(substructure)
    indices2 = molecule.GetSubstructMatches(substructure2) 
    indices3 = molecule.GetSubstructMatches(substructure3)
    return len(indices)+ len(indices2)+ len(indices3)

def fish30(molecule: Chem.rdchem.Mol):
    ''' FISH 30: Multistep or Multiple mechanisms - alpha-chloro benzyldene malononitrile (addizione di Michael su doppio legame o addizione di acqua e formazione di base di schiff)
        '''
    substructure = '[Cl]c1c([Ch1]=C(C#N)(C#N))cccc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)

def fish31(molecule: Chem.rdchem.Mol):
    ''' FISH 31: Multistep or Multiple mechanisms - Tautomerization
        '''
    substructure = 'O=Nc1ccc(N(C)(C))cc1'
    substructure = Chem.MolFromSmarts(substructure)
    indices = molecule.GetSubstructMatches(substructure)       
    return len(indices)
    

def ToxRead0(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'O=C(OCc1cccc(Oc2ccccc2)c1)C3CC3'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def ToxRead1(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'C1=CC2CC1C3CCCC23'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def ToxRead2(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'O=C(OCc1ccccc1)C2CC2'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def ToxRead3(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'O=C(OC)C1C(C=C)C1(C)C'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def ToxRead4(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = '[#7]-1-[#16]-c2ccccc2-[#6]-1=O'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def ToxRead5(molecule: Chem.rdchem.Mol):
	# Alerts Toxic  (< 1mg/l)
	substructure = 'N#Cc1cn(c(c1))'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def ToxRead6(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Cc2ccccc2)cc1'
	substructure2  = 'Cl-C(-Cl)(-Cl)-C(-c1:c:c:c:c:c:1)-c2:c:c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead8(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Cc2ccccc2)cc1'
	substructure2  = 'O-c1:c(:c:c(:c:c:1)-Cl)-C-c2:c:c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead10(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N#C-C(-O-C(-C)=O)-c1:c:c(:c:c:c:1)-O-c2:c:c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead12(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N#C-C(-O-C(-C1-C(-C-1-C)(-C)-C)=O)-c2:c:c(:c:c:c:2)-O-c3:c:c:c:c:c:3'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead14(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N#C-C(-O-C(-C1-C(-C-1-C=C)(-C)-C)=O)-c2:c:c(:c:c:c:2)-O-c3:c:c:c:c:c:3'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead16(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'O=C(Nc1ccccc1)c2ccccc2'
	substructure2  = 'Cl-c1:c:c(:c:c:c:1)-N'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead18(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'C1=CCCCC1'
	substructure2  = 'Cl-C12-C3-C(-C(-C-2(-Cl)-Cl)(-C(=C-1-Cl)-Cl)-Cl)-C-C-C-3'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead20(molecule: Chem.rdchem.Mol):
	# coppie Toxic alerts 
	substructure = 'C1=CC2CCC1C2'
	substructure2  = 'Cl-C12-C3-C(-C(-C-2(-Cl)-Cl)(-C(=C-1-Cl)-Cl)-Cl)-C-C-C-3'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead22(molecule: Chem.rdchem.Mol):
	# Alerts NoTOx  ( > 100 mg/l)
	substructure = 'O=C1NC=CC=C1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def ToxRead23(molecule: Chem.rdchem.Mol):
	# coppie NoTOx alerts 
	substructure = 'c1ccncc1'
	substructure2  = 'O-c1:c:n:c:c:c:1'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead25(molecule: Chem.rdchem.Mol):
	# coppie NoTOx alerts 
	substructure = 'c1ccncc1'
	substructure2  = 'O=C1-N-C=C-C=C-1'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead27(molecule: Chem.rdchem.Mol):
	# coppie NoTOx alerts 
	substructure = 'c1ccncc1'
	substructure2  = 'n1:c:c:c(:c:c:1)-C-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead29(molecule: Chem.rdchem.Mol):
	# coppie NoTOx alerts 
	substructure = 'c1cncnc1'
	substructure2  = 'O=C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead31(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc2ccccc2c1'
	substructure2  = 'O'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead33(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc2ccccc2c1'
	substructure2  = 'O-c1:c2:c(:c:c:c:1):c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead35(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc2ccccc2c1'
	substructure2  = 'O-c1:c:c:c:c:c:1'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead37(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead39(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'N-c1:c:c:c:c:c:1'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead41(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'c1ccc(Oc2ccccc2)cc1'
	substructure2  = 'O(-c1:c:c:c(:c:c:1)-O)-c2:c:c:c:c:c:2'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead43(molecule: Chem.rdchem.Mol):
	# coppie Alerts Toxic 1-10 mg/l 
	substructure = 'O=C1NC(=O)c2ccccc21'
	substructure2  = 'N1(-C(-c2:c(-C-1=O):c:c:c:c:2)=O)-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead45(molecule: Chem.rdchem.Mol):
	# coppie Alerts 10-100 mg/l Toxic
	substructure = 'C1CCCCC1'
	substructure2  = 'N(-C1-C-C-C-C-C-1)-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead47(molecule: Chem.rdchem.Mol):
	# coppie Alerts 10-100 mg/l Toxic
	substructure = 'C1CCCCC1'
	substructure2  = 'N-C1-C-C-C(-C-C-1)-C-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead49(molecule: Chem.rdchem.Mol):
	# coppie Alerts 10-100 mg/l Toxic
	substructure = 'C1CCCCC1'
	substructure2  = 'O=C1-C(-C-C-C-C-1)-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

def ToxRead51(molecule: Chem.rdchem.Mol):
	# coppie Alerts 10-100 mg/l Toxic
	substructure = 'O=C1CCCCC1'
	substructure2  = 'O=C1-C(-C-C-C-C-1)-C'
	substructure = Chem.MolFromSmarts(substructure)
	substructure2 = Chem.MolFromSmarts(substructure2)
	indices = molecule.GetSubstructMatches(substructure)
	indices2 = molecule.GetSubstructMatches(substructure2)
	if  len(indices)>0 and len(indices2)>0: out = 1
	else: out = 0
	return out

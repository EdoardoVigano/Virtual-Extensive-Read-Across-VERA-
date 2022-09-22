import Fish_Alerts as Fish
import pandas as pd
from rdkit import Chem

'''Note:
Every time you add something here you have to add the functions in file.py alerts_index, calculator and alerts
'''

def descriptor(mol: Chem.rdchem.Mol):
    fish1=[]; fish2=[]; fish3=[]; fish4=[]; fish5=[]; fish6=[]; fish7=[]; fish8=[]; fish9=[]; fish10=[]; fish11=[]; fish12=[]; fish13=[]; fish14=[]; fish15=[]; fish16=[]; fish18=[]; fish19=[]; fish20=[]; fish21=[]; fish22=[]; fish23=[]; fish24=[];fish27=[]; fish28=[];fish29=[]; fish30=[]; fish31=[]; ToxRead0=[];	ToxRead1=[];	ToxRead2=[];	ToxRead3=[];	ToxRead4=[];	ToxRead5=[];	ToxRead6=[];	ToxRead8=[];	ToxRead10=[];	ToxRead12=[];	ToxRead14=[];	ToxRead16=[];	ToxRead18=[];	ToxRead20=[];	ToxRead22=[];	ToxRead23=[];	ToxRead25=[];	ToxRead27=[];	ToxRead29=[];	ToxRead31=[];	ToxRead33=[];	ToxRead35=[];	ToxRead37=[];	ToxRead39=[];	ToxRead41=[];	ToxRead43=[];	ToxRead45=[];	ToxRead47=[];	ToxRead49=[];	ToxRead51=[]


    lst_chemical_fish=[fish1, fish2, fish3, fish4, fish5, fish6, fish7, fish8, fish9, fish10, fish11, fish12, fish13, fish14, fish15, fish16, fish18, fish19, fish20, fish21, fish22, fish23, fish24, fish27, fish28, fish29, fish30, fish31, ToxRead0,	ToxRead1,	ToxRead2,	ToxRead3,	ToxRead4,	ToxRead5,	ToxRead6,	ToxRead8,	ToxRead10,	ToxRead12,	ToxRead14,	ToxRead16,	ToxRead18,	ToxRead20,	ToxRead22,	ToxRead23,	ToxRead25,	ToxRead27,	ToxRead29,	ToxRead31,	ToxRead33,	ToxRead35,	ToxRead37,	ToxRead39,	ToxRead41,	ToxRead43,	ToxRead45,	ToxRead47,	ToxRead49,	ToxRead51] 
    
    if mol != 'none':    
        fish1.append(Fish.fish1(mol))
        fish2.append(Fish.fish2(mol))
        fish3.append(Fish.fish3(mol))
        fish4.append(Fish.fish4(mol))
        fish5.append(Fish.fish5(mol))
        fish6.append(Fish.fish6(mol))
        fish7.append(Fish.fish7(mol))
        fish8.append(Fish.fish8(mol))
        fish9.append(Fish.fish9(mol))
        fish10.append(Fish.fish10(mol))
        fish11.append(Fish.fish11(mol))
        fish12.append(Fish.fish12(mol))
        fish13.append(Fish.fish13(mol))
        fish14.append(Fish.fish14(mol))
        fish15.append(Fish.fish15(mol))
        fish16.append(Fish.fish16(mol))
        #fish17.append(Fish.fish17(mol))
        fish18.append(Fish.fish18(mol))
        fish19.append(Fish.fish19(mol))
        fish20.append(Fish.fish20(mol))
        fish21.append(Fish.fish21(mol))
        fish22.append(Fish.fish22(mol))
        fish23.append(Fish.fish23(mol))
        fish24.append(Fish.fish24(mol))
        #fish26.append(Fish.fish26(mol))
        fish27.append(Fish.fish27(mol))
        fish28.append(Fish.fish28(mol))
        fish29.append(Fish.fish29(mol))
        fish30.append(Fish.fish30(mol))
        fish31.append(Fish.fish31(mol))
        ToxRead0.append(Fish.ToxRead0(mol))
        ToxRead1.append(Fish.ToxRead1(mol))
        ToxRead2.append(Fish.ToxRead2(mol))
        ToxRead3.append(Fish.ToxRead3(mol))
        ToxRead4.append(Fish.ToxRead4(mol))
        ToxRead5.append(Fish.ToxRead5(mol))
        ToxRead6.append(Fish.ToxRead6(mol))
        ToxRead8.append(Fish.ToxRead8(mol))
        ToxRead10.append(Fish.ToxRead10(mol))
        ToxRead12.append(Fish.ToxRead12(mol))
        ToxRead14.append(Fish.ToxRead14(mol))
        ToxRead16.append(Fish.ToxRead16(mol))
        ToxRead18.append(Fish.ToxRead18(mol))
        ToxRead20.append(Fish.ToxRead20(mol))
        ToxRead22.append(Fish.ToxRead22(mol))
        ToxRead23.append(Fish.ToxRead23(mol))
        ToxRead25.append(Fish.ToxRead25(mol))
        ToxRead27.append(Fish.ToxRead27(mol))
        ToxRead29.append(Fish.ToxRead29(mol))
        ToxRead31.append(Fish.ToxRead31(mol))
        ToxRead33.append(Fish.ToxRead33(mol))
        ToxRead35.append(Fish.ToxRead35(mol))
        ToxRead37.append(Fish.ToxRead37(mol))
        ToxRead39.append(Fish.ToxRead39(mol))
        ToxRead41.append(Fish.ToxRead41(mol))
        ToxRead43.append(Fish.ToxRead43(mol))
        ToxRead45.append(Fish.ToxRead45(mol))
        ToxRead47.append(Fish.ToxRead47(mol))
        ToxRead49.append(Fish.ToxRead49(mol))
        ToxRead51.append(Fish.ToxRead51(mol))
    else:
        for chem_group in lst_chemical_fish:
            chem_group.append('none')
       
    out = pd.DataFrame(list(zip(fish1, fish2, fish3, fish4, fish5, fish6, fish7, fish8, fish9, fish10, fish11, fish12, fish13, fish14, fish15, fish16, fish18, fish19, fish20, fish21, fish22, fish23, fish24,fish27, fish28, fish29, fish30, fish31, ToxRead0, ToxRead1, ToxRead2, ToxRead3, ToxRead4, ToxRead5, ToxRead6, ToxRead8, ToxRead10, ToxRead12, ToxRead14,	ToxRead16, ToxRead18, ToxRead20, ToxRead22,	ToxRead23,	ToxRead25,	ToxRead27,	ToxRead29,	ToxRead31,	ToxRead33,	ToxRead35,	ToxRead37,	ToxRead39,	ToxRead41,	ToxRead43,	ToxRead45,	ToxRead47,	ToxRead49,	ToxRead51)),
    columns=["fish1", "fish2", "fish3", "fish4", "fish5", "fish6", "fish7", "fish8", "fish9", "fish10", "fish11", "fish12", "fish13", "fish14", "fish15", "fish16", "fish18", "fish19", "fish20", "fish21", "fish22", "fish23", "fish24", "fish27", "fish28", "fish29", "fish30", "fish31", "ToxRead0",	"ToxRead1",	"ToxRead2",	"ToxRead3",	"ToxRead4",	"ToxRead5",	"ToxRead6",	"ToxRead8",	"ToxRead10", "ToxRead12", "ToxRead14", "ToxRead16",	"ToxRead18", "ToxRead20", "ToxRead22", "ToxRead23",	"ToxRead25", "ToxRead27", "ToxRead29", "ToxRead31",	"ToxRead33", "ToxRead35", "ToxRead37", "ToxRead39",	"ToxRead41", "ToxRead43", "ToxRead45", "ToxRead47",	"ToxRead49", "ToxRead51"])
    return out

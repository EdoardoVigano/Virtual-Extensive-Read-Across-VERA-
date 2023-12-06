import pandas as pd
from rdkit import Chem
import Benigni_Bossa as BB

def descriptor(mol: Chem.rdchem.Mol):
    f_SA1=[];f_SA2=[];f_SA3=[];f_SA4=[];f_SA5=[];f_SA6=[];f_SA7=[];f_SA8=[];f_SA9=[];f_SA10=[];f_SA11=[];f_SA12=[];f_SA13=[];
    f_SA14=[];f_SA15=[];f_SA16=[];f_SA18=[];f_SA19=[];f_SA20=[];f_SA17=[];f_SA21=[];f_SA22=[];f_SA23=[];f_SA24=[];f_SA25=[];
    f_SA26=[];f_SA27=[];f_SA28=[];f_SA29=[];f_SA30=[];f_SA31a=[];f_SA31b=[];f_SA31c=[];f_SA37=[];f_SA38=[];f_SA39=[];
    f_SA40=[];f_SA41=[];f_SA42=[];f_SA43=[];f_SA44=[];f_SA45=[];f_SA46=[];f_SA47=[];f_SA48=[];f_SA49=[];f_SA50=[];f_SA51=[];f_SA52=[];f_SA53=[];f_SA54=[];f_SA55=[];f_SA56=[];f_SA57=[]
    lst_chemical_BB = [f_SA1,f_SA2,f_SA3,f_SA4,f_SA5,f_SA6,f_SA7,f_SA8,f_SA9,f_SA10,f_SA11,f_SA12,f_SA13,f_SA14,f_SA15,f_SA16,f_SA17,f_SA18,f_SA19,f_SA21,f_SA22,f_SA23,f_SA24,f_SA25,f_SA26,f_SA27,f_SA28,f_SA29,f_SA30,f_SA31a,f_SA31b,f_SA31c,f_SA37,f_SA38,f_SA39,f_SA40,f_SA41,
    f_SA42,f_SA43,f_SA44,f_SA45,f_SA46,f_SA47,f_SA48,f_SA49,f_SA50,f_SA51,f_SA52, f_SA53, f_SA54,f_SA55,f_SA56,f_SA57]
    if mol != 'none':
        f_SA1.append(BB.f_SA1(mol))
        f_SA2.append(BB.f_SA2(mol))
        f_SA3.append(BB.f_SA3(mol))
        f_SA4.append(BB.f_SA4(mol))
        f_SA5.append(BB.f_SA5(mol))
        f_SA6.append(BB.f_SA6(mol))
        f_SA7.append(BB.f_SA7(mol))
        f_SA8.append(BB.f_SA8(mol))
        f_SA9.append(BB.f_SA9(mol))
        f_SA10.append(BB.f_SA10(mol))
        f_SA11.append(BB.f_SA11(mol))
        f_SA12.append(BB.f_SA12(mol))
        f_SA13.append(BB.f_SA13(mol))
        f_SA14.append(BB.f_SA14(mol))
        f_SA15.append(BB.f_SA15(mol))
        f_SA16.append(BB.f_SA16(mol))
        f_SA17.append(BB.f_SA17(mol))
        f_SA18.append(BB.f_SA18(mol))
        f_SA19.append(BB.f_SA19(mol))
        f_SA20.append(BB.f_SA20(mol))
        f_SA21.append(BB.f_SA21(mol))
        f_SA22.append(BB.f_SA22(mol))
        f_SA23.append(BB.f_SA23(mol))
        f_SA24.append(BB.f_SA24(mol))
        f_SA25.append(BB.f_SA25(mol))
        f_SA26.append(BB.f_SA26(mol))
        f_SA27.append(BB.f_SA27(mol))
        f_SA28.append(BB.f_SA28(mol))
        f_SA29.append(BB.f_SA29(mol))
        f_SA30.append(BB.f_SA30(mol))
        f_SA31a.append(BB.f_SA31a(mol))
        f_SA31b.append(BB.f_SA31b(mol))
        f_SA31c.append(BB.f_SA31c(mol))
        f_SA37.append(BB.f_SA37(mol))
        f_SA38.append(BB.f_SA38(mol))
        f_SA39.append(BB.f_SA39(mol))
        f_SA40.append(BB.f_SA40(mol))
        f_SA41.append(BB.f_SA41(mol))
        f_SA42.append(BB.f_SA42(mol))
        f_SA43.append(BB.f_SA43(mol))
        f_SA44.append(BB.f_SA44(mol))
        f_SA45.append(BB.f_SA45(mol))
        f_SA46.append(BB.f_SA46(mol))
        f_SA47.append(BB.f_SA47(mol))
        f_SA48.append(BB.f_SA48(mol))
        f_SA49.append(BB.f_SA49(mol))
        f_SA50.append(BB.f_SA50(mol))
        f_SA51.append(BB.f_SA51(mol))
        f_SA52.append(BB.f_SA52(mol))
        f_SA53.append(BB.f_SA53(mol))
        f_SA54.append(BB.f_SA54(mol))
        f_SA55.append(BB.f_SA55(mol))
        f_SA56.append(BB.f_SA56(mol))
        f_SA57.append(BB.f_SA57(mol))
    else:
        for chem_group in lst_chemical_BB:
            chem_group.append('none')
        
    out = pd.DataFrame(list(zip(f_SA1,f_SA2,f_SA3,f_SA4,f_SA5,f_SA6,f_SA7,f_SA8,f_SA9,f_SA10,f_SA11,f_SA12,f_SA13,f_SA14,f_SA15,f_SA16,f_SA17,f_SA18,f_SA19,f_SA20,f_SA21,f_SA22,f_SA23,f_SA24,f_SA25,f_SA26,f_SA27,f_SA28,f_SA29,f_SA30,f_SA31a,f_SA31b,f_SA31c,f_SA37,f_SA38,f_SA39,f_SA40,f_SA41,f_SA42,f_SA43,f_SA44,f_SA45,f_SA46,f_SA47,f_SA48,f_SA49,f_SA50,f_SA51,f_SA52, f_SA53, f_SA54,f_SA55,f_SA56,f_SA57)),
    columns=['f_SA1','f_SA2','f_SA3','f_SA4','f_SA5','f_SA6','f_SA7','f_SA8','f_SA9','f_SA10',
             'f_SA11','f_SA12','f_SA13','f_SA14','f_SA15','f_SA16','f_SA18','f_SA19','f_SA20','f_SA17','f_SA21','f_SA22','f_SA23',
             'f_SA24','f_SA25','f_SA26','f_SA27','f_SA28','f_SA29','f_SA30','f_SA31a','f_SA31b','f_SA31c',
             'f_SA37','f_SA38','f_SA39','f_SA40','f_SA41','f_SA42','f_SA43','f_SA44','f_SA45','f_SA46',
             'f_SA47','f_SA48','f_SA49','f_SA50','f_SA51','f_SA52',' f_SA53',' f_SA54','f_SA55','f_SA56','f_SA57'])
    return out

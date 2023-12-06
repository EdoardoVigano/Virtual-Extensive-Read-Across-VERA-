import Function_for_descriptor as fd #MY LIBRARY
import pandas as pd
from rdkit import Chem


def descriptor(mol: Chem.rdchem.Mol):
    
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
    tetrazole=list(); thiazole=list(); thiocyan=list(); thiophene=list(); unbrch_alkane=list(); urea=list()
    
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
    ketone_dehydro =[]
    
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
                         heteroatoms_heterocycles5_list,heteroatoms_heterocycles6_list,steroid_list,hetero5_list,hetero6_list,benzCH2_list,
                         benzaldehyde_list, biphenol_list, Ar_OR_list, Ar_R_list, op_diphenolo_OR_list, Ar_COR_list, Ar_COO_R_H_list, 
                         C_3phenols_list, Ar_Cl_Br_list, Ring_3OH_3OR_list, Sulfoxide_list, CH2_Terminal_list, Alcool_1_list, triple_list, 
                         Ar_COO, Ar_OSO2, CC4, Imidothioesters, anhydrides, carbamate, aniline_term, charge2, charge1, nitro_aliphatic, 
                         ketone_aliphatic, Ar_ketone, ketone_dehydro, Ar_bicycle, Al_bicycle]
    
    if mol != 'none':    
        Al_COO.append(Chem.Fragments.fr_Al_COO(mol))
        Al_OH.append(Chem.Fragments.fr_Al_OH(mol))
        Al_OH_noTert.append(Chem.Fragments.fr_Al_OH_noTert(mol))
        ArN.append(Chem.Fragments.fr_ArN(mol))
        Ar_N.append(Chem.Fragments.fr_Ar_N(mol))
        Ar_NH.append(Chem.Fragments.fr_Ar_NH(mol))
        Ar_OH.append(Chem.Fragments.fr_Ar_OH(mol))
        COO2.append(Chem.Fragments.fr_COO2(mol))
        C_O.append(Chem.Fragments.fr_C_O(mol))
        C_O_noCOO.append(Chem.Fragments.fr_C_O_noCOO(mol))
        C_S.append(Chem.Fragments.fr_C_S(mol))
        HOCCN.append(Chem.Fragments.fr_HOCCN(mol))
        Imine.append(Chem.Fragments.fr_Imine(mol))
        NH0.append(Chem.Fragments.fr_NH0(mol))
        NH1.append(Chem.Fragments.fr_NH1(mol))
        NH2.append(Chem.Fragments.fr_NH2(mol))
        N_O.append(Chem.Fragments.fr_N_O(mol))
        Ndealkylation1.append(Chem.Fragments.fr_Ndealkylation1(mol))
        Ndealkylation2.append(Chem.Fragments.fr_Ndealkylation2(mol))
        Nhpyrrole.append(Chem.Fragments.fr_Nhpyrrole(mol))
        SH.append(Chem.Fragments.fr_SH(mol))
        aldehyde.append(Chem.Fragments.fr_aldehyde(mol))
        alkyl_carbamate.append(Chem.Fragments.fr_alkyl_carbamate(mol))
        #alkyl_halide.append(Chem.Fragments.fr_alkyl_halide(mol))
        allylic_oxid.append(Chem.Fragments.fr_allylic_oxid(mol))
        amide.append(Chem.Fragments.fr_amide(mol))
        amidine.append(Chem.Fragments.fr_amidine(mol))
        aniline.append(Chem.Fragments.fr_aniline(mol))
        aryl_methyl.append(Chem.Fragments.fr_aryl_methyl(mol))
        azide.append(Chem.Fragments.fr_azide(mol))
        azo.append(Chem.Fragments.fr_azo(mol))
        barbitur.append(Chem.Fragments.fr_barbitur(mol))
        benzene.append(Chem.Fragments.fr_benzene(mol))
        benzodiazepine.append(Chem.Fragments.fr_benzodiazepine(mol))
        bicyclic.append(Chem.Fragments.fr_bicyclic(mol))
        diazo.append(Chem.Fragments.fr_diazo(mol))
        dihydropyridine.append(Chem.Fragments.fr_dihydropyridine(mol))
        epoxide.append(Chem.Fragments.fr_epoxide(mol))
        ester.append(Chem.Fragments.fr_ester(mol))
     
        furan.append(Chem.Fragments.fr_furan(mol))
        guanido.append(Chem.Fragments.fr_guanido(mol))
        halogen.append(Chem.Fragments.fr_halogen(mol))
        hdrzine.append(Chem.Fragments.fr_hdrzine(mol))
        hdrzone.append(Chem.Fragments.fr_hdrzone(mol))
        imidazole.append(Chem.Fragments.fr_imidazole(mol))
        imide.append(Chem.Fragments.fr_imide(mol))
        isocyan.append(Chem.Fragments.fr_isocyan(mol))
        isothiocyan.append(Chem.Fragments.fr_isothiocyan(mol))
        ketone.append(Chem.Fragments.fr_ketone(mol))
        ketone_Topliss.append(Chem.Fragments.fr_ketone_Topliss(mol))
        lactam.append(Chem.Fragments.fr_lactam(mol))
        lactone.append(Chem.Fragments.fr_lactone(mol))
        methoxy.append(Chem.Fragments.fr_methoxy(mol))
        morpholine.append(Chem.Fragments.fr_morpholine(mol))
        nitrile.append(Chem.Fragments.fr_nitrile(mol))
        nitro.append(Chem.Fragments.fr_nitro(mol))
        nitro_arom.append(Chem.Fragments.fr_nitro_arom(mol))
        nitro_arom_nonortho.append(Chem.Fragments.fr_nitro_arom_nonortho(mol))
        nitroso.append(Chem.Fragments.fr_nitroso(mol))
        oxazole.append(Chem.Fragments.fr_oxazole(mol))
        oxime.append(Chem.Fragments.fr_oxime(mol))
        para_hydroxylation.append(Chem.Fragments.fr_para_hydroxylation(mol))
        phenol.append(Chem.Fragments.fr_phenol(mol))
        phenol_noOrthoHbond.append(Chem.Fragments.fr_phenol_noOrthoHbond(mol))
        phos_acid.append(Chem.Fragments.fr_phos_acid(mol))
        phos_ester.append(Chem.Fragments.fr_phos_ester(mol))
        piperdine.append(Chem.Fragments.fr_piperdine(mol))
        piperzine.append(Chem.Fragments.fr_piperzine(mol))
        priamide.append(Chem.Fragments.fr_priamide(mol))
        prisulfonamd.append(Chem.Fragments.fr_prisulfonamd(mol))
        pyridine.append(Chem.Fragments.fr_pyridine(mol))
        quatN.append(Chem.Fragments.fr_quatN(mol))
        sulfide.append(Chem.Fragments.fr_sulfide(mol))
        sulfonamd.append(Chem.Fragments.fr_sulfonamd(mol))
        sulfone.append(Chem.Fragments.fr_sulfone(mol))
        term_acetylene.append(Chem.Fragments.fr_term_acetylene(mol))
        tetrazole.append(Chem.Fragments.fr_tetrazole(mol))
        thiazole.append(Chem.Fragments.fr_thiazole(mol))
        thiocyan.append(Chem.Fragments.fr_thiocyan(mol))
        thiophene.append(Chem.Fragments.fr_thiophene(mol))
        unbrch_alkane.append(Chem.Fragments.fr_unbrch_alkane(mol))
        urea.append(Chem.Fragments.fr_urea(mol))
        
        heteroatoms_heterocycles6_list.append(fd.f_heteroatoms_heterocycles6(mol))
        heteroatoms_heterocycles5_list.append(fd.f_heteroatoms_heterocycles5(mol))
        steroid_list.append(fd.f_steroid(mol))
        hetero5_list.append(fd.f_five_ring_hetero(mol))
        hetero6_list.append(fd.f_six_ring_hetero(mol))
        benzCH2_list.append(fd.f_benzCH2(mol))
        benzaldehyde_list.append(fd.f_benzaldehyde(mol))
        biphenol_list.append(fd.f_biphenol(mol))
        
        Ar_OR_list.append(fd.f_Ar_OR(mol))
        Ar_R_list.append(fd.f_Ar_R(mol))
        op_diphenolo_OR_list.append(fd.f_op_diphenolo_OR(mol))
        Ar_COR_list.append(fd.f_Ar_COR(mol))
        Ar_COO_R_H_list.append(fd.f_Ar_COO_R_H(mol))
        C_3phenols_list.append(fd.f_C_3phenols(mol))
        Ar_Cl_Br_list.append(fd.f_Ar_Cl_Br(mol))
        Ring_3OH_3OR_list.append(fd.f_Ring_3OH_3OR(mol))
        Sulfoxide_list.append(fd.f_Sulfoxide(mol))
        CH2_Terminal_list.append(fd.f_CH2_Terminal(mol))
        Alcool_1_list.append(fd.f_Alcool_1(mol))
        
        triple_list.append(fd.f_triple(mol))
        
        ether.append(fd.f_ethere2(mol)) #sostituisce ethere di rdkit
        
        #ti eri scordato Ar_COO?
        Ar_COO.append(Chem.Fragments.fr_Ar_COO(mol))
        
        # nuovo 10/1
        Ar_OSO2.append(fd.f_Ar_OSO2(mol))
        CC4.append(fd.f_CC4(mol))
        
        #erika
        
        Imidothioesters.append(fd.f_Imidothioesters(mol))
        anhydrides.append(fd.f_anhydrides(mol))
        carbamate.append(fd.f_carbamate(mol))
        
        #13/01
        aniline_term.append(fd.f_aniline_term(mol))
        
        #25/01 rdkit alkyl_halide have to be replaced
        alkyl_halide.append(fd.f_alkyl_halide(mol))
        
        #26/01
        charge2.append(fd.f_charge2(mol))
        charge1.append(fd.f_charge2(mol))
        nitro_aliphatic.append(fd.f_nitro_aliphatic(mol))
        ketone_aliphatic.append(fd.f_ketone_aliphatic(mol))
        Ar_ketone.append(fd.f_Ar_ketone(mol))

        # 09/09
        ketone_dehydro.append(fd.f_ketone_dehydro(mol))
        
        #23 Nov 23
        Ar_bicycle.append(fd.f_Ar_bicycle(mol))
        Al_bicycle.append(fd.f_Al_bicycle(mol))
        
    else:
        for chem_group in lst_chemical_groups:
            chem_group.append('none')
            
    out = pd.DataFrame(list(zip(Al_COO, Al_OH, Al_OH_noTert, ArN, Ar_N, Ar_NH, Ar_OH,COO2, C_O, C_O_noCOO, C_S,
                          HOCCN, Imine, NH0, NH1, NH2, N_O, Ndealkylation1, Ndealkylation2, Nhpyrrole, SH,
                          aldehyde, alkyl_carbamate, alkyl_halide, allylic_oxid, amide,amidine, aniline, aryl_methyl,
                          azide, azo, barbitur, benzene, benzodiazepine, bicyclic, diazo, dihydropyridine, epoxide,
                          ester,ether, furan, guanido, halogen, hdrzine, hdrzone, imidazole, imide, isocyan, isothiocyan,
                          ketone, ketone_Topliss, lactam, lactone, methoxy, morpholine, nitrile, nitro, 
                          nitro_arom, nitro_arom_nonortho, nitroso, oxazole, oxime, para_hydroxylation, phenol, phenol_noOrthoHbond, 
                          phos_acid, phos_ester, piperdine, piperzine, priamide, prisulfonamd, pyridine, quatN, 
                          sulfide, sulfonamd, sulfone, term_acetylene, tetrazole, thiazole, thiocyan, thiophene,  unbrch_alkane, urea, 
                          heteroatoms_heterocycles5_list,  heteroatoms_heterocycles6_list, steroid_list, hetero5_list, 
                          hetero6_list, benzCH2_list, benzaldehyde_list, biphenol_list, Ar_OR_list, Ar_R_list, op_diphenolo_OR_list, Ar_COR_list, 
                          Ar_COO_R_H_list, C_3phenols_list, Ar_Cl_Br_list, Ring_3OH_3OR_list, Sulfoxide_list, CH2_Terminal_list, Alcool_1_list, 
                          triple_list,Ar_COO, Ar_OSO2, CC4, Imidothioesters, anhydrides,carbamate,aniline_term,charge2,charge1,nitro_aliphatic,
                          ketone_aliphatic,Ar_ketone, ketone_dehydro,Ar_bicycle,Al_bicycle)),
                       
                  columns=["Al_COO","Al_OH","Al_OH_noTert","ArN","Ar_N","Ar_NH","Ar_OH","COO2","C_O",
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
                           "benzaldehyde","biphenol", "Ar_OR", "Ar_R", "op_diphenolo_OR", "Ar_COR", "Ar_COO_R_H", "C_3phenols", 
                           "Ar_Cl_Br", "Ring_3OH_3OR", "Sulfoxide", "CH2_Terminal", "Alcool_1", "triple_bond", "Ar_COO", "Ar_OSO2", "CC4", 
                           "Imidothioesters","anhydrides","carbamate","aniline_term","charge+","charge-","nitro_aliphatic",'ketone_aliphatic',
                           'Ar_ketone','ketone_dehydro', 'Ar_bicycle','Al_bicycle'])
    return out

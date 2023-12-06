# SCEGLI ALERTS
###############################################################################
# import BB_calculator
# import Benigni_Bossa_index as FAI

import Fish_calculator as BB_calculator
import Fish_Alerts_index as FAI

###############################################################################

# rdkit packages
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True 
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

# Data Analisys library
import pandas as pd
import numpy as np
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error

# standard library
import os
import shutil


# library useful for output
from fpdf import FPDF
import matplotlib.pyplot as plt
from rdkit.Chem.Draw import rdMolDraw2D
from fpdf import FPDF
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ColorConverter 

# personal library and functions
from Function_Vera import*
import Function_descriptor_idx as FDT_indx # bond e atoms for images


class Output:
    def __init__(self, Dataset, df_input, Tox, new2): # new2 output del passaggio prima del codice su notebook   
        self.dataset = Dataset
        self.input = df_input 
        self.tox = Tox
        self.new2 = new2

    def defineToxNotox(self):

        self.dataset_Mol = []
        for mol_input in (self.dataset['SMILES']):
            mol2 = SingleMolecule(mol_input)
            self.dataset_Mol.append(mol2)
        for i, smile in enumerate(self.dataset_Mol):
            if i==0:
                df_groups = SingleMolecule.Mol_group_removeRedundant(smile)
                df_alerts = BB_calculator.descriptor(Chem.MolFromSmiles(smile.Smile))
            else:
                df_groups = df_groups.append(SingleMolecule.Mol_group_removeRedundant(smile))
                df_alerts = df_alerts.append(BB_calculator.descriptor(Chem.MolFromSmiles(smile.Smile)))
        df_groups = df_groups.reset_index(drop=True)
        df_alerts = df_alerts.reset_index(drop=True)
        df_groups = df_groups.assign(SMILES = self.dataset['SMILES'])
        df_alerts = df_alerts.assign(SMILES = self.dataset['SMILES'], Class = self.dataset['Experimental value'])
        df_all = pd.merge(df_groups, df_alerts, how='inner',on='SMILES')

        if len(df_all)>0 and len(set(df_all['Class']))<=2 and all(df_all['Class'].apply(lambda x: x in [0,1])):
           self.carcino = df_all.loc[df_all['Class']==1]
           self.noCarcino = df_all.loc[df_all['Class']==0]
        else:
           self.carcino = df_all.loc[df_all['Class']>-2] # modificato 13/12
           self.noCarcino = df_all.loc[df_all['Class']<-2]
    
    def statistiche(self):
        # take input object from principal class Molecule_input
        std = self.dataset.describe().loc[['std'],:].values[0]
        mean = self.dataset.describe().loc[['mean'],:].values[0]
        top = list(mean+(2*std))[0]
        bottom = list(mean-(2*std))[0]
        limit = ((-2) - bottom) / (top - bottom)# modifica 13/12
        return (std, mean, top, bottom, limit)

    def pdf_output_class(self, mol, num_mol):
        self.defineToxNotox()
        _, _, top, bottom, limit = self.statistiche()
        
        drop_list2 = ['TotGrp','Dissimilarity','GROUPS','Class','Similarity','SMILES',
                'Grp_Similarity_mean','Distance','Grp_Similarity_Emilio','Experimental value']
        
        pdf = FPDF(orientation='P', unit='mm', format='A4') #nuovo pdf
        # Molecule_input = 'input.png'

        #for u, mol in tqdm(enumerate(lst)):
        # grp dell'input 
        grp_input = mol.target_grp_noredundant # valuta la presenza della colonna smiles?
        grp_input = grp_input[grp_input != 0].dropna(axis=1)

        if any(pd.Series(grp_input.keys()).apply(lambda x: x in ['charge-', 'charge+'])):
            if 'charge-' in grp_input.keys():
                grp_input = grp_input.drop('charge-',axis=1)
            if 'charge+' in grp_input.keys():
                grp_input = grp_input.drop('charge+',axis=1)

        #simili e gruppi optimized rdkit
        b = mol.filter_simil    
        
        if len(set(b['Experimental value']))<=2 and all(b['Experimental value'].apply(lambda x: x in [0,1])):
            output_type ='Class'
        else:
            output_type ='Continuo'
            
        if output_type == 'Continuo':
            """ 6/12/2023
            def predictionMode(x,y):
                mode = []
                if x == 'no': mode = 'knn'
                elif y == 'no' and x == 'yes': mode = 'reasoning'
                elif y == 'yes' and x == 'yes': mode = 'linearModel'
                else: mode = 'GB'
                return mode

            mode_ = []
            for x, y in zip(self.new2['reasoning'], self.new2['LocalModel']):
                mode_.append(predictionMode(x, y))
            """
            self.new2['PredictionMode'] = [ i+" model" for i in self.new2['pred_type']]
            out = pd.merge(self.input, self.new2.loc[:,['SMILES', 'Best Pred', 'PredictionMode']], on='SMILES')
            # out['Err'] = abs(out['Experimental value'] - out['Best Pred'])
            
        
        # ci devono essere almeno 1 simile e 1 grp nel target 
        if len(b)>0 and len(grp_input.keys())>0:
            
            # alerts nel target rdkit
            c = mol.alerts_no_zero.drop('SMILES',axis=1)
            # orthogonal research
            SaG = mol.simili6
            b['Class'] = b['Experimental value'] # per questione di nomi di colonne
            
            if len(b)>0:
                if output_type=='Class':
                    # reasoning classificazione
                    if list(grp_input.keys()) == 0:
                    # nel caso non si siano MG (ridondante con il codice di reasoning, mi seriva per velocizzare)
                        prediction = 'No MG find in target'
                    else:
                        if len(b)<2: prediction = 'No similar substance in dataset for make a prediction'
                        else:
                            a = mol.prediction
                            if len(a[0])==0: prediction = a[1]
                            else: prediction = 'NON Active for exception'
                elif output_type=='Continuo':                                                            
                    reas_intero = mol.prediction_intero

            # page pdf                        
            pdf.add_page()
            pdf.set_font("Helvetica", size = 16)
            
            # logo
            if output_type=='Continuo':
                pdf.set_xy(10.0, 19.0)
                pdf.image('Images/PREMIER_logo.png',  link='', type='', w=74, h=20)
            else:
                pdf.set_xy(10.0, 15.0)
                pdf.image('Images/LogConcert.png',  link='', type='', w=70, h=24)
                
            # titolo
            pdf.set_xy(70.0, 10.0)
            pdf.cell(80, 10, 'VERA: Similarity and Grouping:', 1)       
            pdf.set_font("Helvetica", size = 8)
            
            # input image
            Chem.Draw.MolToFile(Chem.MolFromSmiles(mol.smile),'input.png')
            pdf.set_xy(120.0, 30.0)
            pdf.image('input.png',  link='', type='', w=60, h=50)
            os.remove('input.png')
            
            # Mol for Chem
            mol_chem = Chem.MolFromSmiles(mol.smile)
            
            # GasteigerCharges
            AllChem.ComputeGasteigerCharges(mol_chem)
            contribs = [mol_chem.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol_chem.GetNumAtoms())]
            fig = SimilarityMaps.GetSimilarityMapFromWeights(
                mol_chem, contribs, colorMap='jet', contourLines=10);
            fig.savefig('GasteigerCharges.png',bbox_inches='tight')
            plt.close(fig)
            
            # 'GasteigerCharges.png' position
            pdf.set_xy(90.0, 120.0)
            pdf.image('GasteigerCharges.png',  link='', type='', w=110, h=80)
            os.remove('GasteigerCharges.png')        
            pdf.text(10, 145, txt="GasteigerCharges:")
            
            # CrippenContribs
            contribs = rdMolDescriptors._CalcCrippenContribs(mol_chem)
            fig = SimilarityMaps.GetSimilarityMapFromWeights(
                mol_chem,[x for x,y in contribs], colorMap='jet', contourLines=10);
            fig.savefig('CrippenContribs.png', bbox_inches='tight')
            plt.close(fig)
            
            # 'CrippenContribs' position
            pdf.set_xy(90.0, 195.0)
            pdf.image('CrippenContribs.png',  link='', type='', w=110, h=80)
            os.remove('CrippenContribs.png')       
            pdf.text(10, 215, txt="CrippenContribs:")
            
            # Descriptors value
            pdf.set_xy(10.0, 80.0)
            des = Descriptors_Calculator(mol.smile) 
            for nd, d in enumerate(des.keys()):
                if nd == 0:               
                    s = f'{d}: {round(des[d][0],3)} \n'
                else:
                    s = s + f'{d}: {round(des[d][0],3)} \n'
            pdf.multi_cell(100.0, 4, txt= s) 
            # sottotili vari
            pdf_h=297
            # molecule n and input
            pdf.text(10, (pdf_h/5)-10, txt=f"Molecule Input:")
            # Grouppi  
            pdf.text(10, (pdf_h/5)-5, txt="TotGroups: " +str(len(grp_input.keys())-1))
            pdf.set_xy(10,(pdf_h/5))
            pdf.multi_cell(100, 4, txt = str(list(grp_input.keys())[0:-1]))
            
            # alert
            if len(c.keys())>0:
                pdf.set_xy(10,(pdf_h/5)+15)
                for j, string_ in enumerate(c.keys()):
                    if j == 0: alert_string = str(string_)
                    else: alert_string = alert_string +' '+str(string_)
                        
                pdf.multi_cell(100, 4, txt = f"Structural Alers: {alert_string}")
                                
            if b['Similarity'][0]==1:
                if output_type == 'Class':
                    if b['Class'][0]==0:cl='NON-Active'
                    else: cl='Active'
                    pdf.set_xy(120,(pdf_h/5)+25)
                    pdf.multi_cell(100, 4, txt = f"Experimental Value is known: {cl}")
                else:
                    pdf.set_xy(120,(pdf_h/5)+25)
                    pdf.multi_cell(100, 4, txt = f"Experimental Value is known: {b['Class'][0]} -log(mg/l)")
            
            # predizione
            pdf.set_xy(120,(pdf_h/5)+30)     
            if output_type == 'Class':          
                pdf.multi_cell(90, 4, txt = f"Prediction: {prediction}")
            else:
                gb = out.drop(out[out['SMILES'].apply(lambda x: x != reas_intero['SMILES'][0])].index)
                gb.reset_index(drop=True, inplace=True)
                txt = gb['Best Pred'][0]
                type_pred = gb['PredictionMode'][0]
                if type(txt)!=str:
                    pdf.multi_cell(90, 4, txt = f"Prediction from {gb['PredictionMode'][0]}: {round(txt,3)} -log(mg/l)")
                else:
                    pdf.multi_cell(90, 4, txt = f"Prediction: No prediction")
        
            # disegni dei sei simili
            pdf.add_page()
            name_png = []
            name_barImg = []
            for j, smile in enumerate(b['SMILES']):
                    if j < 6:
                        name_png.append('Simile'+str(j)+'.png')
                        if len(c.keys())>0:
                            mol_id = Chem.MolFromSmiles(smile)
                            ind = FAI.all_index_atoms(mol_id)
                            ind_bonds = FAI.all_index_bonds(mol_id)
                            img = Draw.MolToImage(mol_id, highlightAtoms=ind, highlightBonds=ind_bonds) # highlightBonds= ,highlightColor=ColorConverter().to_rgb("aqua"))
                            img.save(name_png[j])                        
                        else: Chem.Draw.MolToFile(Chem.MolFromSmiles(smile),name_png[j])
                        if output_type == 'Continuo':
                            plt.rcParams["figure.figsize"] = (3,5)
                            value = (b['Experimental value'][j] - bottom) / (top - bottom)
                            objects = (['Toxicity'])
                            x_pos = np.arange(len(objects))
                            performance = [value]
                            if value > top: value = top
                            elif value < bottom: value = bottom
                            try:
                                if value<=limit: 
                                    plt.bar(x_pos, performance, align='center', color='g', alpha=(1-value))
                                else: 
                                    plt.bar(x_pos, performance, align='center', color='r', alpha=value)
                            except: 
                                if value<=limit: 
                                    plt.bar(x_pos, performance, align='center', color='g')
                                else: 
                                    plt.bar(x_pos, performance, align='center', color='r')
                            
                            plt.xticks(x_pos, objects)
                            plt.ylim(0, 1)
                            plt.axhline(y=limit)
                            plt.yticks([])
                            plt.savefig(f'bar_img{j}.png', bbox_inches='tight')
                            plt.close()
                            name_barImg.append(f'bar_img{j}.png')                        
                    else: break


            # specifiche dei simili
            adj = 10
            for i in range(len(name_png)):
                if i<3:

                    #Immagine
                    x_imm = 115.0 
                    y_imm = 30 + (90*i)

                    #scritte varie
                    x_write=10
                    if i == 0: 
                        y_write = 10 + (90*i) + adj
                    else: 
                        y_write = 10 + (90*i) + adj

                    #gruppi
                    x_grp=10
                    if i == 0: 
                        y_grp = 40 + (90*i) + adj
                    else: 
                        y_grp = 40 + (90*i) + adj
                        
                    # bar infoself.tox intensity ( continuo values )
                    if output_type == 'Continuo':
                        x_barInt = 180 
                        y_barInt = 36 + (90*i)                 

                else:
                    if i == 3:
                        pdf.add_page()

                    #Immagine
                    x_imm = 115.0 
                    y_imm = 15 + (90*(i-3))

                    #scritte varie
                    x_write = 10
                    y_write = 10+(90*(i-3))

                    #gruppi
                    x_grp = 10
                    y_grp = 37.8 + (90*(i-3))
                    
                    # bar infoself.tox intensity ( continuo values )
                    if output_type == 'Continuo':
                        x_barInt = 180 
                        y_barInt = 21 + (90*(i-3))   
                
                # Imagine
                pdf.set_xy(x_imm,y_imm)
                pdf.image(name_png[i],  link='', type='', w=50, h=40)                                    
                
                # bar Intensity
                if output_type == 'Continuo':
                    pdf.set_xy(x_barInt, y_barInt)
                    pdf.image(name_barImg[i],  link='', type='', w=15, h=35)
                
                # scritte
                # titolo pag 1
                if i == 0:
                    pdf.set_font("Helvetica", size = 16) 
                    pdf.set_xy(10, 10)
                    pdf.cell(160, 10, f"Molecules with the highest group similarity\n", 1)
                
                # similarity etc
                pdf.set_xy(x_write, y_write)
                pdf.set_font("Helvetica", size = 8) 
                if output_type == 'Class':           
                    if b['Class'][i] == 1: cl = 'Active'
                    else: cl='NON-Active'
                elif output_type == 'Continuo': 
                    cl = round(b['Class'][i],3)         
                if output_type == 'Class': 
                    testo = f"{i+1}) SMILES: {b['SMILES'][i]}\nClass: {cl}\nSimilarityVega: {b['Similarity'][i]} \nSimilarityGrouping: {b['Grp_Similarity_Emilio'][i]} \nGrp_Similarity_Mean: {b['Grp_Similarity_mean'][i]} \nGroups in common: {b['GROUPS'][i]} \nTotGroups {b['TotGrp'][i]}\n"
                else: 
                    testo = f"{i+1}) SMILES: {b['SMILES'][i]}\nExpValue: {cl} -log(mg/l)\nSimilarityVega: {b['Similarity'][i]} \nSimilarityGrouping: {b['Grp_Similarity_Emilio'][i]} \nGrp_Similarity_Mean: {b['Grp_Similarity_mean'][i]} \nGroups in common: {b['GROUPS'][i]} \nTotGroups {b['TotGrp'][i]}\n"
                pdf.multi_cell(200, 4, 
                    txt = testo)
                
                # gruppi
                k = pd.DataFrame(b.iloc[i,:]).copy()
                k = k.transpose()
                k = k.drop(drop_list2, axis=1)
                k = k[k != 0].dropna(axis=1)
                for n_key, key in enumerate(k.keys()):
                    if n_key == 0:
                        s1 = self.tox.drop(self.tox[self.tox['GROUPS']!=key].index)
                        s1 = s1.reset_index()
                        s2 = s1['Tox'][0]
                        s = f"{key}   Stat: {str(s2)}\n" #Prev%: {str(s3)};\n"
                    else:
                        s1 =self.tox.drop(self.tox[self.tox['GROUPS']!=key].index)
                        s1 = s1.reset_index()
                        s2 = s1['Tox'][0]
                        s = s + f"{key}   Stat: {str(s2)}\n" #f"{key}   Stat: {str(s2)}  Prev%: {str(s3)};\n"

                #gruppi
                pdf.set_xy(x_grp, y_grp)
                pdf.multi_cell(90, 3, txt = f'GROUPS:\n{s}',border=0)

                #alerts                
                pdf.multi_cell(90, 3, txt = f'Alerts: {str(list(c))}',border=0)
            
            for j in name_png:
                os.remove(j)
            for j in name_barImg:
                os.remove(j)
                
                
            # Barchart GRP 
            # questo lo puoi fare solo se sono due classi
            
            pdf.add_page()
            # titolo
            pdf.set_font("Helvetica", size = 16) 
            pdf.set_xy(10, 10)
            pdf.cell(160, 10, "Molecule Input: Bar Chart for groups in Dataset find in target", 1)
            
            pdf.set_font("Helvetica", size = 8)
            data1 = []
            data2 = []
            grp_list = list(grp_input.keys())[:-1]

            for i, grp in enumerate(grp_list):
                data1.append(self.carcino[grp].sum())
                data2.append(self.noCarcino[grp].sum())

            X = np.arange(len(grp_list))
            fig = plt.figure(figsize=(10,5))
            ax = fig.add_axes([0,0,1,1])
            ax.bar(X + 0.00, data1[:], color = 'r', width = 0.25)
            ax.bar(X + 0.25, data2[:], color = 'b', width = 0.25)
            ax.set_ylabel('N°_groups', size = 20)
            ax.set_title('Chemicals groups in Active-NoActive', size = 30)
            ax.set_xticks(X+0.10)
            ax.set_xticklabels(grp_list, rotation=45)
            ax.tick_params(axis='both', which='major', labelsize=20)
            ax.legend(labels=['Active', 'NoActive'])

            #figura
            x = 20
            y = 160 
            plt.savefig('fig.png',dpi=300, bbox_inches='tight') 
            pdf.set_xy(x, y)
            pdf.image('fig.png',  link='', type='', w=150, h=100)
            os.remove("fig.png")
            plt.close(fig=None)
            
            #Barchart Alert
            if len(list(c.keys()))>0:
                data1 = []
                data2 = []
                grp_list = list(c.keys())
                for i, grp in enumerate(grp_list):
                    data1.append(self.carcino[grp].sum())
                    data2.append(self.noCarcino[grp].sum())
                X = np.arange(len(grp_list))
                fig = plt.figure(figsize=(10,5))
                ax = fig.add_axes([0,0,1,1])
                ax.bar(X + 0.00, data1[:], color = 'r', width = 0.25)
                ax.bar(X + 0.25, data2[:], color = 'b', width = 0.25)
                ax.set_ylabel('N°_groups', size=20)
                ax.set_title('Alert in Active/NoActive', size = 30)
                ax.set_xticks(X+0.10)
                ax.set_xticklabels(grp_list, rotation=45)
                ax.tick_params(axis='both', which='major', labelsize=20)
                ax.legend(labels=['Active', 'NoActive'])

                #figura
                x = 20
                y = 30 
                plt.savefig('fig2.png', bbox_inches='tight') 
                pdf.set_xy(x, y)
                pdf.image('fig2.png',  link='', type='', w=150, h=100)
                os.remove("fig2.png")
                plt.close(fig=None)
            
        
            #Reasonig

            if len(SaG)>0: #and len(list(c.keys()))>1: # c'è alert 

                for i in range(len(SaG)):  
                    if len(SaG[i]['SMILES'])>1:

                        pdf.add_page()
                        pdf.set_font("Arial", size = 8)
                        pdf.set_line_width(0.0)
                        pdf.rect(5.0, 5.0, 200.0,287.0)
                        pdf.rect(8.0, 8.0, 194.0,282.0)
                        
                        # titolo
                        pdf.set_font("Arial", size = 16) 
                        pdf.set_xy(10, 10)
                        pdf.cell(50, 10, "Reasoning", 1)

                        pdf.set_font("Arial", size = 8)
                        pdf.set_xy(10, 25)
                        if len(list(c.keys()))>0:
                            pdf.multi_cell(160, 4, f"Ortogonal Research based on copresence of alerts {str(list(c.keys()))} and GRP {SaG[i]['Target_Group'][0]}",
                                        1)
                        else:
                            pdf.multi_cell(160, 4, f"Ortogonal Research based on copresence of group with high prevalence {SaG[i]['Tox_reason'][0]} and GRP {SaG[i]['Target_Group'][0]}",
                                        1)

                        for z in range(len(SaG[i]['SMILES'])):

                            if z%2==0:
                                x_imm = 20 
                                y_imm = 50 + (40*(z))
                            else:
                                x_imm = 120.0 
                                y_imm = 50 + (40*(z-1)) 

                            pdf.set_xy(x_imm,y_imm)
                            
                            mol_id = Chem.MolFromSmiles(SaG[i]['SMILES'][z])
                            dic = FDT_indx.inx_groups(mol_id)
                            grp_name = SaG[i]['Target_Group'][0]
                            if type(dic[grp_name][0])!=int:
                                ind_target_grp_atoms = dic[grp_name][0][1]
                                ind_target_grp_bonds = dic[grp_name][0][0]
                            else:
                                ind_target_grp_atoms = []
                                ind_target_grp_bonds =  []
                                ind_target_grp_atoms.append(dic[grp_name][0])
                                ind_target_grp_bonds.append(dic[grp_name][0])
                                print(grp_name)
                            
                            if len(c.keys())>0:
                                
                                ind_tox_grp_atoms = FAI.all_index_atoms(mol_id)
                                ind_tox_grp_bonds = FAI.all_index_bonds(mol_id)
                                
                            else:
                                
                                grp_name = SaG[i]['Tox_reason'][0]
                                ind_tox_grp_atoms = dic[grp_name][0][1]
                                ind_tox_grp_bonds = dic[grp_name][0][0] 
                                
                            atom_cols= {}
                            for at in ind_tox_grp_atoms:
                                atom_cols[at] = (1, 0.5, 0.5)

                            bond_cols = {}
                            for bd in ind_tox_grp_bonds:
                                bond_cols[bd] = (1, 0.5, 0.5)

                            for at in ind_target_grp_atoms:
                                atom_cols[at] = (0.5, 0.6, 1)

                            for bd in ind_target_grp_bonds:
                                bond_cols[bd] = (0.5, 0.6, 1)

                            hit_ats = ind_target_grp_atoms + ind_tox_grp_atoms
                            hit_bonds = ind_target_grp_bonds + ind_tox_grp_bonds

                            drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
                            drawer.DrawMolecule(mol_id, highlightAtoms=hit_ats,
                                                    highlightAtomColors=atom_cols,
                                                    highlightBonds=hit_bonds,
                                                    highlightBondColors=bond_cols)
                            drawer.FinishDrawing()
                            with open('fig2'+str(i)+str(z)+'.png','wb') as f:
                                f.write(drawer.GetDrawingText())
                                
                            pdf.image('fig2'+str(i)+str(z)+'.png',  link='', type='', w=50, h=40)
                            pdf.set_xy(x_imm-5, y_imm-15)
                            
                            if output_type == 'Class':
                                if SaG[i]['Class'][z] == 1: cl = 'Active'
                                else: cl='NON-Active'
                            elif output_type == 'Continuo': cl = round(SaG[i]['Class'][z],3)
                            
                            if output_type == 'Class': 
                                testo = f"Smile: {SaG[i]['SMILES'][z]}\nClass: {cl}\nSimilarity: {SaG[i]['Similarity'][z]}"
                            else: 
                                testo = f"Smile: {SaG[i]['SMILES'][z]}\nExpValue: {cl} -log(mg/l)\nSimilarity: {SaG[i]['Similarity'][z]}"                         
                            
                            pdf.multi_cell(80, 4, txt = testo)
                            os.remove('fig2'+str(i)+str(z)+'.png')    
            
            
            # linear model for continuo:
            if output_type == 'Continuo': 
                pdf.add_page()
                pdf.set_font("Arial", size = 8)
                pdf.set_line_width(0.0)
                pdf.rect(5.0, 5.0, 200.0,287.0)
                pdf.rect(8.0, 8.0, 194.0,282.0)

                # titolo
                pdf.set_font("Arial", size = 16) 
                pdf.set_xy(10, 10)
                pdf.cell(50, 10, "Local Linear Model", 1)
                pdf.set_font("Arial", size = 8)
                pdf.set_xy(10, 25)
                gb = out.drop(out[out['SMILES'].apply(lambda x: x != reas_intero['SMILES'][0])].index)
                gb.reset_index(drop=True, inplace=True)
                txt = gb['Best Pred'][0]
                type_pred = gb['PredictionMode'][0]
                
                if type_pred == 'linearModel':
                    pdf.multi_cell(160, 4, "The target assessment was made by this local model: " )
                else:
                    pdf.multi_cell(160, 4, f"The target assessment was made by using {type_pred} model not by this Local model")
                
                if 'no data' in reas_intero['exp value local'][0]:
                    reas_intero.drop(reas_intero[reas_intero['exp value local']=='no data'].index, inplace=True)
                
                if len(reas_intero['exp value local'])>0:
                    x = reas_intero['exp value local'][0]
                    y = reas_intero['predictions'][0]
                    reas_intero_frame = pd.DataFrame(zip(x, y), columns=['exp value local', 'predictions'])
                    reas_intero_frame.drop(reas_intero_frame[reas_intero_frame['exp value local']=='no data'].index, inplace=True)
                    x = reas_intero_frame['exp value local'].tolist()
                    y = reas_intero_frame['predictions'].tolist()
                    max_ = max([max(x),max(y)])
                    min_ = min([min(x),min(y)])
                    
                    plot = sns.scatterplot(x=x, y=y)
                    plt.title('Linear Model')
                    plt.xlabel('Experimental Value')
                    plt.ylabel('Predicted Value')
                    plot.set_xlim(min_-0.5, max_+0.5)
                    plot.set_ylim(min_-0.5, max_+0.5)
                    fig = plot.get_figure()  
                    fig.savefig("LinearModel.png", bbox_inches='tight') 
                    plt.close()     
                    pdf.set_xy(30, 40)
                    pdf.image("LinearModel.png",  link='', type='', w=75, h=75) 
                    os.remove("LinearModel.png")            
                    mean_square_error_value_LM = round(mean_squared_error(x, y),2)
                    # r2_score_value = round(r2_score(x, y), 2)
                    pdf.set_font("helvetica", size = 12) 
                    try:
                        if type_pred == 'linearModel': 
                            testo = f"Mean square error: {mean_square_error_value_LM}\nPrediction Target: {str(txt)}"
                        else:testo = f"Mean square error: {mean_square_error_value_LM}" 
                    except: testo = f"Mean square error: {mean_square_error_value_LM}"
                    pdf.set_xy(30, 150)
                    pdf.multi_cell(80, 4, txt = testo)
                
                        
            # Scaffold
            
            scaffold = MoleculeInput.reasoning_on_scaffold(mol)
            scaffold.reset_index(drop=True, inplace=True)
            
            scaffold_ = scaffold_index(mol.smile)
            if len(scaffold)>0 and len(scaffold_[0])>0:
                
                pdf.add_page()
                pdf.set_font("Arial", size = 8)
                pdf.set_line_width(0.0)
                pdf.rect(5.0, 5.0, 200.0,287.0)
                pdf.rect(8.0, 8.0, 194.0,282.0)
                
                # titolo
                pdf.set_font("Arial", size = 16) 
                pdf.set_xy(10, 13)
                pdf.cell(40, 15, "Scaffold", 1)
                pdf.set_font("Arial", size = 8)
                pdf.set_xy(10, 30)          
                pdf.multi_cell(100, 4, txt = f"Scaffold analysis")
                
                            
                for z in range(len(scaffold)):
                    
                    scaffold_idx = scaffold_index(scaffold['SMILES'][z])
                    
                    if z==0: 
                        Chem.Draw.MolToFile(scaffold_idx[2],'scaffold.png')
                        x_imm = 70
                        y_imm = 30
                        pdf.set_xy(x_imm,y_imm)
                        pdf.image('scaffold.png',  link='', type='', w=70, h=60)
                        os.remove('scaffold.png')
                        
                    if z<4:
                        if z%2==0:
                            x_imm = 20 
                            y_imm = 50 + (40*(z+2))
                        else:
                            x_imm = 120.0 
                            y_imm = 50 + (40*(z+1)) 

                        pdf.set_xy(x_imm,y_imm)
                        
                        
                        # Chem.Draw.MolToFile(Chem.MolFromSmiles(scaffold['SMILES'][z]),'fig2'+str(z)+'.png')

                        img = Draw.MolToImage(Chem.MolFromSmiles(scaffold['SMILES'][z]), 
                                            highlightAtoms=scaffold_idx[1],
                                            highlightBonds=scaffold_idx[0],
                                            highlightColor=ColorConverter().to_rgb('aqua')
                                            ) 

                        img.save('fig2'+str(z)+'.png') 
                        pdf.image('fig2'+str(z)+'.png',  link='', type='', w=50, h=40)

                        pdf.set_xy(x_imm-5, y_imm-15)
                        if output_type == 'Class':
                            if scaffold['Class'][z] == 1:cl = 'Active'
                            else: cl='NON-Active'
                        elif output_type == 'Continuo':
                            cl = round(scaffold['Class'][z],3)


                        if output_type == 'Class': 
                            testo = f"Smile: {scaffold['SMILES'][z]}\nClass: {cl}\nSimilarity: {scaffold['Similarity'][z]}"
                        else: 
                            testo = f"Smile: {scaffold['SMILES'][z]}\nExpValue: {cl} -log(mg/l)\nSimilarity: {scaffold['Similarity'][z]}"                        

                        pdf.multi_cell(80, 4, txt = testo)
                        os.remove('fig2'+str(z)+'.png')
                            
                        
            # descriptor similarity
            
            desc_similarity = MoleculeInput.descriptor_similarity(mol)
            
            if len(desc_similarity)>0:
                
                pdf.add_page()
                pdf.set_font("Arial", size = 8)
                pdf.set_line_width(0.0)
                pdf.rect(5.0, 5.0, 200.0,287.0)
                pdf.rect(8.0, 8.0, 194.0,282.0)
                
                # titolo
                pdf.set_xy(12, 13)
                pdf.set_font("Arial", size = 16)             
                pdf.cell(40, 15, "Descriptors", 1)
                pdf.set_font("Arial", size = 8)
                pdf.set_xy(12, 30)
                pdf.multi_cell(100, 4, txt = f"Calculator Similarity descriptor")
                
                            
                for j, z in enumerate(range(len(desc_similarity))):
                    
                    if j<6:
                        if z%2==0:
                            x_imm = 20 
                            y_imm = 60 + (40*(z))
                        else:
                            x_imm = 120.0 
                            y_imm = 60 + (40*(z-1)) 

                        pdf.set_xy(x_imm,y_imm)
                        Chem.Draw.MolToFile(Chem.MolFromSmiles(desc_similarity['SMILES'][z]),'figdescr'+str(z)+'.png')

                        pdf.image('figdescr'+str(z)+'.png',  link='', type='', w=40, h=30)
                        pdf.set_xy(x_imm-5, y_imm-25)
                        if output_type == 'Class':
                            if desc_similarity['Class'][z] == 1:cl = 'Active'
                            else: cl='NON-Active'
                        elif output_type == 'Continuo':
                            cl = round(desc_similarity['Class'][z],3)
                        
                        
                        if output_type == 'Class': 
                            testo = f"Smile: {desc_similarity['SMILES'][z]}\nClass: {cl}\nDistanceCityBlock: {round(desc_similarity['descriptor_similarity'][z],3)}\nLogP:{desc_similarity['LogP'][z]}\nMW: {desc_similarity['ExactMolWt'][z]}"
                        else: 
                            testo = f"Smile: {desc_similarity['SMILES'][z]}\nExpValue: {cl} -log(mg/l)\nDistanceCityBlock: {round(desc_similarity['descriptor_similarity'][z],3)}\nLogP:{desc_similarity['LogP'][z]}\nMW: {desc_similarity['ExactMolWt'][z]}"                      
                            
                        pdf.multi_cell(80, 4, txt = testo)
                        os.remove('figdescr'+str(z)+'.png')

            pdf.output(f"Mol{num_mol}.pdf") 

    def makepdf(self, lst):
        '''make pdf: list of pdf output for every element of the obj in lst'''
        path = os.getcwd()
        path_result = os.path.join(path, "result")
        shutil.rmtree(path_result)
        os.mkdir(path_result)
        '''
        try: 
            name = os.path.join(path_result, "result/prediction.xlsx")
            os.rename('result/prediction.xlsx', name)
        except: 
            name = os.path.join(path_result, "result/prediction_continuo.xlsx")
            os.rename('result/prediction_continuo.xlsx', name)
        '''
        for i in range(len(lst)):
            self.pdf_output_class(lst[i],i)
            mol_dis = os.path.join(path_result, f"Mol{i}.pdf")
            try: 
                os.rename(f"Mol{i}.pdf", mol_dis)
                print(f"Mol{i}: has done!")
            except: 
                print(f"Mol{i}: hasn't similar compound!")
            
                                                    
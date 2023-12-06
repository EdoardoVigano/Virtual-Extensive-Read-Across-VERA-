# import BB_calculator
import Fish_calculator as BB_calculator


import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors 
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.Scaffolds import MurckoScaffold as MS

rdDepictor.SetPreferCoordGen(True)
IPythonConsole.drawOptions.minFontSize=20

import descriptor_calculator as dc
import Function_for_descriptor_tuple as Fdt

from sklearn.preprocessing import scale
from scipy.spatial import distance

# internal pack
import os

# statistic packages
from scipy.spatial import distance
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer


global list_out 
global Dataset_class
Dataset_class = [1,0]

# lista di gruppi funzionali "poco interessanti" quindi eliminati post
list_out = ["COO2","C_O","C_O_noCOO","Ar_N","Ar_NH","Al_OH_noTert",
            "Ar_N","benzaldehyde", 'halogen','piperdine','Alcool_1','nitro', 
             'nitro_arom_nonortho'] #'hetero5', 'hetero6' - rimessi nell'output il 24 Nov 23
    


#  classe per gli oggetti utilizzati in tutto il codice
# funzioni della classe utili per i calcoli dei MG e ALerts

class SingleMolecule:
   
    def __init__(self, SMILES): 
        
        self.Smile = SMILES
                
    def Mol_group(self):
        a = pd.DataFrame(dc.descriptor(Chem.MolFromSmiles(self.Smile)))
        a = a.assign(SMILES=self.Smile)
        return a
    
    def Mol_group_no_zero_removeRedundant(self):
        a = Chem.MolFromSmiles(self.Smile)
        k = pd.DataFrame(dc.descriptor(a))
        k1 = remove_Redundant(a, k)
        k1 = k1[k1 != 0].dropna(axis=1)     
        k1['SMILES']=self.Smile        
        return k1

    def Mol_group_removeRedundant(self):
        a = Chem.MolFromSmiles(self.Smile)
        b = pd.DataFrame(dc.descriptor(a))
        b = b.assign(SMILES=self.Smile)
        return remove_Redundant(a, b)

    def Mol_group_removeRedundant_no_zero(self):
        k = self.Mol_group_removeRedundant()
        k = k[k != 0].dropna(axis=1)
        k = k.assign(SMILES=self.Smile)
        return k

    def Alerts(self):
        a = pd.DataFrame(BB_calculator.descriptor(Chem.MolFromSmiles(self.Smile)))
        a = a.assign(SMILES=self.Smile)
        return a
    
    def Alerts_no_zero(self):
        a = SingleMolecule(self.Smile)
        k = SingleMolecule.Alerts(a)
        k = k[k != 0].dropna(axis=1)
        k = k.assign(SMILES=self.Smile)
        return k



# classe di lavoro principale
class MoleculeInput:
    
    # definisci il path per istsimilarity
    path_program = os.getcwd()
    path_program = f"{path_program}/CustomFeaturesCLI/CustomFeaturesCLI.jar" 
    # path_program = "/home/evigano/Scrivania/CustomFeaturesCLI/CustomFeaturesCLI.jar" #work
    # path_program = f"/Users/edoar/Desktop/CustomFeaturesCLI/CustomFeaturesCLI.jar" #home
    # path_program = "/Users/erika/Desktop/CustomFeaturesCLI/CustomFeaturesCLI.jar"

    def __init__(self, smile, simili):

        '''simili, must be a pandas dataframe, the output from Generator_Molecule'''

        self.smile = smile
        self.simili = simili # smiles e vega simililarity dei composti con similarità maggiore di 0.650 al target

        # utilizzo la classe SingleMolecule per calcolare MG e alerts
        self.target_grp_noredundant = SingleMolecule.Mol_group_removeRedundant(SingleMolecule(self.smile))
        self.alerts = SingleMolecule.Alerts(SingleMolecule(self.smile))
        self.alerts_no_zero = SingleMolecule.Alerts_no_zero(SingleMolecule(self.smile))

        # istanze utili successivamente:
        # grp.simili: MG dei simili con 0.650 di threashold
        self.grp_simili = pd.DataFrame()
        self.filter_simil = pd.DataFrame()
        self.simili6 = []
        self.prediction = []
        self.prediction_intero = pd.DataFrame()
        self.scaffold = pd.DataFrame()
        # self.scaffold_index = []
        self.simil_descriptor_calculated = pd.DataFrame()


    # calcolo dei gruppi molecolari 
    def Calculate_grp(self):

        for j, smile in enumerate(self.simili['SMILES']):

            a = SingleMolecule(smile)
            grp = SingleMolecule.Mol_group(a)
            al_simili = SingleMolecule.Alerts(a)
            grp_no_red = SingleMolecule.Mol_group_removeRedundant(a)

            if j == 0:
                out = pd.merge(grp_no_red.copy(), al_simili, on='SMILES')
            else:
                out1 = pd.merge(grp_no_red, al_simili, on='SMILES')
                out = out.append(out1)
            
            out.reset_index(drop=True, inplace= True)

        data = pd.merge(out, self.simili, on= 'SMILES')
        self.grp_simili = data.copy()
        return data
    
    # filtro tramite presenza degli stessi alert (genera subset per ogni molecola target)
    # i simili da considerare devono avere gli stessi alerts del target

    def filter_simil(self):
        
        # alert simili

        # define target molecule 
        a = SingleMolecule(self.smile)

        # define alerts
        al_target = self.alerts.copy().drop('SMILES', axis=1)
        key_alerts = list(al_target.keys())
        al_target_noZero = self.alerts_no_zero.copy().drop('SMILES', axis=1)
        key_alerts_noZero = list(al_target_noZero.keys())

        # define MG presenti nel target
        MG_target = self.target_grp_noredundant.copy()
        key_MG = list(MG_target.drop('SMILES', axis=1).keys())
    
        # alert simili != alert target: eliminati
        
        idx = []
        for i in range(len(self.grp_simili)):
            r = al_target.loc[0,key_alerts] == self.grp_simili.loc[i,key_alerts]
            s = r.all()
            if s == False:
                    idx.append(i) # indici dei simili da eliminare che non hanno gli stessi alert del target 
        
        data = self.grp_simili.drop(idx).reset_index(drop=True)

        # conteggio dei gruppi in comune tra target e ogni simile del subset
        # conteggio dei gruppi NON in comune tra target e ogni simile del subset
        #  calcolo della similarità (Emilio_distance) di grouping (Emilio)
        grp_tot = []
        dissimilarity = []
        tot = [] 
        Emilio_distance = []

        # non considerare il conteggio della carica come MG
        grp_targhet = MG_target.loc[0,key_MG].sum() - MG_target['charge+'][0] - MG_target['charge-'][0]

        for i in range(len(data)):
            tot.append(data.loc[i,key_MG].sum() - data['charge+'][i] - data['charge-'][i])         
            group = []
            diss = []
            for key in key_MG:
                if key != 'charge+' and key !='charge-': 
                    if MG_target[key][0]>0:
                        if MG_target[key][0] == data[key][i]:group.append(data[key][i])
                        else:
                            if data[key][i]>0:
                                if data[key][i]>MG_target[key][0]:
                                    group.append(MG_target[key][0])
                                    diss.append(data[key][i] - MG_target[key][0])
                                else:
                                    group.append(data[key][i])
                                    diss.append(MG_target[key][0] - data[key][i])

                    else:
                        if data[key][i] > 0:
                            diss.append(data[key][i])
                            
                           
            dissimilarity.append(sum(diss))     
            grp_tot.append(sum(group))

        # Gouping similarity
        for i in range(len(grp_tot)):
            pt1 = (1 + (grp_tot[i]/grp_targhet))/2
            pt2 = (dissimilarity[i]/tot[i])/8
            edist = round((pt1-pt2),3)
            Emilio_distance.append(edist)
          
        c = pd.DataFrame({'GROUPS':grp_tot, 'Dissimilarity': dissimilarity, 'TotGrp':tot})
        result = pd.concat([data, c], axis=1, join="inner")
        
        result = result.reset_index(drop=True)
        result.drop(key_alerts, axis=1, inplace=True) # adesso alerts non servono più dopo il drop

        # questo non serve per te Alessio, sono implementazioni nuove:
        ##################################################################################
        desc_s = Descriptors_Calculator(self.smile)
        for smile in result['SMILES']:
            desc_s = desc_s.append(Descriptors_Calculator(smile))
        desc_s.reset_index(drop=True, inplace = True)
        
        # scale data
        x = desc_s.values
        try: 
            X = scale(x) 
            desc_dist = []
            n1 = X[0,:]
            for arr in range(len(X)):
                if arr>0:
                    n2 = X[arr,:]
                    dist = distance.cityblock(n1,n2)
                    desc_dist.append(dist)                        
            result['Distance'] = desc_dist
        except:
            result['Distance'] = 0
        ###################################################################################    
        # Emilio distance
        
        result['Grp_Similarity_Emilio'] = Emilio_distance 
        result['Grp_Similarity_mean'] = round((result['Grp_Similarity_Emilio'] 
                                               + result['Similarity'])/2, 3)#non definitivo
        
        
        #forse è meglio chiedere all'utente come ordinare
        result = result.sort_values(['Grp_Similarity_Emilio','GROUPS','Similarity','Grp_Similarity_mean'], 
                                    ascending=False) 
        result = result.reset_index(drop=True)
        
        # eliminazione tra i simili delle forme cariche di molecole già presenti in forma neutra: 
        """
        if any(result['charge+'][0:20])>0 and any(result['charge-'][0:20])>0:
            
            charge = result.drop(result[result['charge+']+result['charge-']==0].index)
            result2 = result.drop(result[result['charge+']+result['charge-']!=0].index)
            smile_remove = []
            for n_charge_mol, charge_mol in enumerate(charge['SMILES']):
                if n_charge_mol<7:
                    for n_result2_mol, result2_mol in enumerate(result2['SMILES']):
                        if n_result2_mol<7:
                            command = ["java", "-jar", "-Xmx8000m", MoleculeInput.path_program, "-similarity", charge_mol, result2_mol]
                            subprocess.call(command)
                            similarity_charge  =  pd.read_csv("similarity.csv", delimiter= '\t', encoding='cp1252')
                            os.remove('similarity.csv')
                            if similarity_charge['Similarity'][0]>0.99:
                                smile_remove.append(similarity_charge['SMILES 1'][0])
                    
            result = result.drop(result[result['SMILES'].apply(lambda x: x in smile_remove)].index)
            result = result.reset_index(drop=True)
        """    
        result = result.drop(['charge+','charge-'], axis=1)
        if 'level_0' in result.keys(): del result['level_0']
        
        
        ##aggiunta 21/03, per evitare NA
        if all(result['Grp_Similarity_Emilio'].isnull()): result['Grp_Similarity_mean']=result['Similarity']

        # definisci l'istanza subset filtrato
        self.filter_simil = result.copy()
        
        return result
        
    # function per cercare simili che condividono con il target la presenza di un preciso MG
    def ortogonal_research(self, Tox):

        mol = SingleMolecule(self.smile)
        target = SingleMolecule.Mol_group_no_zero_removeRedundant(mol).drop('SMILES',axis=1)
        alert = SingleMolecule.Alerts_no_zero(mol).drop('SMILES',axis=1)


        df_simil_rdkit_noRedundant = self.filter_simil.copy()
        df_simil_rdkit_noRedundant.rename(columns={'Experimental value':'Class'}, inplace=True)

        # rimuovi il grp "cariche", non interessa ai fini del read across
        if any(pd.Series(target.keys()).apply(lambda x: x in ['charge-', 'charge+'])):
            if 'charge-' in target.keys():
                target = target.drop('charge-',axis=1)
            if 'charge+' in target.keys():
                target = target.drop('charge+',axis=1)
                
        simili6 = []
        gruppi_target = target.copy()
        grp = list(gruppi_target.keys())
        
        if len(grp)>0: # il target deve avere almeno un gruppo funzionale
            
            # se c'è almeno un alerts
            if len(list(alert.keys()))>0:

                # genera n- subset dei simili che possiedono alert e gruppo funzionale i-esimo
                # un subset per ogni gruppo funzionale trovato nel target:
                
                # alert + grp1: subset 1
                # alert + grp2: subset 2
                
                for key in target.keys():
                    k = df_simil_rdkit_noRedundant.copy()
                    k = k.drop(k[k[key]==0].index).sort_values(['Similarity'], ascending=False)
                    k = k.reset_index(drop=True)
                    k['Target_Group'] = key
                    if len(k)>1:
                        # 29/09
                        if k['Similarity'][0]==1: out = k[1:7]
                        else: out = k[0:6]
                        # out = k[0:6]
                        out = out.reset_index(drop=True)
                        simili6.append(out)

            else:

                # se alert non presente nel target devi cercare il motivo di tox più rilevante quindi
                # ti basi sul gruppo a prevalenza tox maggiore nel database
                
                # dalla tabella delle prevalenze (funzione prevalence, riga 1133), quali grp sono presenti nel target?
                Tox2 = Tox.drop(Tox[Tox['GROUPS'].apply(lambda x: x not in grp)].index).reset_index(drop=True)

                # qui sotto ci sono un pò di condizioni di comportamento dell'algoritmo,
                # ma dipendono dalla natura dell'endpoint se continuo o meno
                
                # classificatore
                if len(set(df_simil_rdkit_noRedundant['Class']))==2 and all(df_simil_rdkit_noRedundant['Class'].apply(lambda x: x in [0,1])):
                    alert_gruppi = Tox2.drop(Tox2[Tox2['Prevalence%']!=max(Tox2['Prevalence%'])].index).reset_index(drop=True)
                
                elif len(set(df_simil_rdkit_noRedundant['Class']))==1 and all(df_simil_rdkit_noRedundant['Class'].apply(lambda x: x in [0,1])):
                    alert_gruppi = Tox2.drop(Tox2[Tox2['Prevalence%']!=max(Tox2['Prevalence%'])].index).reset_index(drop=True)
                
                # se fish 
                else: alert_gruppi = Tox2.drop(Tox2[Tox2['Prevalence%']!=max(Tox2['Prevalence%'])].index).reset_index(drop=True)
                '''    
                # se continuo l'endpoint (questa roba in teoria solo se non ho treshold)
                else:
                    alert_gruppi = Tox2.drop(Tox2[Tox2['Coeff']!=max(Tox2['Coeff'])].index).reset_index(drop=True)
                '''
                
                # filtra i simili analogamente alle fasi precedenti ma considera come alert il MG a prevalenza maggiore,
                # --> concettualmente se non hai un alert di riferimento ti basi sulla stat dei gruppi
             
                k1 = df_simil_rdkit_noRedundant.drop(
                    df_simil_rdkit_noRedundant[df_simil_rdkit_noRedundant[alert_gruppi['GROUPS'][0]]==0].index)
                gruppi_target_combine = gruppi_target.drop([alert_gruppi['GROUPS'][0]],axis=1)

                for key in gruppi_target_combine.keys():
                    if key not in ['charge+','charge-']: 
                        k = k1.copy()
                        k = k.drop(k[k[key]==0].index).sort_values(['Similarity'], ascending=False)
                        k = k.reset_index(drop=True)

                    if len(k)>0:

                        k['Target_Group'] = key
                        k['Tox_reason'] = alert_gruppi['GROUPS'][0]
                        # modifica 29/09
                        if k['Similarity'][0]==1: out = k[1:7]
                        else: out = k[0:6]
                        # out = k[0:6]
                        out = out.reset_index(drop=True)
                        simili6.append(out)
        self.simili6 = simili6                 
        return simili6

    # procedura di reasoning per la classificazione
    def reasoning(self):

        classification = 'not define'
        exception = []
        out1 = []
        # simili a 0.650 filtrati dagli alert
        df_simil_rdkit_noRedundant = self.filter_simil.copy()
        df_simil_rdkit_noRedundant.rename(columns={'Experimental value':'Class'},inplace=True)
        # simili dei grouping
        simili6 = self.simili6
        alert = self.alerts_no_zero.copy().drop('SMILES', axis=1)  
        # se non ci sono alemno 2 composti simili non si fa l'assessment
        if len(df_simil_rdkit_noRedundant)<3:  classification = 'No data'

        # verifica per quale similarity è sorted df_simil_rdkit_noRedundant
        # verifica per quale similarity è sorted df_simil_rdkit_noRedundant
        # df_simil_rdkit_noRedundant.sort_values(['Similarity'], ascending = False, inplace=True)  test
        
        # rimuovi l'identità se presente.. questo serviva a noi per i test, in teoria se trova l'identità andrebbe classificato e stop.
        

        # 22/09 PROVA A TOGLIERRE STA COSA PER BOTANICALS
        '''
        if df_simil_rdkit_noRedundant['Similarity'][0]==1: 
            df_simil_rdkit_noRedundant = df_simil_rdkit_noRedundant[1:]
        '''

        # Fase 1: 
        # Fase 1:
        
        # considera i primi dieci più simili secondo la similarità di grouping
        
        # se i simili sono praticamente tutti di un valore sperimentale allora classifica in accordo

        if all(df_simil_rdkit_noRedundant['Class'][:10]==Dataset_class[0]) and len(df_simil_rdkit_noRedundant)>9:
            # classification = 'Active***(1)'
            classification = 'Active***: the 10 more similar compounds are all active.'

        elif all(df_simil_rdkit_noRedundant['Class'][:10]==Dataset_class[1]) and len(df_simil_rdkit_noRedundant)>9: 
            # classification = 'NON-Active***(1)'
            classification = 'NON-Active***: the 10 more similar compounds are all Inactive.'
        
        # se su i 10 più simili si può già dire qualcosa (con una certa sicurezza) a livello statistico allora 
        # la procedura di reasoning sui gruppi non è necessaria

        elif len(list(alert.keys()))==0: #23/02 10 similar only for target without alert

            simili10 = df_simil_rdkit_noRedundant[:10]
            if len(simili10)==10:
                simili10_canc = simili10.drop(simili10[simili10['Class']==Dataset_class[1]].index)
                simili10_nocanc = simili10.drop(simili10[simili10['Class']==Dataset_class[0]].index)
                carci_value = round(sum(simili10_canc['Grp_Similarity_mean'])/len(simili10_canc),3)
                nocarci_value = round(sum(simili10_nocanc['Grp_Similarity_mean'])/len(simili10_nocanc),3)
                if len(simili10_canc)>len(simili10_nocanc) and carci_value>nocarci_value:
                    if len(simili10_canc) > 8 : 
                        classification = 'Active***: 9 of 10 most similar compounds are active'
                        #  classification = 'Active***(2)'
                    elif len(simili10_canc)<=8 and len(simili10_canc)>=7: 
                        # classification =  'Active**(2)'
                        classification = 'Active**: 8 of 10 most similar compounds are active'
                    elif len(simili10_canc)>6: 
                        # classification = 'Active*(2)'
                        classification =  'Active**: 7 of 10 most similar compounds are active and the average similarity of the active compounds is greater than the inactive'
                elif len(simili10_canc)>len(simili10_nocanc) and carci_value<nocarci_value: classification = 'equivocal2.0'
                else:
                    if len(simili10_nocanc) > 8 : 
                        # classification = 'NON-Active**(2)'
                        classification = 'NON-Active*** at least 9 of 10 most similar compounds are Inactive'
                    elif len(simili10_nocanc)<=8 and len(simili10_nocanc)>7: 
                        classification =  'NON-Active** 8 of 10 most similar compounds are Inactive'
                    elif len(simili10_nocanc)>=5 and len(simili10_nocanc)<=7: classification = 'equivocal2.0'
            else:
                classification = 'NO 10 data'

            
        # se la classificazione non è esaustiva:
        
        # seconda fase:
        # seconda fase:
        
        # se non è stato classificato in precedenza basandosi sui 10 più simili
        if classification == 'equivocal2.0' or classification == 'not define' or classification == 'NO 10 data':
            
            # se tutti i più simili Grouping emilio sono active ( magari un pò ridondante con la prima.. ma è una casistica rara)

            if all(df_simil_rdkit_noRedundant['Class'][:6]==Dataset_class[0]) and len(list(alert.keys()))>0:
                classification = 'Active***: 6 more similar compounds are all Active.'
                # classification = 'Active***(3)'
            elif all(df_simil_rdkit_noRedundant['Class'][:6]==Dataset_class[1]): 
                # classification = 'NON-Active***(3)'
                classification = 'NON-Active***: 6 more similar compounds are all Inactive.'

            # se i 6 più simili Grouping emilio non hanno tutti lo stesso vaore sperimentale 
            # allora procedura di reasoning:
            
            else: 
                if len(simili6) != 0:
                    
                    carcinogen1 = []
                    Nocarcinogen1 = []                   

                    # define active - no active cluster
                    # ogni cluster generato dall'ortogonal research va classificato in toto
                    # come cluster Active o NON active: serve per definire la caratteristica che contraddistingue
                    # il cluster come possibile caratteristica che influenza la scelta per la classificazione
                    
                    for grp_n, grp in enumerate(simili6):
                        # se il cluster ha almeno due elementi viene considerato
                        if len(grp)>1:

                            # quanti elementi nel cluster hanno valore sperimentale attivo o non attivo?
                            carcinogen1.append(len(grp.loc[grp['Class']==Dataset_class[0]]))
                            Nocarcinogen1.append(len(grp.loc[grp['Class']==Dataset_class[1]]))

                            # se in percentuale più del 50% è attivo allora clusterActive
                            if round((len(grp.loc[grp['Class']==Dataset_class[0]])/(len(grp)))*100,0)>50:
                                simili6[grp_n]['Active/Noactive'] = 'active'

                            # se in percentuale più del 65% è Inattivo allora clusterNoactive 
                            elif round((len(grp.loc[grp['Class']==Dataset_class[1]])/(len(grp)))*100,0)>65:
                                simili6[grp_n]['Active/Noactive'] = 'Noactive'

                            # se non rientra in nessuna casistica non lo si considera influente    
                            else:
                                simili6[grp_n]['Active/Noactive'] = 'No data'
                        else:
                            carcinogen1.append(-1)
                            Nocarcinogen1.append(-1)
                            simili6[grp_n]['Active/Noactive'] = 'No data'

                    
                    # ricerca regole di eccezione (vedi deliverable)
                    # se un cluster contiene solo elementi Non-active sperimentalmente
                    # allora la caratteristica (MG) considerata da questo cluster potrebbe annullare l'effetto dell'alert (o MG alta prevalenza)

                    if 0 in carcinogen1 and len(list(alert.keys()))>0:
                        idx = [ i for i in range(len(carcinogen1)) if carcinogen1[i] == 0 ]
                        exception = []
                        for i in idx:
                            exception.append(simili6[i]['Target_Group'][0])

                    # se active>50% per ogni grp (questo sarà rindondante con condizioni sotto)
                    number = sum(d1 > 2 for d1 in carcinogen1)
                    if number == len(simili6) and number>0: # Modifica: number>1 aggiunto il 9/09
                        # classification = 'Active***(4)'
                        # modifica reliability 5/settembre 23
                        if number==1: classification = 'Active*: Each cluster has half of the compounds labelled as Active.'
                        elif number==2: classification = 'Active**: Each cluster has half of the compounds labelled as Active.'
                        else: classification = 'Active***: Each cluster has half of the compounds labelled as Active.'
                    else: 
                        
                        # merge dei cluster portandosi dietro la relativa etichetta (Active, NON-Active),
                        for j1, i1 in enumerate(simili6): 
                            if j1 == 0: 
                                db = i1.copy()
                            else:
                                db = db.append(i1)

                                
                        # si generano due CLUSTER che dovranno essere confrontati:
                        # uno contiene tutti i subset definiti attivi e l'altro tutti quelli inattivi:
                        
                        # valutazione di quanti sono active o inactive e quanti sono in comune ai due CLUSTER:      
                              
                        db = db.reset_index(drop=True)

                        # rimuovi i cluster non definiti
                        db = db.drop(db[db['Active/Noactive'] == 'No data'].index).reset_index(drop=True)

                        # cerca quelli che nel CLUSTER active hanno dato sperimentale active
                        # cerca quelli che nel CLUSTER Inactive hanno dato sperimentale Inactive

                        db_active = db.drop(db[db['Active/Noactive'] == 'Noactive'].index)
                        db_Noactive = db.drop(db[db['Active/Noactive'] == 'active'].index)

                        db_active_smile = db_active[['SMILES','Active/Noactive','Class']]
                        db_Noactive_smile = db_Noactive[['SMILES','Active/Noactive','Class']]

                        # elimina i duplicati nei due CLUSTER
                        db_active_smile = db_active_smile.drop_duplicates().reset_index(drop=True)
                        db_Noactive_smile = db_Noactive_smile.drop_duplicates().reset_index(drop=True)

                        # trova elementi in comune ai due CLUSTER
                        db_merg = pd.merge(db_active_smile, db_Noactive_smile, how='inner', on='SMILES')

                        # considerazioni per classificazione:

                        # se sono presenti solo cluster active:
                        if len(db_Noactive_smile)==0 and len(db_active_smile)>1: 
                            # se siamo nel caso di equivocal:
                            if classification == 'equivocal2.0': 
                                # classification = 'Active*(5)'
                                classification = 'Active*: All cluster has defined as Active'
                            # altrimenti
                            else: 
                                # classification = 'Active***(5)'
                                classification = 'Active*: All cluster has defined as Active'

                        # se è presente solo 1 cluster ed è Active: 
                        elif len(db_Noactive_smile)==0 and len(db_active_smile)==1:
                            
                            # considera anche i più simili
                            simili10 = df_simil_rdkit_noRedundant[:10]
                            simili10_canc = simili10.drop(simili10[simili10['Class']==Dataset_class[1]].index)
                            
                            # tra i simili se predominano quelli experimental tox allora active
                            if len(simili10_canc)>(len(simili10)-len(simili10_canc)):
                                # classification = 'Active**(5)'
                                classification = 'Active**: Is present only one Active cluster and the 10 more similar compounds are predominantly toxic'
                            else:
                                classification = 'Equivocal 3.0'

                        # se sono presenti solo clusters Inactive:            
                        if len(db_active_smile)==0 and len(db_Noactive_smile)>1: 
                            if classification == 'equivocal2.0': 
                                # classification = 'NON-Active***(6)'
                                classification = 'NON-Active*: All cluster has defined as Inactive'
                            else: 
                                # classification = 'NON-Active*(6)'
                                if len(db_Noactive_smile)>2: 
                                    classification = 'NON-Active***: All cluster has defined as Inactive'
                                elif len(db_Noactive_smile)>1: 
                                    classification = 'NON-Active**: All cluster has defined as Inactive'
                                else:
                                    classification = 'NON-Active*: All cluster has defined as Inactive'
                                

                        # se è presente solo 1 cluster ed è Inctive:        
                        elif len(db_Noactive_smile)==1 and len(db_active_smile) == 0:
                            simili10 = df_simil_rdkit_noRedundant[:10]
                            simili10_nocanc = simili10.drop(simili10[simili10['Class']==Dataset_class[0]].index)
                            # tra i simili se predominano quelli experimental notox allora non-active
                            if len(simili10_nocanc)>(len(simili10)-len(simili10_nocanc)):
                                # classification = 'NON-Active*(7)'
                                classification = 'NON-Active**: Is present only one Inactive cluster and the 10 more similar compounds are predominantly NON-toxic'
                           
                            else:
                                classification = 'Equivocal 3.0'  


                        # modifica: va considerato il rateo generale?
                        
                        # questa considerazione non è definitiva al momento
                        # questa considerazione non è definitiva al momento
                        # però è parte del codice usato ad oggi
                        
                        # considerazioni un pò più generali sul rateo di composti dei CLUSTERS active e inactive
                       
                        if len(db_active_smile)+len(db_Noactive_smile) > 0:       
                            if len(db_active_smile)>0 and len(db_Noactive_smile) > 0:
                                rateo = round((len(db_active_smile)/(len(db_active_smile)+len(db_Noactive_smile)))*100,0)
                            elif len(db_active_smile)>0 and len(db_Noactive_smile) == 0: rateo = 100
                            elif len(db_active_smile) == 0 and len(db_Noactive_smile) > 0: rateo = 1
                        else: rateo = 0

                        # se non è molto sbilanciato il numero di elementi nei due CLUSTERS 
                        # active e inactive, si valutano gli elementi in comune                        
                        if rateo<75 and rateo>25:

                            # se ci sono elementi in comune ai due CLUSTERS e non abbiamo ancora 
                            # definito l'output secondo le casistiche precedenti
                           
                            if len(db_merg) !=0 and classification == 'not define':

                                # Ci si focalizza su gli elementi in comune ai due
                                # ClUSTERS per evidenziare le caratteristiche comuni e definire quali sono 
                                # predominanti secondo il dato sperimentale:
                                # esempio: CLUSTER ACTIVE CONTIENE 3 SOSTANZE IN COMUNE AL CLUSTER INACTIVE
                                # ALLORA QUESTI COMPOSTI CHE DATO SPERIMENTALE HANNO? SONO DAVVERO ATTIVI  O INATTIVI?
                                # QUESTE CARATTERISTICHE COMUNI AI CLUSTERS ANNULLANO EFFETTO TOX O NO?
            
                                # definisco treshold che riguardano gli elementi in comune
                                # per cercare di definire quale dei due cluster predomina.
                    
                                treshold_carc= round((sum(db_merg['Class_x']==Dataset_class[0])/len(db_merg))*100,0)
                                treshold_Nocarc = round((sum(db_merg['Class_x']==Dataset_class[1])/len(db_merg))*100,0)
                                merg_treshold = round(len(db_merg)/((len(db_active_smile)+len(db_Noactive_smile)))*100,0)


                                if  merg_treshold<=16:

                                    if treshold_carc>=83:
                                        if round((sum(db_active_smile['Class']==Dataset_class[0])/len(db_active_smile))*100,0)>=83:
                                            # classification = 'Active***(8)' 
                                            classification = 'Active***: The compounds having the characteristics of several clusters are mainly toxic'
                                        elif round((sum(db_active_smile['Class']==Dataset_class[0])/len(db_active_smile))*100,0)<83 and round((sum(db_active_smile['Class']==Dataset_class[0])/len(db_active_smile))*100,0)>=66:
                                            # classification = 'Active**(8)'
                                            classification = 'Active**: The compounds having the characteristics of several clusters are mainly toxic'
                                        else: 
                                            # classification = 'Active*(8)'
                                            classification = 'Active**: The compounds having the characteristics of several clusters are mainly toxic'
                                    elif treshold_Nocarc>=83:
                                        if round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)>=83:
                                            # classification = 'NON-Active***(8)'
                                            classification = 'NON-Active***: The compounds having the characteristics of several clusters are mainly NO-toxic'
                                        elif round((sum(db_Noactive_smile['Class'] ==Dataset_class[1])/len(db_Noactive_smile))*100,0)<83 and round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)>=66:
                                            # classification = 'NON-Active**(8)'
                                            classification = 'NON-Active**: The compounds having the characteristics of several clusters are mainly NO-toxic'
                                        else: 
                                            classification = 'NON-Active*: The compounds having the characteristics of several clusters are mainly NO-toxic'
                                            # classification = 'NON-Active*(8)'

                                elif  merg_treshold<=50 and merg_treshold>16:
                                    if treshold_carc>=83:
                                        # classification = 'Active***(9)' 
                                        classification = 'Active***: The compounds having the characteristics of several clusters are mainly toxic'
                                    if treshold_carc>=50 and treshold_carc<83:
                                        # classification = 'Active*(9)'
                                        classification = 'Active*: The compounds having the characteristics of several clusters are mainly toxic'
                                    if treshold_carc>=33 and treshold_carc<50:
                                        # classification = 'Active*(9)E'
                                        classification = 'Active*: The compounds having the characteristics of several clusters are mainly toxic'
                                    if treshold_Nocarc>=83: 
                                        if round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)>=83:  
                                            # classification = 'NON-Active***(9)'
                                             classification = 'NON-Active***: The compounds having the characteristics of several clusters are mainly NO-toxic'
                                        if  round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)>=66 and round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)<83:
                                            # classification = 'NON-Active**(9)'
                                             classification = 'NON-Active**: The compounds having the characteristics of several clusters are mainly NO-toxic'
                                        if  round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)<66:   
                                            # classification = 'NON-Active*(9)'
                                             classification = 'NON-Active*: The compounds having the characteristics of several clusters are mainly NO-toxic'


                                elif  merg_treshold>50:
                                    if  treshold_Nocarc>=83: 
                                        if round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)>83:  
                                            # classification = 'NON-Active***(10)'
                                             classification = 'NON-Active***: The compounds having the characteristics of several clusters are mainly NO-toxic'
                                        if  round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)>=66 and round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)<83:
                                            # classification = 'NON-Active**(10)'
                                             classification = 'NON-Active**: The compounds having the characteristics of several clusters are mainly NO-toxic'
                                        if  round((sum(db_Noactive_smile['Class']==Dataset_class[1])/len(db_Noactive_smile))*100,0)<66:   
                                            # classification = 'NON-Active*(10)'
                                             classification = 'NON-Active*: The compounds having the characteristics of several clusters are mainly NO-toxic'
                                    if  treshold_Nocarc<83 and treshold_Nocarc>50: classification = 'Uncertain' 
                                    if treshold_Nocarc>=16 and treshold_Nocarc<=50: 
                                        # classification = 'Active*(10)'
                                         classification = 'Active*: The compounds having the characteristics of several clusters are mainly toxic'
                             
                            # se non hanno composti in comune i CLUSTERS
                            elif classification == 'not define':

                                Noactive_nocarc = round(((sum(db_Noactive_smile['Class']==Dataset_class[1]))/(len(db_Noactive_smile)))*100,0)
                                active_carc = round(((sum(db_active_smile['Class']==Dataset_class[0]))/(len(db_active_smile)))*100,0)
                                
                                if active_carc>=83: 
                                    # classification = 'Active***(11)'
                                    classification = 'Active***: The population of compounds in active clusters is numerically predominant'
                                if active_carc>=66 and Noactive_nocarc>=83: classification = 'Equivocal'
                                if active_carc>=66 and Noactive_nocarc>=66: 
                                    # classification = 'Active*(11)'
                                    classification = 'Active*: The population of compounds in active clusters is numerically predominant'
                                if active_carc>=50 and active_carc<83 and Noactive_nocarc>=83: 
                                    # classification = 'NON-Active*'
                                    classification = 'NON-Active*: The population of compounds in Inctive clusters is numerically predominant'
                                if active_carc>=50 and active_carc<83 and Noactive_nocarc>=66 and Noactive_nocarc<83: 
                                    classification = 'Equivocal'                                                 
                        else:
                            
                            if len(db_merg) != 0:
                                if rateo != 100 or rateo != 1: # queste casistiche sarebbero già considerate
                                    if rateo>75: 
                                        # classification = 'Active**(12)'
                                        classification = 'Active**: The population of compounds in active clusters is numerically predominant'
                                    if rateo<25 and rateo>0: 
                                        # classification = 'NON-Active**(12)'
                                        classification = 'NON-Active**: The population of compounds in active clusters is numerically predominant'
                                    if rateo == 0: classification = 'no group'


                    out = [exception,classification]

                else:
                    exception = []
                    classification = 'No data for groups' 
                    out = [exception,classification]
                
                if 'SMILES' in alert.keys(): alert = alert.drop(['SMILE'],axis=1)
                if len(exception)>0 and len(list(alert.keys()))>0:

                    df_all = self.grp_simili.copy().rename(columns={'Experimental value':'Class'})
                    alt = list(alert.keys()) 
                    
                    for n_alt, alt1 in enumerate(alt):
                            if n_alt == 0:
                                exception_df = df_all.drop(df_all[df_all[alt1]==0].index)
                            else:
                                exception_df = exception_df.drop(exception_df[exception_df[alt1]==0].index)
                    
                    exception_df_lst = []
                    for nexception in exception:
                        df_exp = exception_df.drop(exception_df[exception_df[nexception]==0].index).reset_index(drop=True)
                        exception_df_lst.append(df_exp)

                    exception_df_Carcino_lst = []
                    exception_df_NoCarcino_lst = []

                    for exception_df1 in exception_df_lst:

                        exception_df_Carcino = len(exception_df1.loc[exception_df1['Class']==Dataset_class[0]])
                        exception_df_NoCarcino = len(exception_df1.loc[exception_df1['Class']==Dataset_class[1]])
                        exception_df_Carcino_lst.append(exception_df_Carcino)#,exception_df_NoCarcino])
                        exception_df_NoCarcino_lst.append(exception_df_NoCarcino)
                        out1.append([exception_df_Carcino_lst, exception_df_NoCarcino_lst])

        # se il reasoning non ha portato a nulla: knn               
        if classification in ['not define', 'NO 10 data', 'No data for groups', 'No data']: # 'No data aggiunto 29/09
            df_simil_rdkit_noRedundant.sort_values('Similarity', inplace=True)
            df_simil_rdkit_noRedundant.reset_index(drop=True, inplace = True)
            if len(df_simil_rdkit_noRedundant)>1 and np.mean(df_simil_rdkit_noRedundant['Similarity'][:2])>=0.80:
                if sum(df_simil_rdkit_noRedundant['Class'][:2])>0:
                    classification = 'Active by knn'
                else:
                    classification = 'NON Active by knn'
        
        out = [exception,classification,simili6]  
        if len(out1)>0: out.append(out1) 
        self.prediction = out
        return out

    def reasoning_intero(self, Tox):
        
        self.Calculate_grp()
        # simili a 0.650
        df_simil_rdkit_noRedundant = MoleculeInput.filter_simil(self)
        # calcola gruppi sul target
        target = SingleMolecule.Mol_group_removeRedundant_no_zero(SingleMolecule(self.smile)).drop('SMILES', axis=1)
        # valuta i gruppi di reasoning
        if len(self.simili6)>0: simili6 = self.simili6
        else: 
            MoleculeInput.ortogonal_research(self, Tox)
            simili6 = self.simili6
        # drop caariche, non sono un MG
        if any(pd.Series(target.keys()).apply(lambda x: x in ['charge-', 'charge+'])):
            if 'charge-' in target.keys():
                target = target.drop('charge-',axis=1)
            if 'charge+' in target.keys():
                target = target.drop('charge+',axis=1)
                

        # similarità choice
        # simil_type = 'Grp_Similarity_Emilio'
        # simil_type = 'Grp_Similarity_mean'
        simil_type = 'Similarity'
        
        # riordina in relazione della similarità scelta
        df_simil_rdkit_noRedundant.sort_values(simil_type, ascending=False, inplace = True)
        df_simil_rdkit_noRedundant.reset_index(drop=True, inplace=True)
        value = []
        
        # knn, se serve il reasoning metti yes: se non valogono i criteri per il knn (yes)
        if len(df_simil_rdkit_noRedundant)>2:

            value.append(np.inner(np.power(np.array(df_simil_rdkit_noRedundant[simil_type][:3]),3),
                                np.array(df_simil_rdkit_noRedundant['Experimental value'][:3]))
                        /(np.sum(np.power(np.array(df_simil_rdkit_noRedundant[simil_type][:3]),3))))

            # criteri per il reasoning #0.85, 3 valori pre 6-luglio
            if np.mean(df_simil_rdkit_noRedundant[simil_type][:3])>=0.85 and abs(max(df_simil_rdkit_noRedundant['Experimental value'][:3])-min(df_simil_rdkit_noRedundant['Experimental value'][:3]))<1.5:        
                reasoning = 'no'
            else: 
                reasoning = 'yes'
                          
        else: value.append('no data')

        # valore medio pesato dei singoli grp del reasoning
        simil_value = []
        if len(simili6)>1: # aggiunto 20/10/2022
            for grp in simili6:
                if len(grp)>1 and np.mean(grp[simil_type][:2])>=0.75: # criteri per il reasoning #0.80, 3 valori pre 6-luglio
                    pred_grp = (np.inner(np.power(np.array(grp[simil_type][:3]),3),np.array(grp['Class'][:3])))/(np.sum(np.power(np.array(grp[simil_type][:3]),3)))
                    simil_value.append(pred_grp)

        ########################## modelli lineari ##################################à
        # se i dati sono sufficienti
        if value[0]!='no data':
            
            # aggiungi alert a gruppi
            cv = LeaveOneOut()
            descr1 = SingleMolecule.Alerts(SingleMolecule(self.smile))
            descr = SingleMolecule.Mol_group_removeRedundant(SingleMolecule(self.smile))
            descr.drop(['charge-','charge+'], axis=1, inplace=True)
            descriptor = pd.merge(descr,descr1, on = 'SMILES')


            #modello sui simili a 0,650 cerca fattore di correzione
            X_train = df_simil_rdkit_noRedundant.iloc[:,:-9] # dovrebbe togliere tutte le cose inutili tipo class etc
            y_train = df_simil_rdkit_noRedundant['Experimental value'] # -value[0]# voglio un fattore di correzione

            #modello sui simili con descrittori grp target
            if len(list(target.keys()))>0:
                modello2 = LinearRegression()
                modello2.random_state = 42             
                modello2.fit(X_train[list(target.keys())],y_train)
                scores = cross_val_score(modello2,X_train[list(target.keys())],y_train, 
                                        scoring='neg_mean_absolute_error',
                                        cv=cv, n_jobs=-1)
                modello2_err = np.mean(np.absolute(scores))
                p2=modello2.predict(descriptor[list(target.keys())])
            else: 
                p2 = 'no target grp'
                modello2_err = 'no target grp'
            
            #modello sui simili con descrittori tutti
            modello3 = LinearRegression()
            modello3.random_state = 42           
            modello3.fit(X_train, y_train)
            scores = cross_val_score(modello3, X_train, y_train, scoring='neg_mean_absolute_error',
                        cv = cv, n_jobs=-1)
            modello3_err = np.mean(np.absolute(scores))
            p3 = modello3.predict(descr.drop('SMILES', axis=1))

            # modello sui simili con descrittori molecolari
            X_ = MoleculeInput.descriptor_similarity(self)
            X_train2 = X_.iloc[:, -12:-1]
            modello4 = LinearRegression()
            modello4.random_state = 42
            process = [('num', SimpleImputer(strategy='median')),
                       ('scaling', StandardScaler()),
                       ('model', modello4)]

            my_pipeline = Pipeline(steps = process)                
            my_pipeline.fit(X_train2, y_train)
            target_y_descr = Descriptors_Calculator(self.smile)
            scores = cross_val_score(my_pipeline, X_train2, y_train, scoring = 'neg_mean_absolute_error',
                        cv = cv, n_jobs=-1)
            modello4_err = np.mean(np.absolute(scores))           
            p4 = my_pipeline.predict(target_y_descr)

            err = [modello2_err, modello3_err, modello4_err] 

            if 'no target grp' in err:              
                err.remove('no target grp')  
                min_err = err.index(min(err))             
                if min_err == 0:
                    best_value = p3
                    predictions = modello3.predict(X_train) 
                elif min_err == 1:
                    best_value = p4
                    predictions = my_pipeline.predict(X_train2) 
            else:
                min_err = err.index(min(err))
                if min_err == 0:
                    best_value = p2
                    predictions = modello2.predict(X_train[list(target.keys())])
                elif min_err == 1:
                    best_value = p3
                    predictions = modello3.predict(X_train)           
                elif min_err == 2:
                    best_value = p4
                    predictions = my_pipeline.predict(X_train2) 
                    
            err = round(min(err),2)
            # p = round(p[0],3)
            best_value = round(best_value[0],3)
            #and type(simil_value_pt2[0])==str   
            if simil_value != []: 
                simil_value_medio = round(sum(simil_value)/len(simil_value),3)
                Linear = 'no'
            else: 
                simil_value_medio = 'no data'
                Linear = 'yes'
            if isinstance(err, float):
                if reasoning == 'no': pred_type = 'Knn'
                elif reasoning == 'yes':
                    if Linear == 'no': pred_type = 'Reasoning'
                    elif Linear == 'yes':
                        # Modifica 19/1 err = 0.5 
                        if err < 1: pred_type = 'Linear'
                        else: pred_type = 'Global'
                else:             
                    if Linear == 'no': pred_type = 'Reasoning'
                    elif Linear == 'yes':
                        # Modifica 19/1 prima err = 0.5
                        if err < 1: pred_type = 'Linear'
                        else: pred_type = 'Global' # no data o global?
            else: pred_type = 'Global' # Modifica 19/1 prima era 'no data'

        else:
            reasoning = 'no data'
            err = 'no data'
            best_value = 'no data'
            Linear = 'no data'
            predictions = 'no data' 
            pred_type = 'no data' # Modifica 19/1 prima era 'no data'
            simil_value_medio = 'no data'
            y_train = 'no data'


        out_frame = {'SMILES': self.smile, 
                    'Knn': value[0], 
                    'Grouping Value': simil_value_medio,
                    'LinearModel': best_value,
                    'err': err, 
                    'reasoning': reasoning,
                    'LocalModel': Linear,
                    'pred_type': pred_type,
                    'predictions': predictions,
                    'exp value local': y_train}


        out_frame1 = pd.DataFrame([out_frame])
        self.prediction_intero = out_frame1.copy()  
        return out_frame1

    def reasoning_on_scaffold(self):

        '''valuta se lo scaffold del target è uguale a quello dei simili filter, 
        prima devi calcolarli quindi (filter_simil)'''

        simili_same_scaffold = self.filter_simil.copy()
        simili_same_scaffold['Scaffold'] = simili_same_scaffold['SMILES'].apply(lambda x: scaffold(self.smile, x))
        simili_same_scaffold.drop(simili_same_scaffold[simili_same_scaffold['Scaffold']==0].index, inplace=True)
        simili_same_scaffold.reset_index(drop=True, inplace=True)
        self.scaffold = simili_same_scaffold.copy()
        return simili_same_scaffold

    def descriptor_similarity(self):
        
        target_descr = Descriptors_Calculator(self.smile)

        for j, smi in enumerate(self.filter_simil['SMILES']):

            d = Descriptors_Calculator(smi)
            if j == 0:
                data = d.copy()               
            else:
                data = pd.concat([data, d])
                data.reset_index(drop=True, inplace=True)
            
        
        data_target = pd.concat([target_descr, data])
        data_target.reset_index(drop=True, inplace=True)    
        X = scale(data_target.values)
        X_target_scaled = X[0]

        desc_dist = []
        for i in range(len(X)):
            if i > 0:
                desc_dist.append(distance.cityblock(X_target_scaled, X[i]))
                
        if len(desc_dist) == len(data):
            data['descriptor_similarity'] = desc_dist

        self.simil_descriptor_calculated = data.copy()
        data_out = pd.concat([self.filter_simil, data], axis=1)
        data_out.sort_values('descriptor_similarity', inplace=True)
        data_out.reset_index(drop=True, inplace=True)
        return data_out

# generatore di oggetti della classe principale
# ogni oggetto è definito con la smiles del target e il subset
# di simili VEGA a 0.650  

def Generator_Molecule(df_input, df):

    '''
    Codice per gestire il risultato della matrix similarity,
    non penso ti serva Alessio..
    Filtra per ogni target i simili a 0.650 e con questo subset + la smile del target genero un oggetto della classe
    principale MoleculeInput
    '''
    
    mol = []
    mol_outAD_index = [] #nessun elemento del database è simile oltre il valore soglia
    mol_identity = [] #similarità 1

    for j, i in enumerate(df.keys()):

        if j>0 and j<(len(df.columns)-1): 

            a = df[[i, 'SMILES', 'Experimental value']]
            b = a.drop(a[a[i]<0.65].index).reset_index(drop=True)
            c = a.drop(a[a[i]<1].index).reset_index(drop=True) #cerca identità

        #occhio alcune colonne potrebbero essere vuote e occhio a quelle che contengono l'identità
            if len(b)>0:
                b.rename(columns={i:'Similarity'}, inplace=True)
                mol.append(b.sort_values('Similarity',ascending=False).reset_index(drop=True))
            else: mol_outAD_index.append(i) #append smile

            if len(c)>0: mol_identity.append([c,i])


    if len(mol_outAD_index)>0:strange_molecule=pd.DataFrame([mol_outAD_index])           
    if 'SMILES' not in df_input.keys():df_input.rename(columns={0:'SMILES'}, inplace = True)
    df_input = df_input.drop(df_input[df_input['SMILES'].apply(lambda x: x in mol_outAD_index)].index).reset_index(drop=True)
    
    lst = []
    for j, mol_input in enumerate(df_input['SMILES']):
        mol2 = MoleculeInput(mol_input, mol[j]) # genero oggetto della classse MoleculeInput
        lst.append(mol2)
        
    out = [df_input, lst, mol_outAD_index]
    return out

# function per rimuovere le ridondanze dei MG:
# es: non contare benzene se è un fenolo.. etc
def remove_Redundant(a:Chem.rdchem.Mol, b:pd.DataFrame): #dove b è grp_rdkit_target o Mol_group
            
    #BENZENE 
    c = Fdt.f_benzene(a)
    
    b0 = Fdt.f_Ar_OH(a)
    b1 = Fdt.f_Ar_Cl_Br(a)
    b2 = Fdt.f_aryl_methyl(a)
    b3 = Fdt.f_benzodiazepine(a)
    b4 = Fdt.f_benzCH2(a)
    b5 = Fdt.f_biphenol(a)
    b6 = Fdt.f_Ar_OR(a)
    b7 = Fdt.f_Ar_R(a)
    b8 = Fdt.f_op_diphenolo_OR(a)
    b9 = Fdt.f_Ar_COR(a)
    b10 = Fdt.f_Ar_COO_R_H(a)
    b11 = Fdt.f_C_3phenols(a)
    b12 = Fdt.f_Ring_3OH_3OR(a)
    b13 = Fdt.f_aniline(a)
    
    t = [b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13] 
    
    del_benz = []
    for benz in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(benz))>0:
                        del_benz.append(benz)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(benz))>0:
                                del_benz.append(benz)
    del_benz_set = set(del_benz)
    b['benzene'] = b['benzene']-len(del_benz_set)
    
    #amide urea
    if 'amide' in b.keys():
        c = Fdt.f_amide(a) 
        
        b0 = Fdt.f_urea(a)
        
        t = [b0]
        del_amide = []
        for amide in c:
            for j, k in enumerate(t):
                if len(k)>0 and type(k) is tuple:
                    for i in range(len(k)):
                        if len(set(k[i]).intersection(amide))>0:
                            del_amide.append(amide)
                elif len(k)>0 and type(k) is list:
                    for i in range(len(k)): 
                        if len(k[i])>0:
                            for z in k[i]:
                                if len(set(z).intersection(amide))>0:
                                    del_amide.append(amide)
        del_amide_set = set(del_amide)
        b['amide'] = b['amide']-len(del_amide_set)
    
                    
    #NH0 NH1 NH2
    c = Fdt.f_NH0(a)
    c1 = Fdt.f_NH1(a)
    c2 = Fdt.f_NH2(a)

    b1 = Fdt.f_pyridine(a) 
    #b2 = Fdt.f_piperdine(a)
    b3 = Fdt.f_nitro(a)
    b4 = Fdt.f_nitro_arom_nonortho(a)
    b5 = Fdt.f_amide(a)
    b6 = Fdt.f_urea(a)
    b7 = Fdt.f_hdrzine(a)
    b8 = Fdt.f_aniline_term(a)
    b9 = Fdt.f_imidazole(a)


    #t = [b1,b2,b3,b4,b5,b6,b7,b8,b9]
    t = [b1,b3,b4,b5,b6,b7,b8,b9]

    del_NH0 = []
    del_NH1 = []
    del_NH2 = []

    for NH0 in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(NH0))>0:
                        del_NH0.append(NH0)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(NH0))>0:
                                del_NH0.append(NH0)
    del_NH0_set = set(del_NH0)
    b['NH0'] = b['NH0']-len(del_NH0_set)
    
    for NH1 in c1:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(NH1))>0:
                        del_NH1.append(NH1)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(NH1))>0:
                                del_NH1.append(NH1)
    del_NH1_set = set(del_NH1)
    b['NH1'] = b['NH1']-len(del_NH1_set)
    
    for NH2 in c2:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(NH2))>0:
                        del_NH2.append(NH2)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(NH2))>0:
                                del_NH2.append(NH2)
    del_NH2_set = set(del_NH2)
    b['NH2'] = b['NH2']-len(del_NH2_set)
    
    # aniline

    c = Fdt.f_aniline(a)

    b1 = Fdt.f_amide(a)
    b2 = Fdt.f_urea(a)
    b3 = Fdt.f_carbamate(a)
    b4 = Fdt.f_aniline_term(a)

    t = [b1,b2,b3,b4]
    
    del_anilina = []
    
    for anilina in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(anilina))>0:
                        del_anilina.append(anilina)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(anilina))>0:
                                del_anilina.append(anilina)
    del_anilina_set = set(del_anilina)
    b['aniline'] = b['aniline']-len(del_anilina_set)

    # priamide

    c = Fdt.f_priamide(a)

    b1 = Fdt.f_amide(a)
    b2 = Fdt.f_urea(a)
    b3 = Fdt.f_carbamate(a)
    b4 = Fdt.f_aniline_term(a)

    t = [b1,b2,b3,b4]
    
    del_priamide = []
    
    for priamide in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(priamide))>0:
                        del_priamide.append(priamide)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(priamide))>0:
                                del_priamide.append(priamide)
    del_priamide_set = set(del_priamide)

    b['priamide'] = b['priamide']-len(del_priamide_set)



    # ketone 1/16

    c = Fdt.f_ketone(a)

    b1 = Fdt.f_ketone_aliphatic(a)
    b2 = Fdt.f_Ar_ketone(a)


    t = [b1,b2]
    
    del_ketone = []
    
    for ketone in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(ketone))>0:
                        del_ketone.append(ketone)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(ketone))>0:
                                del_ketone.append(ketone)
    del_ketone_set = set(del_ketone)

    b['ketone'] = b['ketone']-len(del_ketone_set)

    # Ar_COR
    c = Fdt.f_ketone(a)
    c1 = Fdt.f_ketone_Topliss(a)

    b1 = Fdt.f_pyridine(a)


    t = [b1]

    del_ketone = []
    del_ketone_Topliss = []
    
    for ketone in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(ketone))>0:
                        del_ketone.append(ketone)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(ketone))>0:
                                del_ketone.append(ketone)
    del_ketone_set = set(del_ketone)
    b['ketone'] = b['ketone']-len(del_ketone_set)
    
    for ketone_Topliss in c1:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(ketone_Topliss))>0:
                        del_ketone_Topliss.append(ketone_Topliss)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(ketone_Topliss))>0:
                                del_ketone_Topliss.append(ketone_Topliss)
    del_ketone_Topliss_set = set(del_ketone_Topliss)
    b['ketone_Topliss'] = b['ketone_Topliss']-len(del_ketone_Topliss_set)

    
    #quatN == NH2?
    b['NH2'] = b['NH2'] + b['quatN']

    # togli methoxy se estere 
    # 29/09
    c = Fdt.f_methoxy(a)

    b1 = Fdt.f_ester(a)

    t = [b1]
    del_methoxy = []
    for methoxy in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(methoxy))>0:
                        del_methoxy.append(methoxy)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(methoxy))>0:
                                del_methoxy.append(methoxy)                 
    del_methoxy_set = set(del_methoxy)
    b['methoxy'] = b['methoxy']-len(del_methoxy_set)


    
    # togli etere se estere e methoxy 

    c = Fdt.f_ether(a)

    b1 = Fdt.f_ester(a)
    # 29/09
    b2 = Fdt.f_methoxy(a)

    # t = [b1]
    t = [b1, b2]
    del_ether = []
    for ether in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(ether))>0:
                        del_ether.append(ether)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(ether))>0:
                                del_ether.append(ether)                 
    del_ether_set = set(del_ether)
    b['ether'] = b['ether']-len(del_ether_set)
    
    #phenolo e Ar_OH
    
    c = Fdt.f_Ar_OH(a) #sottrai a c bn

    b1 = Fdt.f_phenol(a)

    t = [b1]
    
    del_Ar_OH = []
    for Ar_OH in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(Ar_OH))>0:
                        del_Ar_OH.append(Ar_OH)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(Ar_OH))>0:
                                del_Ar_OH.append(Ar_OH)                 
    del_Ar_OH_set = set(del_Ar_OH)
    b['Ar_OH'] = b['Ar_OH']-len(del_Ar_OH_set)
    
    # Biclycle 24 Nov 23 ERIKA

    c = Fdt.f_bicyclic(a)

    b1 = Fdt.f_Ar_bicycle(a)
    b2 = Fdt.f_Al_bicycle(a)


    t = [b1,b2]
    
    del_bicycle = []
    
    for bicycle in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(bicycle))>0:
                        del_bicycle.append(bicycle)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(bicycle))>0:
                                del_bicycle.append(bicycle)
    del_bicycle_set = set(del_bicycle)

    b['bicyclic'] = b['bicyclic']-len(del_bicycle_set)

    
    #charge remove in N03
    c = Fdt.f_charge2(a)
    c1 = Fdt.f_charge1(a)

    b1 = Fdt.f_nitro(a)
    b2 = Fdt.f_nitro_arom_nonortho(a)


    t = [b1,b2]

    del_charge1 = []
    del_charge2 = []  
    


    for charge1 in c:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(charge1))>0:
                        del_charge1.append(charge1)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(charge1))>0:
                                del_charge1.append(charge1)
    del_charge1_set = set(del_charge1)
    b['charge-'] = b['charge-']-len(del_charge1_set)
    
    for charge2 in c1:
        for j, k in enumerate(t):
            if len(k)>0 and type(k) is tuple:
                for i in range(len(k)):
                    if len(set(k[i]).intersection(charge2))>0:
                        del_charge2.append(charge2)
            elif len(k)>0 and type(k) is list:
                for i in range(len(k)): 
                    if len(k[i])>0:
                        for z in k[i]:
                            if len(set(z).intersection(charge2))>0:
                                del_charge2.append(charge2)
    del_charge2_set = set(del_charge2)
    b['charge+'] = b['charge+']-len(del_charge2_set)

    #if b.keys().any() in list_out: 
    b.drop(list_out, axis=1, inplace=True)
    # test 29/09  prima non c'era 
    num = b._get_numeric_data()
    num[num < 0] = 0
    
    return b 
        
def cassification(a): 
    if a > 65: 
        classification='Tox'
    elif a<35 and a !=-1:
        classification='NoTox'
    else:
        if a != -1:
            classification='Neutro'
        else:
            classification='No Data'
    return classification

def prevalence_threshold(Dataset, *args):

    '''def prevalence for endpoint with continous value'''
    
    if len(args)>1:
        thr1 = args[0]
        thr2 = args[1]

    Dataset_Mol = []
    for mol_input in (Dataset['SMILES']):
        mol2 = SingleMolecule(mol_input)
        Dataset_Mol.append(mol2)
    for i, smile in enumerate(Dataset_Mol):
        if i==0:
            df_groups = SingleMolecule.Mol_group_removeRedundant(smile)
        else:
            df_groups = df_groups.append(SingleMolecule.Mol_group_removeRedundant(smile))
    df_groups = df_groups.reset_index(drop=True)
    df_groups = df_groups.assign(SMILES = Dataset['SMILES'], Class = Dataset['Experimental value'])
    df_all = df_groups.copy()
    Dataset_active = df_all.drop(df_all[df_all['Class']<thr2].index)
    Dataset_active['Class'] = 1
    Dataset_Noactive = df_all.drop(df_all[df_all['Class']>thr1].index)
    Dataset_Noactive['Class'] = 0
    dic = {'GROUPS':list(Dataset_active.keys()),'Active':list(Dataset_active.sum()),
           'NoActive':list(Dataset_Noactive.sum())}
    frame = pd.DataFrame(dic)
    frame = frame.drop(frame[frame['GROUPS'].apply(lambda x: x in ['Class','SMILES'])].index)
    frame['Prevalence%'] = (frame['Active']/(frame['Active']+frame['NoActive']))*100
    frame['Prevalence%'] = round(frame['Prevalence%'].fillna(-1),3)
    frame['Tox'] =  frame['Prevalence%'].apply(lambda x: cassification(x))

# Aggiunta del 24 Nov 23
    Tox = frame.copy()
    ac = len(Dataset_active)
    noac = len(Dataset_Noactive)
    ToxClu = Tox.copy()
    ToxClu['ActiveScale'] = ToxClu['Active']/ac
    ToxClu['NOActiveScale'] = ToxClu['NoActive']/noac
    ToxClu['Prevalence%Scale'] = (ToxClu['ActiveScale']/(ToxClu['ActiveScale']+ToxClu['NOActiveScale']))*100
    ToxClu['ToxScale'] = ToxClu['Prevalence%Scale'].apply(lambda x: cassification(x))
    
    return ToxClu


def prevalence(Dataset):

    Dataset_Mol = []
    for mol_input in (Dataset['SMILES']):
        mol2 = SingleMolecule(mol_input)
        Dataset_Mol.append(mol2)
    for i, smile in enumerate(Dataset_Mol):
        if i==0:
            df_groups = SingleMolecule.Mol_group_removeRedundant(smile)       
        else:
            df_groups = df_groups.append(SingleMolecule.Mol_group_removeRedundant(smile))

    df_groups = df_groups.reset_index(drop=True)
    df_groups = df_groups.assign(SMILES = Dataset['SMILES'],Class = Dataset['Experimental value'])
    df_all = df_groups.copy()

    Dataset_active = df_all.loc[df_all['Class']==Dataset_class[0]].reset_index(drop=True)
    Dataset_Noactive = df_all.loc[df_all['Class']==Dataset_class[1]].reset_index(drop=True)

    dic = {'GROUPS':list(Dataset_active.keys()),'Active':list(Dataset_active.sum()),'NoActive':list(Dataset_Noactive.sum())}
    frame = pd.DataFrame(dic)
    frame = frame.drop(frame[frame['GROUPS'].apply(lambda x: x in ['Class','SMILES'])].index)
    frame['Prevalence%'] = (frame['Active']/(frame['Active']+frame['NoActive']))*100
    frame['Prevalence%'] = round(frame['Prevalence%'].fillna(-1),3)
    frame['Tox'] =  frame['Prevalence%'].apply(lambda x: cassification(x))
    Tox = frame.copy()
    ac = len(Dataset_active)
    noac = len(Dataset_Noactive)
    ToxClu = Tox.copy()
    ToxClu['ActiveScale'] = ToxClu['Active']/ac
    ToxClu['NOActiveScale'] = ToxClu['NoActive']/noac
    ToxClu['Prevalence%Scale'] = (ToxClu['ActiveScale']/(ToxClu['ActiveScale']+ToxClu['NOActiveScale']))*100
    ToxClu['ToxScale'] = ToxClu['Prevalence%Scale'].apply(lambda x: cassification(x))
    
    return ToxClu

def prevalence_alerts(Dataset):

    '''def prevalence for alerts'''
    
    Dataset_Mol = []
    
    for mol_input in (Dataset['SMILES']):
        mol2 = SingleMolecule(mol_input)
        Dataset_Mol.append(mol2)
        
    for i, smile in enumerate(Dataset_Mol):
            if i==0:
                df_groups = SingleMolecule.Alerts(smile)       
            else:
                df_groups = df_groups.append(SingleMolecule.Alerts(smile) )

    df_groups = df_groups.reset_index(drop=True)
    df_groups = df_groups.assign(SMILES = Dataset['SMILES'],Class = Dataset['Experimental value'])
    df_all = df_groups.copy()

    Dataset_active = df_all.loc[df_all['Class']==Dataset_class[0]].reset_index(drop=True)
    Dataset_Noactive = df_all.loc[df_all['Class']==Dataset_class[1]].reset_index(drop=True)

    dic = {'GROUPS':list(Dataset_active.keys()),'Active':list(Dataset_active.sum()),'NoActive':list(Dataset_Noactive.sum())}
    frame = pd.DataFrame(dic)
    frame = frame.drop(frame[frame['GROUPS'].apply(lambda x: x in ['Class','SMILES'])].index)
    frame['Prevalence%'] = (frame['Active']/(frame['Active']+frame['NoActive']))*100
    frame['Prevalence%'] = round(frame['Prevalence%'].fillna(-1),3)
    frame['Tox'] =  frame['Prevalence%'].apply(lambda x: cassification(x))
    Tox = frame.copy()
    ac = len(Dataset_active)
    noac = len(Dataset_Noactive)
    ToxClu = Tox.copy()
    ToxClu['ActiveScale'] = ToxClu['Active']/ac
    ToxClu['NOActiveScale'] = ToxClu['NoActive']/noac
    ToxClu['Prevalence%Scale'] = (ToxClu['ActiveScale']/(ToxClu['ActiveScale']+ToxClu['NOActiveScale']))*100
    ToxClu['ToxScale'] = ToxClu['Prevalence%Scale'].apply(lambda x: cassification(x))

    return ToxClu

def prevalence_threshold_alerts(Dataset, *args):

    '''def prevalence for endpoint with continous value'''
    
    if len(args)>1:
        
        thr1 = args[0]
        thr2 = args[1]

    Dataset_Mol = []
    
    for mol_input in (Dataset['SMILES']):
        mol2 = SingleMolecule(mol_input)
        Dataset_Mol.append(mol2)
        
    '''    
    for i, smile in enumerate(Dataset_Mol):
        if i==0:
            df_groups = SingleMolecule.Alerts(smile)
        else:
            df_groups = df_groups.append(SingleMolecule.Alerts(smile))
    '''    
    for i, smile in enumerate(Dataset_Mol):
        d = SingleMolecule.Alerts(smile)
        if i==0:
            df_groups = d.copy()
        else:
            df_groups = pd.concat([df_groups, d])
            
    df_groups = df_groups.reset_index(drop=True)
    
    # aggiunta del 13 luglio
    Dataset.reset_index(drop=True, inplace=True)
    
    df_groups = df_groups.assign(SMILES = Dataset['SMILES'], Class = Dataset['Experimental value'])
    df_all = df_groups.copy()
    
    # minori di -1 tolti: rimangono i maggiori di -1 che sono tox
    Dataset_active = df_all.drop(df_all[df_all['Class']<thr2].index)
    Dataset_active['Class'] = 1
    # maggiori di -2 tolti: rimangono i minori di -2 che sono no tox
    Dataset_Noactive = df_all.drop(df_all[df_all['Class']>thr1].index)
    Dataset_Noactive['Class'] = 0
    
    dic = {'GROUPS':list(Dataset_active.keys()),'Active':list(Dataset_active.sum()),
           'NoActive':list(Dataset_Noactive.sum())}
    
    frame = pd.DataFrame(dic)
    frame = frame.drop(frame[frame['GROUPS'].apply(lambda x: x in ['Class','SMILES'])].index)
    frame['Prevalence%'] = (frame['Active']/(frame['Active']+frame['NoActive']))*100
    frame['Prevalence%'] = round(frame['Prevalence%'].fillna(-1),3)
    frame['Tox'] =  frame['Prevalence%'].apply(lambda x: cassification(x))
    return frame

# questa funzione non ti serve
def view_difference(smile1, smile2):

    mol1 = Chem.MolFromSmiles(smile1)
    mol2 = Chem.MolFromSmiles(smile2)
    
    
    mcs = rdFMCS.FindMCS([mol1,mol2])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    
    #evidenzia atomi della molecola 1 che non fanno parte del MCS
    match1 = mol1.GetSubstructMatch(mcs_mol)
    target_atm1 = []
    for atom in mol1.GetAtoms():
        if atom.GetIdx() not in match1:
            target_atm1.append(atom.GetIdx())
            
    #evidenzia atomi della molecola 2 che non fanno parte del MCS        
    match2 = mol2.GetSubstructMatch(mcs_mol)    
    target_atm2 = []
    for atom in mol2.GetAtoms():
        if atom.GetIdx() not in match2:
            target_atm2.append(atom.GetIdx())
            
    return Draw.MolsToGridImage([mol1, mol2],highlightAtomLists=[target_atm1, target_atm2])

def Descriptors_Calculator(smile):
    
    ''' some specific descriptors'''

    m = Chem.MolFromSmiles(smile)
    LogP = Descriptors.MolLogP(m)
    MaxPartialCharge = Descriptors.MaxPartialCharge(m)
    MinPartialCharge = Descriptors.MinPartialCharge(m)
    NumHAcceptors = Descriptors.NumHAcceptors(m)
    NumHDonors = Descriptors.NumHDonors(m)
    NumRotatableBonds = Descriptors.NumRotatableBonds(m)
    NumAliphaticRings = Descriptors.NumAliphaticRings(m)
    NumAromaticCarbocycles = Descriptors.NumAromaticCarbocycles(m)
    NumAromaticHeterocycles = Descriptors.NumAromaticHeterocycles(m)
    NumAromaticRings = Descriptors.NumAromaticRings(m)
    ExactMolWt = Descriptors.ExactMolWt(m)

    NameDes = ['LogP','MaxPartialCharge','MinPartialCharge','NumHAcceptors','NumHDonors',
                'NumRotatableBonds','NumAliphaticRings','NumAromaticCarbocycles',
                'NumAromaticHeterocycles','NumAromaticRings','ExactMolWt']

    Desr = [LogP,MaxPartialCharge,MinPartialCharge,NumHAcceptors,NumHDonors,
                NumRotatableBonds,NumAliphaticRings,NumAromaticCarbocycles,
                NumAromaticHeterocycles,NumAromaticRings,ExactMolWt]
    
    dic = {'LogP':LogP,
            'MaxPartialCharge':MaxPartialCharge,
            'MinPartialCharge':MinPartialCharge,
            'NumHAcceptors':NumHAcceptors,
            'NumHDonors':NumHDonors,
            'NumRotatableBonds':NumRotatableBonds,
            'NumAliphaticRings':NumAliphaticRings,
            'NumAromaticCarbocycles':NumAromaticCarbocycles,
            'NumAromaticHeterocycles':NumAromaticHeterocycles,
            'NumAromaticRings':NumAromaticRings,
            'ExactMolWt':ExactMolWt
            }
    
    Descriptors_result = pd.DataFrame([dic])
    
    return Descriptors_result

def all_rdkit_descriptors(smile):

    ''' calculation of all descriptors present in rdkit'''

    name = []

    for nm,calc in Descriptors.descList:
        name.append(nm)
        
    alldescrs = []
    for nm,calc in Descriptors.descList:
        try: alldescrs.append(calc(Chem.MolFromSmiles(smile)))
        except: alldescrs.append(-99999)
    data = pd.DataFrame(alldescrs).T
    data.columns = name
    return data

def DataCuration(Dataset: pd.DataFrame):

    ''' take smiles in input and canonize, kekulaize and define tautomers. Delate duplicate smiles'''
    
    #This variable counts how many chemicals have been treated 
    number_molecule = 0
    #In order to prevent a combinatorial explosion do limit the number of tautomers
    MAX_TAUTOMERS = 500
    file_output = []
    error = []
    for line in Dataset['SMILES']:   
        input_line = line
        try: 
            number_molecule = number_molecule + 1
            input_line = input_line.strip()
            columns = input_line
            import standardiser
            mol = Chem.MolFromSmiles(columns, sanitize=True)
            from molvs import tautomer
            function_canonicaliser = tautomer.TautomerCanonicalizer(max_tautomers=MAX_TAUTOMERS)
            mol4 = function_canonicaliser(mol)       
            Chem.Kekulize(mol4)
            clean_smiles = Chem.MolToSmiles(mol4,isomericSmiles=False, kekuleSmiles=True)
            clean_smiles = str(clean_smiles)
            file_output.append(clean_smiles)
        except: 
            error.append(input_line)

    Data_error = Dataset.drop(Dataset[Dataset['SMILES'].apply(lambda x: x not in error)].index)
    Dataset.drop(Dataset[Dataset['SMILES'].apply(lambda x: x in error)].index, inplace=True)
    
    if len(Dataset['SMILES']) == len(file_output):
        Dataset['SMILES_clean'] = file_output
    else:print('error')

    # Dataset.drop_duplicates('SMILES_clean', inplace=True)
    Dataset.reset_index(drop=True, inplace=True)
    
    smile_can = []
    for smile in Dataset['SMILES']:
        smile_can.append(Chem.CanonSmiles(smile))

    if len(Dataset['SMILES']) == len(smile_can):
        Dataset['SMILES_old'] =  Dataset['SMILES']
        Dataset['SMILES'] = smile_can

    Dataset.drop_duplicates('SMILES', inplace=True)
    Dataset.drop_duplicates('SMILES_clean', inplace=True)
    Dataset.reset_index(drop=True, inplace=True)
    return [Dataset,Data_error]

def DataCurationTarget(Dataset: pd.DataFrame):

    '''curation like Datacurtation functions but not delate duplicates'''
    
    #This variable counts how many chemicals have been treated 
    number_molecule = 0
    #In order to prevent a combinatorial explosion do limit the number of tautomers
    MAX_TAUTOMERS = 500
    file_output = []
    error = []
    for line in Dataset['SMILES']:   
        input_line = line
        try: 
            number_molecule = number_molecule + 1
            input_line = input_line.strip()
            columns = input_line
            import standardiser
            mol = Chem.MolFromSmiles(columns, sanitize=True)
            from molvs import tautomer
            function_canonicaliser = tautomer.TautomerCanonicalizer(max_tautomers=MAX_TAUTOMERS)
            mol4 = function_canonicaliser(mol)       
            Chem.Kekulize(mol4)
            clean_smiles = Chem.MolToSmiles(mol4,isomericSmiles=False, kekuleSmiles=True)
            clean_smiles = str(clean_smiles)
            file_output.append(clean_smiles)
        except: 
            error.append(input_line)

    Data_error = Dataset.drop(Dataset[Dataset['SMILES'].apply(lambda x: x not in error)].index)
    Dataset.drop(Dataset[Dataset['SMILES'].apply(lambda x: x in error)].index, inplace=True)
    if len(Dataset['SMILES']) == len(file_output):
        Dataset['SMILES_clean'] = file_output
    else:print('error')
    # Dataset.drop_duplicates('SMILES_clean', inplace=True)
    Dataset.reset_index(drop=True, inplace=True)
    
    smile_can = []
    for smile in Dataset['SMILES']:
        smile_can.append(Chem.CanonSmiles(smile))

    if len(Dataset['SMILES']) == len(smile_can):
        Dataset['SMILES_old'] =  Dataset['SMILES']
        Dataset['SMILES'] = smile_can
        # Dataset.drop_duplicates('SMILES', inplace=True)
        Dataset.reset_index(drop=True, inplace=True)
    return [Dataset,Data_error]

def scaffold(smile1, smile2):

    'If 2 molecule have the same scaffold return 1'
    try:
        mol1 = Chem.MolFromSmiles(smile1)
        mol2 = Chem.MolFromSmiles(smile2)
        target_scaffold = MS.GetScaffoldForMol(mol1)
        similar_scaffold = MS.GetScaffoldForMol(mol2)
        
        scaffold_target_smiles = Chem.CanonSmiles(Chem.MolToSmiles(target_scaffold))
        scaffold_similar_smiles = Chem.CanonSmiles(Chem.MolToSmiles(similar_scaffold))
        
        if scaffold_target_smiles == scaffold_similar_smiles:
            out = 1
        else:
            out =0
    except:
        out =0 
     
def scaffold_index(smile):
    '''atoms and bonds index for output, index of atoms and bond 
        of the scaffold in molecule similar to target
    '''
    # mol = Chem.MolFromSmiles(self.smile)
    
    mol1 = Chem.MolFromSmiles(smile)
    scaffold = MS.GetScaffoldForMol(mol1)
    index = mol1.GetSubstructMatches(scaffold)

    hit_atss = list(index)
    hit_bondss = []
    for hit_ats in hit_atss:
        hit_bonds = []
        for bond in scaffold.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(mol1.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hit_bondss.append(hit_bonds)
    list_hit_atss = np.array(hit_atss).reshape(1,-1).tolist()[0]
    highlightBonds = np.array(hit_bondss).reshape(1,-1).tolist()[0]
    
    return (highlightBonds, list_hit_atss, scaffold)
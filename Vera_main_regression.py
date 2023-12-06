# rdkit packages
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True 

# Data Analisys library
import pandas as pd

# standard library
import subprocess
import os
from tqdm import tqdm

# personal library and functions
from Function_Vera import*
import VeraGlobalModelForContinuoRF as RF
import reliability
from pdf_out import*

import threading
import warnings
warnings.filterwarnings('ignore')

class myThread (threading.Thread):
   
    def __init__(self, threadID, name, Dataset, df_input, path_program):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.dataset = Dataset
        self.input = df_input
        self.path_program = path_program

    def make_matrix(self):        
        self.dataset['SMILES'].to_csv('dataset1.csv',header = None,index=False)
        self.input['SMILES'].to_csv('dataset2.csv',header = None,index=False)
        command = ["java", "-jar", self.path_program, "-similarity_ds", "dataset1.csv","dataset2.csv"]
        subprocess.call(command)
        df = pd.read_table('similarity_dataset1.csv_dataset2.csv.csv').drop('#', axis=1)
        os.remove('similarity_dataset1.csv_dataset2.csv.csv')
        os.remove('dataset1.csv')
        os.remove('dataset2.csv')
        df.rename(columns={'Mol':'SMILES'}, inplace=True)
        a = pd.merge(df, self.dataset, how='inner', on='SMILES')
        df = a.copy()
        return df

    def run(self):
        print ("Starting " + self.name)
        global df
        df = myThread.make_matrix(self)
        print ("Exiting " + self.name)

def load_input():
    # Ask the user to input the file path
    file_path = input("\n\nPlease enter the TARGET file name (.xlsx)\n\n\tNote the test file must be in Data_Test_ReadAcross: ")

    # Check if the file path is valid
    path = os.path.join(os.getcwd(), "Data_Test_ReadAcross", file_path)
    if os.path.exists(path) and path.lower().endswith('.xlsx'):
        # Load the CSV file using pandas
        df = pd.read_excel(path)

        # Display the loaded data or perform further processing
        print("CSV file loaded successfully. Here are the first few rows:")
        print(df.head())
        return df
    else:
        print("Invalid file path or file is not a CSV.")

def output_name():
    # Ask the user to input the file path
    file_path = input("\n\nPlease enter the name of output file (.xlsx)\n\tNote the output file will be located in result folder: ")

    # Check if the file path is valid
    path = os.path.join(os.getcwd(), "result", file_path)
    if os.path.exists(path) and path.lower().endswith('.xlsx'):
        # Display the loaded data or perform further processing
        print(f"{file_path} file is already present in folder. The old file will be overwrite..)")
    elif not path.lower().endswith('.xlsx'): 
        path = os.path.join(os.getcwd(), "result", f"{file_path}.xlsx")
    
    return path

def get_user_choice():
    print("Do you want a pdf report for each target molecule?\n Note: Time demanding.")
    pdf = None
    while not pdf:
        choice = input("y/n: ")
        if choice == 'y' or  choice == 'yes':
            pdf = 'yes'
        elif choice == 'n' or  choice == 'no':
            pdf = 'no'
        else:
            print("Invalid choice. Please enter y or n.")
    return pdf
    
if __name__ == "__main__":

    df_input = load_input()
    path = output_name()
    pdf = get_user_choice()

    # AcuteTox
    Dataset_path = 'Data_Training_ReadAcross\FishAcuteTox.xlsx'
    Dataset = pd.read_excel(Dataset_path, index_col='Unnamed: 0')

    

    path_program = os.getcwd()
    path_program = os.path.join(path_program, 'CustomFeaturesCLI')
    path_program = os.path.join(path_program, 'CustomFeaturesCLI.jar')

    thread1 = myThread(1, "Thread-1", Dataset, df_input, path_program)
    thread1.start()
    thread1.join()

    if 'SMILES_clean' in df.keys(): del df['SMILES_clean']
    if 'SMILES_Can' in df.keys(): del df['SMILES_Can']
    if 'SMILES_old' in df.keys(): del df['SMILES_old'] # 7/04 
    if 'Unnamed: 0' in df.keys(): del df['Unnamed: 0']


    if len(set(Dataset['Experimental value']))<=2 and all(Dataset['Experimental value'].apply(lambda x: x in [0,1])):
        target_type ='class'
    else:
        target_type ='continuo'

    # prediction 
    """
    if target_type == 'class':
        print("ClassificationProblem")
        obj = Generator_Molecule(df_input, df)
        #global tox
        Tox = prevalence(Dataset)
        prediction = []
        lst = obj[1]
        for j, i in enumerate(tqdm(lst)):
            grp = SingleMolecule.Mol_group_no_zero_removeRedundant(SingleMolecule(i.smile)).drop('SMILES', axis=1)
            if list(grp.keys()) == 0:
                # nel caso non si siano MG (ridondante con il codice di reasoning, mi seriva per velocizzare)
                prediction.append('No MG find in target')
            else:
                c1 = MoleculeInput.Calculate_grp(i)
                c = MoleculeInput.filter_simil(i)
                if len(c)<2: prediction.append('No similar substance in dataset for make a prediction')
                else:
                    b = MoleculeInput.ortogonal_research(i ,Tox)
                    a = MoleculeInput.reasoning(i)
                    if len(a[0])==0:prediction.append(a[1])
                    else: prediction.append('NON Active for exception')

        b = pd.Series(prediction)
        df_input2 = obj[0].copy()
        df_input2['Pred'] = b
        mol_outAD_index = obj[2]
        strange_molecule = pd.Series(mol_outAD_index)
        if len(mol_outAD_index)>0:
            strange_molecule = pd.DataFrame(strange_molecule.transpose(), columns=['SMILES'])
            strange_molecule['SMILES_old'] = 'na'
            strange_molecule['SMILES_clean'] = 'na'
            strange_molecule['Pred'] = 'no similar in database'
            df_input2 = df_input2.append(strange_molecule)
        df_input2.to_excel('result/prediction.xlsx')
        new2 = df_input2.copy()

    """
    if target_type == 'continuo':
        Threshold1 = -2.01 # exp <-2 no tox
        Threshold2 = -2 # exp >-1 tox
        obj = Generator_Molecule(df_input, df)
        # global Tox
        Tox = prevalence_threshold(Dataset, Threshold1, Threshold2)
        lst = obj[1]

        for j, mol in enumerate(tqdm(lst)):
            d =  MoleculeInput.reasoning_intero(mol,Tox)   
            if j == 0:
                data = d.copy()
            else:
                data = pd.concat([data, d])

        data.reset_index(drop=True, inplace=True)
        # data.rename(columns={'pred_medio_grp1':'Grouping Value'}, inplace=True)
        GB = RF.GlobalModel(Dataset, data)
        GB_pred = GB.pipeline_rf()

        data['GBPred'] = GB_pred
        # knn grp(se linear model is not top)
        data_pred = data.drop(data[data['Knn']=='no data'].index) # togli quelli che non hanno almeno 3 simili

        # data_GB = data.drop(data[data['Knn']!='no data'].index)

        new = data_pred.drop(data_pred[data_pred['reasoning']=='yes'].index)
        new['Best Pred'] = new['Knn']
        # new.rename(columns={'Knn':'Best Pred'},inplace = True)

        # need reasoning?
        correction = data_pred.drop(data_pred[(data_pred['reasoning']=='no')].index)

        # correction by grp
        correction_grp = correction.drop(correction[correction['err']<=0.5].index)
        correction_grp = correction.copy()
        correction_grp.drop(correction_grp[correction_grp['Grouping Value']=='no data'].index, inplace=True)
        correction_grp.drop(correction_grp[correction_grp['LocalModel']=='yes'].index, inplace=True)
        correction_grp['Best Pred'] = correction_grp['Grouping Value']

        # correction by linear best value
        correction.drop(correction[correction['err']>0.5].index, inplace=True)
        correction.drop(correction[correction['LocalModel']=='no'].index, inplace=True)
        std_Exp = Dataset['Experimental value'].std()
        mean_Exp = Dataset['Experimental value'].mean()
        correction.drop(correction[correction['LinearModel']>mean_Exp+std_Exp].index,inplace=True)
        correction.drop(correction[correction['LinearModel']<mean_Exp-std_Exp].index,inplace=True)
        correction['Best Pred'] = correction['LinearModel']

        new1 = pd.concat([new, correction, correction_grp])
        data_GB = data.drop(new1.index)
        data_GB['Best Pred'] = data_GB['GBPred']
        # data_GB['pred_type'] = 'Global'

        new2 = pd.concat([new1, data_GB])
        new2['class'] = 0
        new2['class'].loc[list(correction.index)]=1
        new2['class'].loc[list(correction_grp.index)]=2
        new2['class'].loc[list(data_GB.index)]=3
        new2.sort_index(inplace=True)
        stat_obj = Output(Dataset, df_input, Tox, new2)
        std, mean, top, bottom, limit = stat_obj.statistiche()
        reliability = reliability.reliabilityContinuo(new2, top, bottom)
        new2['Reliability'] = reliability
        # con riserva sto pezzo
        new2['pred_type'] = ['Global' if i == 'no data' else i for i in new2['pred_type']]
        new2['pred_type'] = ['Global' if round(row['Best Pred'], 3) == round(row['GBPred'], 3) else row['pred_type'] for _, row in new2.iterrows()] 
        
        # new2.loc[:, ['SMILES', 'Best Pred', 'Reliability']].to_excel(path)
        new2.to_excel(path)
    
    if pdf=="yes":
        print("start pdf generation...")
        pdf_generator = Output(Dataset, df_input, Tox, new2)
        pdf_generator.makepdf(lst)
        # new2.loc[:, ['SMILES', 'Best Pred', 'Reliability']].to_excel(path)
        new2.to_excel(path)
        print("Done: you can find the results in result folder.")
# importing packages

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from rdkit import Chem 


class GlobalModel:
    def __init__(self, data, target):

        ''' self.dataset from self.datacuration '''

        self.data = data
        self.target = target

    def all_rdkit_descriptors(self, smile):

        ''' calculation of all descriptors present in rdkit'''

        name = []

        for nm,calc in Descriptors.descList[:-85]: # tolgo gruppi funzionali
            name.append(nm)       
        alldescrs = []
        for nm,calc in Descriptors.descList[:-85]:
            alldescrs.append(calc(Chem.MolFromSmiles(smile)))
        d = dict([(key, value) for key, value in zip(name, alldescrs)])
        return d

    def pipeline_rf(self):

        df = pd.DataFrame([self.all_rdkit_descriptors(x) for x in self.data['SMILES']])
        X = df.copy()
        df['Experimental value'] = self.data['Experimental value']
        y = df['Experimental value']
        # numerical self.data
        numerical_transformer = SimpleImputer(strategy='median')
        # categorical self.data
        process = [('imputer', SimpleImputer(strategy='most_frequent')),
                ('onehot', OneHotEncoder(handle_unknown='ignore'))]
        categorical_transformer = Pipeline(steps=process)

        # find column cat and num
        categoricals_cols = [cname for cname in X.columns if
                            X[cname].nunique()<10 and
                            X[cname].dtype == 'object']

        numerical_cols = [cname for cname in X.columns if
                        X[cname].dtype in  ['int64', 'float64']]

        # numerical and categorical
        col_process = [('num', numerical_transformer, numerical_cols),
                    ('cat', categorical_transformer, categoricals_cols)]

        preprocessor = ColumnTransformer(transformers=col_process)
        model = RandomForestRegressor(n_estimators=500, max_depth=7, max_features=100, random_state=42)
        final_porcess = [('preprocessor', preprocessor), ('scaling', StandardScaler()), ('model', model)]
        my_pipeline = Pipeline(steps = final_porcess)

        my_pipeline.fit(X, y)

        # predict target
        df_target = pd.DataFrame([self.all_rdkit_descriptors(x) for x in self.target['SMILES']])
        X_target = df_target.copy()
        preds = my_pipeline.predict(X_target)
        return preds
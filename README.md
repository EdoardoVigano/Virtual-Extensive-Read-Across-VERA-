# Virtual-Extensive-Read-Across-VERA-

A new open-access software for chemical read-across using structural alerts and chemical molecular groups
Read-across applies the principle of similarity to identify the most similar substances to represent a given target substance in data poor situations. However, differences between the target and the source substances exist. The present study aims to screen and assess the effect of the key components in a molecule which may escape the evaluation for read-across based only on the most similar substance(s) using a new open-access software: Virtual Extensive Read-Across (VERA). VERA provides a means to assess similarity between chemicals using structural alerts specific to the property, pre-defined molecular groups and structural similarity. The software finds the most similar compounds with a certain feature e.g. structural alerts and molecular groups and provides clusters of similar substances while comparing these similar substances within different clusters. Carcinogenicity is a complex endpoint with several mechanisms, requiring resource intensive experimental bioassays, a large number of animals and the use of read-across as part of New Approach Methodologies would support carcinogenicity assessment. To test the VERA software, carcinogenicity was selected as the endpoint of interest for a range of botanicals. VERA correctly labelled 70% of the botanicals, indicating the most similar substances and the main features associated with carcinogenicity.

## INSTALLATION
## HOW TO INSTALL 
1. python and Anaconda are both required
2. install all packages in your env:
```python

        from rdkit import Chem
        from rdkit.Chem import Fragments
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

        # standard library
        import subprocess
        import os
        from tqdm import tqdm
        import shutil

        # library useful for output
        from fpdf import FPDF
        from PyPDF2 import PdfFileMerger
        import xlsxwriter
        from tqdm import tqdm
        import matplotlib.pyplot as plt
        from rdkit.Chem.Draw import rdMolDraw2D
        from fpdf import FPDF
        from PyPDF2 import PdfFileMerger
        import matplotlib.pyplot as plt
        import seaborn as sns
        from matplotlib.colors import ColorConverter 
        
```


3. install jupyter notebook in the same env
4. Download this repo 

# HOW TO USE THE VERA TOOL
1. open with notebook Vera.ipynb
2. choose the file names with target molecules.
  note: must be a .xlsx files with column name: SMILES
3. Run the code.
4. you can find a folder named **result** that contain pdf output and excel file
  

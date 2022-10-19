# Virtual-Extensive-Read-Across-VERA-

A new open-access software for chemical read-across using structural alerts and chemical molecular groups
Read-across applies the principle of similarity to identify the most similar substances to represent a given target substance in data poor situations. However, differences between the target and the source substances exist. The present study aims to screen and assess the effect of the key components in a molecule which may escape the evaluation for read-across based only on the most similar substance(s) using a new open-access software: Virtual Extensive Read-Across (VERA). VERA provides a means to assess similarity between chemicals using structural alerts specific to the property, pre-defined molecular groups and structural similarity. The software finds the most similar compounds with a certain feature e.g. structural alerts and molecular groups and provides clusters of similar substances while comparing these similar substances within different clusters. Carcinogenicity is a complex endpoint with several mechanisms, requiring resource intensive experimental bioassays, a large number of animals and the use of read-across as part of New Approach Methodologies would support carcinogenicity assessment. To test the VERA software, carcinogenicity was selected as the endpoint of interest for a range of botanicals. VERA correctly labelled 70% of the botanicals, indicating the most similar substances and the main features associated with carcinogenicity.

##### REFERENCE

Further details on the algorithms used to performe read across assessment can be found in the reference pubblication:

Viganò, E.L.; Colombo, E.; Raitano, G.; Manganaro, A.; Sommovigo, A.; Dorne, J.L.C.; Benfenati, E. Virtual Extensive Read-Across: A New Open-Access Software for Chemical Read-Across and Its Application to the Carcinogenicity Assessment of Botanicals. Molecules 2022, 27, 6605. https://doi.org/10.3390/molecules27196605 


## INSTALLATION
### HOW TO INSTALL 
1. _python_ and _Anaconda_ are both required
- Download and install _Anaconda_ from https://www.anaconda.com/distribution/. Choose Anaconda with Python 3
- Open the Anaconda Prompt (anaconda3).
- Create and activate a new RDKit Python environment (my-rdkit-env) with the following commands:
        
       > conda create -c https://conda.anaconda.org/rdkit -n my-rdkit-env rdkit
       > conda activate my-rdkit-env
        
2. Install _jupyter notebook_ in the same env: the procedure are described here https://jupyter.org/install
        
3. Install all _packages_ in your env: use the following command to install python libraries required for VERA. 
        
        pip install package_name
        #pip install pandas

You can find all packages required in the first chunk of VERA.ipynb file:
        
```python
# Data Analisys library
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
import matplotlib.pyplot as plt
from PyPDF2 import PdfFileMerger
import seaborn as sns
from matplotlib.colors import ColorConverter 
        
```
4. Download this repository 

## HOW TO USE THE VERA TOOL
1. Open with notebook Vera.ipynb
2. Choose the file names with target molecules.
  note: must be a .xlsx files with column name: SMILES
3. Run the code
4. You can find a folder named **result** that contain pdf output and excel file

## CONTACT

Edoardo Luca Viganò

Laboratory of Environmental Chemistry and Toxicology

Department of Environmental Health Sciences

Istituto di Ricerche Farmacologiche Mario Negri IRCCS

Via Mario Negri 2, 20156 Milano, Italy

e-mail: edoardo.vigano@marionegri.it


# Calculations
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from scipy.stats import mannwhitneyu
#from padelpy import padeldescriptor # https://github.com/ecrl/padelpy

# Statistics
import numpy as np # necessary?
import pandas as pd

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns

# ML model
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold

# File management
import os
import io
import base64
    
def preprocess(df):

    # Remove NaN values from standard_value
    df = df.dropna(subset=['standard_value'])

    # Simply to contain ChEMBL IDs, SMILES, standard val
    df = df[['molecule_chembl_id', 'canonical_smiles', 'standard_value']]

    bioactivity_class = pd.Series(name='bioactivity_class')

    # Iterate through standard_value
    for i, val in df['standard_value'].items():

        # Set class based on condition
        if float(val) >= 10000:
            bioactivity_class.loc[i] = 'inactive' 
        elif float(val) <= 1000:
            bioactivity_class.loc[i] = 'active'
        else:
            bioactivity_class.loc[i] = 'intermediate'

    # Append bioactivity classification
    bioactivity = pd.Series(bioactivity_class, name='bioactivity_class')
    df = pd.concat([df, bioactivity], axis=1)

    return df

def lipinski(smiles):
    # Compute Lipinski descriptors for druglikeness

    # Convert SMILES to RDKit Mol
    moldata= []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile) 
        moldata.append(mol)
       
    # Compute Lipinski descriptors
    descriptors = []
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol) # Molecular weight
        desc_MolLogP = Descriptors.MolLogP(mol) # log of octanol:water coefficient (P)
        desc_NumHDonors = Lipinski.NumHDonors(mol) # Total H bond donors
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol) # Total H bond acceptors

        row = [desc_MolWt, desc_MolLogP, desc_NumHDonors, desc_NumHAcceptors]
        descriptors.append(row)

    columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]

    # Build df for descriptor data
    results_df = pd.DataFrame(data=descriptors, columns=columnNames)

    return results_df

def pIC50(df, max):
    # Cap IC50 value & Convert IC50 to pIC50

    pIC50 = []
    for i in df['standard_value']:
        if float(i) > max:
            i = max
        i = float(i) * (1e-9)
        pIC50.append(-np.log10(i)) # -log10(M) -> could change to import math, -math.log(i, 10)

    df['pIC50'] = pIC50
    df = df.drop('standard_value', axis=1)
        
    return df

def mannwhitney(df, descriptor):
    # Mann-Whitney U Test
    # https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python/

    # Extract data for active and inactive compounds
    active = df[df.bioactivity_class == 'active'][descriptor]
    inactive = df[df.bioactivity_class == 'inactive'][descriptor]

    # Compare active and inactive samples
    stat, p = mannwhitneyu(active, inactive)
    print(stat, p)

    # Interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'

    results = {'Descriptor':descriptor, 'Statistics':stat, 'p':p, 'alpha':alpha, 'Interpretation':interpretation}
    
    # Build df for Mann-Whitney U data (only necessary for saving data)
    results_df = pd.DataFrame({
        'Descriptor':descriptor,
        'Statistics':stat,
        'p':p,
        'alpha':alpha,
        'Interpretation':interpretation
        },index=[0])

    return results

def plotToImage(fig):
    
    # Converts plt plot to .png to display on webapp
    
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png')
    buffer.seek(0)
    plt.close(fig)
    
    img_base64 = base64.b64encode(buffer.read()).decode()

    return img_base64

def padeldescriptors(df):

    # Converts SMILES to PubChem Fingerprints w/ padel software
    
    df.to_csv('molecule.smi', sep='\t', index=False, header=False)

    try:
        exit_code = os.system('bash padel.sh') 
        if exit_code != 0:
            print(f"Error: Script exited with non-zero code: {exit_code}")
    except Exception as e:
        print(f"Error: {e}")

    df_X = pd.read_csv('descriptors_output.csv').drop(columns=['Name'])
    df_Y = df['pIC50']

    df_out = pd.concat([df_X,df_Y], axis=1)

    # TODO: Delete both csvs that get created

    return df_out

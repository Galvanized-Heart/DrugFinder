
# Flask app
from flask import Flask, redirect, render_template, request, session
from flask_session import Session

# ChEMBL REST API
from chembl_webresource_client.new_client import new_client # https://github.com/chembl/chembl_webresource_client/

# Statistics
import numpy as np
import pandas as pd

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns

# Helper functions
from helpers import *



# Configure app
app = Flask(__name__)

# Configure session
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_TYPE"] = "filesystem"
Session(app)

# Configure ChEMBL queries
client = new_client.target
activity = new_client.activity

# Configure matplotlib to Agg backend
plt.switch_backend('Agg')



@app.route('/', methods=['GET'])
def index():

    session.clear()

    return render_template('index.html')

@app.route('/targets', methods=['POST'])
def targets():
        
    target = request.form.get('target')

    if len(target) < 3:
        target_query = [{'target_chembl_id': 'Query Too Short'}]

    else:
        # Query ChEMBL for targets with user search
        target_query = client.search(target)

    return render_template('targets.html', target_query=target_query)

@app.route('/drugs', methods=['POST'])
def drugs():

    selected_target = request.form.get('selected_target')

    # Query ChEMBL for activity with user target selection
    activity_query = activity.filter(target_chembl_id=selected_target).filter(standard_type='IC50') # Query may return [] if there are not bioactivities!!!

    if not activity_query:
        drug_query = [{'molecule_chembl_id': 'No Data'}]
        stats_buttons = False
    else:
        # Create data frames
        df_raw = pd.DataFrame.from_dict(activity_query)
        df_pre = preprocess(df_raw)

        # Set current df
        session['curr_df'] = df_pre
        drug_query = df_pre.to_dict('records')
        stats_buttons = True

    return render_template('drugs.html', drug_query=drug_query, stats_buttons=stats_buttons)

@app.route('/descriptors', methods=['POST'])
def descriptors():

    # Fetch current df
    df = session['curr_df']

    # Add Lipinski descriptors to df
    df_lipinski = lipinski(df.canonical_smiles)
    df = pd.concat([df, df_lipinski], axis=1)
    
    # Convert IC50 to pIC50 (cap at IC50 = 100,000,000 nM)
    df = pIC50(df, max=100_000_000.0)

    # Remove values that have intermediate from df
    df = df[df.bioactivity_class != 'intermediate'] # must be at end for data mapping purposes

    images = []
    # ** PLOT 1: frequency of active and inactive compounds countplot**
    fig = plt.figure(figsize=(5.5, 5.5))
    sns.countplot(x='bioactivity_class', data=df, edgecolor='black')
    plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
    plt.ylabel('Frequency', fontsize=14, fontweight='bold')
    plt.title('Frequency of Bioactivity Classes', fontsize=16, fontweight='bold')
    images.append(plotToImage(fig))

    # ** PLOT 2: distribution of activity classes MW vs logP scatterplot**
    fig = plt.figure(figsize=(10, 5.5))
    sns.scatterplot(x='MW', y='LogP', data=df, hue='bioactivity_class', size='pIC50', edgecolor='black', alpha=0.7)
    plt.xlabel('MW', fontsize=14, fontweight='bold')
    plt.ylabel('LogP', fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.title('Bioactivity Distribution by MW by logP', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    images.append(plotToImage(fig))

    # ** PLOT 3-7: pIC50, MW, logP, NumHDon, NumHAcc by activity class boxplots **
    utest = []
    descriptors = ['pIC50', 'MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
    for descriptor in descriptors:

        fig = plt.figure(figsize=(5.5, 5.5))
        sns.boxplot(x='bioactivity_class', y=descriptor, data=df)
        plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
        plt.ylabel(descriptor, fontsize=14, fontweight='bold')
        plt.title(f'Class vs. {descriptor}', fontsize=16, fontweight='bold')
        images.append(plotToImage(fig))

        utest.append(mannwhitney(df, descriptor))
        
    
    

    return render_template('descriptors.html', images=images, utest=utest)

@app.route('/models', methods=['POST'])
def models():
    
    # Fetch current df
    df = session['curr_df']

    # Convert IC50 to pIC50 (cap at IC50 = 100,000,000 nM)
    df = pIC50(df, max=100_000_000.0)

    df.to_csv('molecule.smi', sep='\t', index=False, header=False)

if __name__ == '__main__':
    app.run(debug=True)

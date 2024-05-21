import requests
import pandas as pd
import numpy as np
import pickle


def fetch_ensembl_ids(gene_names):
    '''
    Takes a list of gene names and returns a dictionary where the keys are gene names and the values are ENSEMBL gene IDs.
    
    Parameters:
    - 
    
    Returns:
    - 
    
    '''
    # Need to make sure getting the correct ENSEMBL gene version!
    #server = "https://rest.ensembl.org"
    server = "https://grch37.rest.ensembl.org"  # Example for using the GRCh37 assembly
    ext = "/lookup/symbol/homo_sapiens/{}?expand=0"
    headers = {"Content-Type": "application/json"}
    
    gene_id_map = {}
    
    for gene in gene_names:
        r = requests.get(f"{server}{ext.format(gene)}", headers=headers)
        if r.status_code == 200:
            decoded = r.json()
            gene_id = decoded.get('id')
            if gene_id:
                gene_id_map[gene] = gene_id
            else:
                gene_id_map[gene] = 'ID not found'
        else:
            gene_id_map[gene] = 'Error fetching ID'
    
    return gene_id_map


def fetch_gene_names_from_ids(ensembl_ids):
    '''
    Takes a list of ENSEMBL gene IDs and returns a dictionary with ENSEMBL IDs as keys and gene names as values.
    
    Parameters:
    - ensembl_ids: A list of ENSEMBL gene IDs.
    
    Returns:
    - gene_name_map: A dictionary housing the ENSEMBL IDs as keys and gene names as values.
    
    '''
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/{}"
    headers = {'Content-Type': 'application/json'}
    
    gene_name_map = {}
    
    for ensembl_id in ensembl_ids:
        response = requests.get(f"{server}{ext.format(ensembl_id)}", headers=headers)
        if response.status_code == 200:
            data = response.json()
            gene_name = data.get('display_name')
            if gene_name:
                gene_name_map[ensembl_id] = gene_name
            else:
                gene_name_map[ensembl_id] = 'Gene name not found'
        else:
            gene_name_map[ensembl_id] = 'Error fetching gene name'
    
    return gene_name_map


def sample_rna_seq_data(normalized_scaled_bulk_rna, used_sample_ids, num_samples, random_seed=None):
    """
    Randomly sample rows from normalized_scaled_bulk_rna avoiding the samples already used.

    Parameters:
    - normalized_scaled_bulk_rna: DataFrame to sample from.
    - used_sample_ids: List of sample IDs to avoid.
    - num_samples: Number of samples to draw.
    - random_seed: Random seed for reproducibility (default: None).

    Returns:
    - DataFrame with the sampled rows.
    
    """
    # Set the random seed for reproducibility
    if random_seed is not None:
        np.random.seed(random_seed)
    
    # Get available sample IDs by excluding used sample IDs
    available_sample_ids = normalized_scaled_bulk_rna.index.difference(used_sample_ids)
    
    # If there are fewer available samples than requested, adjust the number of samples
    if len(available_sample_ids) < num_samples:
        print(f"Only {len(available_sample_ids)} available samples, adjusting the number of samples to {len(available_sample_ids)}")
        num_samples = len(available_sample_ids)
    
    # Randomly sample from the available sample IDs
    sampled_ids = np.random.choice(available_sample_ids, size=num_samples, replace=False)
    
    # Filter the DataFrame to include only the sampled IDs
    sampled_rna_seq_data = normalized_scaled_bulk_rna.loc[sampled_ids]
    
    return sampled_rna_seq_data


# Function to load a model from a pickle file
def load_model(directory_name, gene_name):
    model_filename = directory_name + f"{gene_name}_logistic_ridge_model.pkl"
    with open(model_filename, 'rb') as f:
        model = pickle.load(f)
    return model


# Functions to extract the necessary information from gene dictionaries
def get_rna_seq_data(gene_data, gene_name):
    gene_info = gene_data[gene_name]
    mut_data = gene_info['mut_sequencing_data']
    non_mut_data = gene_info['non-mut_sequencing_data']
    
    X = pd.concat([mut_data, non_mut_data], ignore_index=True)
    y_mut = np.ones(mut_data.shape[0])
    y_non_mut = np.zeros(non_mut_data.shape[0])
    y = np.concatenate([y_mut, y_non_mut])
    
    return X, y



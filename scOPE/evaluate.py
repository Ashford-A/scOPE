import pandas as pd
import numpy as np
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    confusion_matrix,
    roc_curve,
    precision_recall_curve,
    average_precision_score,
    classification_report
)
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

# Import custom modules
sys.path.append('../')
from scOPE import preprocessing
from scOPE import utilities
from scOPE import models
from scOPE import evaluate


def predict_single_cell_mutation(single_cell_data, model):
    '''Predicts the mutation presence or absence in single-cell data using the provided model.'''
    # Ensure the single-cell data is in the correct shape (transpose if necessary)
    #if single_cell_data.shape[0] > single_cell_data.shape[1]:
        #single_cell_data = single_cell_data.T
    # Use the model to predict the mutation status
    predictions = model.predict(single_cell_data)
    return predictions


def predict_mutations_in_single_cells_logistic(gene_models_dir, single_cell_data):
    '''
    Applies each gene mutation prediction model to the single-cell RNA-seq data.
    
    Parameters:
    - gene_models_dir: Directory where the trained models are saved.
    - single_cell_data: DataFrame containing the single-cell RNA-seq data.
    
    Returns:
    - A dictionary with gene names as keys and prediction results as values.
    '''
    predictions_dict = {}
    
    # Iterate through the files in the gene models directory
    for model_file in os.listdir(gene_models_dir):
        if model_file.endswith('_logistic_ridge_model.pkl'):
            
            # Extract gene name from the file name
            gene_name = model_file.split('_logistic_ridge_model.pkl')[0]
            
            # Load the model
            model_path = os.path.join(gene_models_dir, model_file)
            model = utilities.load_model(gene_models_dir, gene_name)
            
            # Predict mutations in single-cell data
            predictions = predict_single_cell_mutation(single_cell_data, model)
            
            # Store predictions in the dictionary
            predictions_dict[gene_name] = predictions
            
            print(f"Predictions for {gene_name} completed.")
    
    return predictions_dict


def predict_mutations_in_single_cells_random_forest(gene_models_dir, single_cell_data):
    '''
    Applies each gene mutation prediction model to the single-cell RNA-seq data.
    
    Parameters:
    - gene_models_dir: Directory where the trained models are saved.
    - single_cell_data: DataFrame containing the single-cell RNA-seq data.
    
    Returns:
    - A dictionary with gene names as keys and prediction results as values.
    '''
    predictions_dict = {}
    
    # Iterate through the files in the gene models directory
    for model_file in os.listdir(gene_models_dir):
        if model_file.endswith('_random_forest_model.pkl'):
            
            # Extract gene name from the file name
            gene_name = model_file.split('_random_forest_model.pkl')[0]
            
            # Load the model
            model_path = os.path.join(gene_models_dir, model_file)
            model = utilities.load_model(gene_models_dir, gene_name)
            
            # Predict mutations in single-cell data
            predictions = predict_single_cell_mutation(single_cell_data, model)
            
            # Store predictions in the dictionary
            predictions_dict[gene_name] = predictions
            
            print(f"Predictions for {gene_name} completed.")
    
    return predictions_dict


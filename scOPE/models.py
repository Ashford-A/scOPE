import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
import pickle



def logistic_ridge_regression(gene_data, output_dir, test_size=0.2):
    '''
    Takes in a dictionary of gene data (an example of which is shown below in Parameters) and trains logistic ridge 
    regression models to predict the presence or absence of a mutation within the specific gene using the samples'
    corresponding RNA-seq data.
    
    Parameters:
    - gene_data: This is a dictionary where the genes of interest are the keys, and they hold necessary information and 
                 RNA-seq data to train the classification model. An example format of the dictionary is shown below:
                 
                 # Sample data structure
                 gene_data = {
                    'FLT3': {
                        'gene_id': 'gene_1',
                        'positive_sample_count': 10,
                        'mut_sequencing_data': pd.DataFrame(np.random.rand(10, 100)),
                        'non-mut_sequencing_data': pd.DataFrame(np.random.rand(10, 100))
                    },
                    'RUNX1': {
                        'gene_id': 'gene_2',
                        'positive_sample_count': 15,
                        'mut_sequencing_data': pd.DataFrame(np.random.rand(15, 100)),
                        'non-mut_sequencing_data': pd.DataFrame(np.random.rand(15, 100))
                    },
                    # Add other genes similarly...
                 }
                 
    - output_dir: A directory string in the following format: '/home/groups/users/ashforda/chosen_output_dir/'
                  # NOTE: The slash at the end is important!
    
    - test_size: A float representing the proportion of the dataset to include in the test split. Default is 0.2.
    
    Returns:
    - The function saves a trained model to the directory specified as output_dir. Function does not "return" anything.
    
    '''

    # Iterate through the dictionary
    for gene, data in gene_data.items():
        # Get the sequencing data for mutated and non-mutated samples
        mut_data = data['mut_sequencing_data']
        non_mut_data = data['non-mut_sequencing_data']
    
        # Create labels (1 for mutated, 0 for non-mutated)
        mut_labels = np.ones(mut_data.shape[0])
        non_mut_labels = np.zeros(non_mut_data.shape[0])
    
        # Combine the data and labels
        X = pd.concat([mut_data, non_mut_data], ignore_index=True)
        y = np.concatenate([mut_labels, non_mut_labels])
    
        # Split the data into training and test sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
    
        # Train the logistic ridge regression model
        model = make_pipeline(StandardScaler(), LogisticRegression(penalty='l2', solver='liblinear'))
        model.fit(X_train, y_train)
    
        # Save the trained model to a pickle file
        model_filename = f"{gene}_logistic_ridge_model.pkl"
        with open(output_dir + model_filename, 'wb') as f:
            pickle.dump(model, f)
    
        print(f"Trained and saved model for gene: {gene} as {model_filename}")
        print(f"Training set size: {X_train.shape[0]}, Test set size: {X_test.shape[0]}")
    
        # Load the saved model
        with open(output_dir + model_filename, 'rb') as f:
            loaded_model = pickle.load(f)
    
        # Evaluate the model using the test set
        y_pred = loaded_model.predict(X_test)
    
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred)
        conf_matrix = confusion_matrix(y_test, y_pred)
    
        print(f"Evaluation for gene: {gene}")
        print(f"Accuracy: {accuracy}")
        print(f"Precision: {precision}")
        print(f"Recall: {recall}")
        print(f"F1 Score: {f1}")
        print(f"Confusion Matrix:\n{conf_matrix}")



def random_forest_classification(gene_data, output_dir, test_size=0.2):
    '''
    Takes in a dictionary of gene data (an example of which is shown below in Parameters) and trains random forest 
    classification models to predict the presence or absence of a mutation within the specific gene using the samples'
    corresponding RNA-seq data.
    
    Parameters:
    - gene_data: This is a dictionary where the genes of interest are the keys, and they hold necessary information and 
                 RNA-seq data to train the classification model. An example format of the dictionary is shown below:
                 
                 # Sample data structure
                 gene_data = {
                    'FLT3': {
                        'gene_id': 'gene_1',
                        'positive_sample_count': 10,
                        'mut_sequencing_data': pd.DataFrame(np.random.rand(10, 100)),
                        'non-mut_sequencing_data': pd.DataFrame(np.random.rand(10, 100))
                    },
                    'RUNX1': {
                        'gene_id': 'gene_2',
                        'positive_sample_count': 15,
                        'mut_sequencing_data': pd.DataFrame(np.random.rand(15, 100)),
                        'non-mut_sequencing_data': pd.DataFrame(np.random.rand(15, 100))
                    },
                    # Add other genes similarly...
                 }
                 
    - output_dir: A directory string in the following format: '/home/groups/users/ashforda/chosen_output_dir/'
                  # NOTE: The slash at the end is important!
    
    - test_size: A float representing the proportion of the dataset to include in the test split. Default is 0.2.
    
    Returns:
    - The function saves a trained model to the directory specified as output_dir. Function does not "return" anything.
    
    '''

    # Iterate through the dictionary
    for gene, data in gene_data.items():
        # Get the sequencing data for mutated and non-mutated samples
        mut_data = data['mut_sequencing_data']
        non_mut_data = data['non-mut_sequencing_data']
    
        # Create labels (1 for mutated, 0 for non-mutated)
        mut_labels = np.ones(mut_data.shape[0])
        non_mut_labels = np.zeros(non_mut_data.shape[0])
    
        # Combine the data and labels
        X = pd.concat([mut_data, non_mut_data], ignore_index=True)
        y = np.concatenate([mut_labels, non_mut_labels])
    
        # Split the data into training and test sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
    
        # Train the random forest model
        model = make_pipeline(StandardScaler(), RandomForestClassifier(n_estimators=100, random_state=42))
        model.fit(X_train, y_train)
    
        # Save the trained model to a pickle file
        model_filename = f"{gene}_random_forest_model.pkl"
        with open(output_dir + model_filename, 'wb') as f:
            pickle.dump(model, f)
    
        print(f"Trained and saved model for gene: {gene} as {model_filename}")
        print(f"Training set size: {X_train.shape[0]}, Test set size: {X_test.shape[0]}")
    
        # Load the saved model
        with open(output_dir + model_filename, 'rb') as f:
            loaded_model = pickle.load(f)
    
        # Evaluate the model using the test set
        y_pred = loaded_model.predict(X_test)
    
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred)
        conf_matrix = confusion_matrix(y_test, y_pred)
    
        print(f"Evaluation for gene: {gene}")
        print(f"Accuracy: {accuracy}")
        print(f"Precision: {precision}")
        print(f"Recall: {recall}")
        print(f"F1 Score: {f1}")
        print(f"Confusion Matrix:\n{conf_matrix}")
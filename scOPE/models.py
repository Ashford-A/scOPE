import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
import pickle



def logistic_ridge_regression(gene_data, output_dir, test_size=0.2):
    '''
    Takes in a dictionary of gene data and trains logistic ridge regression models to predict the presence or absence 
    of a mutation within the specific gene using the samples' corresponding RNA-seq data.
    '''
    # Define the parameter grid for logistic ridge regression
    param_grid = {
        'logisticregression__C': [0.1, 1.0, 10.0, 100.0],
        'logisticregression__solver': ['liblinear']
    }

    # Iterate through the dictionary
    for gene, data in gene_data.items():
        mut_data = data['mut_sequencing_data']
        non_mut_data = data['non-mut_sequencing_data']
        mut_labels = np.ones(mut_data.shape[0])
        non_mut_labels = np.zeros(non_mut_data.shape[0])
        X = pd.concat([mut_data, non_mut_data], ignore_index=True)
        y = np.concatenate([mut_labels, non_mut_labels])

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)

        model = make_pipeline(StandardScaler(), LogisticRegression(penalty='l2'))
        
        grid_search = GridSearchCV(model, param_grid, cv=4, scoring='accuracy', n_jobs=-1, refit=True)
        
        # Perform the grid search with cross-validation
        grid_search.fit(X_train, y_train)

        # Get the best model
        best_model = grid_search.best_estimator_

        # Save the trained model to a pickle file
        model_filename = f"{gene}_logistic_ridge_model.pkl"
        with open(output_dir + model_filename, 'wb') as f:
            pickle.dump(best_model, f)

        print(f"Trained and saved model for gene: {gene} as {model_filename}")
        print(f"Training set size: {X_train.shape[0]}, Test set size: {X_test.shape[0]}")

        # Evaluate the model using the test set
        y_pred = best_model.predict(X_test)

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
    Takes in a dictionary of gene data and trains random forest classification models to predict the presence or absence 
    of a mutation within the specific gene using the samples' corresponding RNA-seq data.
    '''
    # Define the parameter grid for random forest
    param_grid = {
        'randomforestclassifier__n_estimators': [100, 200, 300],
        'randomforestclassifier__max_depth': [10, 20, 30],
        'randomforestclassifier__min_samples_split': [2, 5, 10],
        'randomforestclassifier__min_samples_leaf': [1, 2, 4],
        'randomforestclassifier__bootstrap': [True, False]
    }

    # Iterate through the dictionary
    for gene, data in gene_data.items():
        mut_data = data['mut_sequencing_data']
        non_mut_data = data['non-mut_sequencing_data']
        mut_labels = np.ones(mut_data.shape[0])
        non_mut_labels = np.zeros(non_mut_data.shape[0])
        X = pd.concat([mut_data, non_mut_data], ignore_index=True)
        y = np.concatenate([mut_labels, non_mut_labels])

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)

        model = make_pipeline(StandardScaler(), RandomForestClassifier(random_state=42))

        grid_search = GridSearchCV(model, param_grid, cv=4, scoring='accuracy', n_jobs=-1, refit=True)
        
        # Perform the grid search with cross-validation
        grid_search.fit(X_train, y_train)

        # Get the best model
        best_model = grid_search.best_estimator_

        # Save the trained model to a pickle file
        model_filename = f"{gene}_random_forest_model.pkl"
        with open(output_dir + model_filename, 'wb') as f:
            pickle.dump(best_model, f)

        print(f"Trained and saved model for gene: {gene} as {model_filename}")
        print(f"Training set size: {X_train.shape[0]}, Test set size: {X_test.shape[0]}")

        # Evaluate the model using the test set
        y_pred = best_model.predict(X_test)

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


        
        
       
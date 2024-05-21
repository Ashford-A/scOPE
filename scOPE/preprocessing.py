# Import modules
import pandas as pd
from sklearn.preprocessing import Normalizer, StandardScaler


# Define functions
def tsv_to_df(file_path, delimiter='\t', index_col=None):
    '''
    Function that takes in a .tsv file and stores it in a pandas dataframe
    
    Delimiter defaults to '\t'
    index_col defaults to none, in the case of wanting to name the rows after the first column, use index_col=0
    otherwise, row names default to the row count index
    
    Returns Pandas dataframe of tsv file contents
    '''
    
    dataframe = pd.read_csv(file_path, delimiter=delimiter, index_col=index_col)
    
    return dataframe
    
    
def preprocess_bulk_RNA(bulk_transcript_df, normalize=True, scale=True):
    '''
    
    '''
    if normalize == True:
        print('Normalizing data..')
        # Normalize data
        normalizer = Normalizer()
        df_normalized = normalizer.fit_transform(bulk_transcript_df)

    if scale == True:
        print('Scaling data..')
        # Scale data
        scaler = StandardScaler()
        df_scaled = scaler.fit_transform(bulk_transcript_df)

    # Optionally convert back to DataFrame
    df_scaled = pd.DataFrame(df_scaled, index=df_transposed.index, columns=df_transposed.columns)
    
    return df_scaled

 


    
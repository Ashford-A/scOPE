# Import modules
import pandas as pd
from sklearn.preprocessing import Normalizer, StandardScaler
import scanpy as sc


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
    if normalize == False and scale == False:
        print('Function aborting, all options set to False.')
        return None
    
    out_df = bulk_transcript_df
    #print(out_df)
    
    if normalize == True:
        print('Normalizing data..')
        # Normalize data
        normalizer = Normalizer()
        out_df = normalizer.fit_transform(out_df)

    if scale == True:
        print('Scaling data..')
        # Scale data
        scaler = StandardScaler()
        out_df = scaler.fit_transform(out_df)

    # Optionally convert back to DataFrame
    preprocessed_df = pd.DataFrame(out_df, index=bulk_transcript_df.index, columns=bulk_transcript_df.columns)
    
    return preprocessed_df


def preprocess_single_cell_RNA(adata, normalize=True, scale=True):
    '''
    
    '''
    if normalize == False and scale == False:
        print('Function aborting, all options set to False.')
        return None
    
    # Assuming 'adata' is your AnnData object loaded with scRNA-seq data
    if normalize == True:
        print('Normalizing total counts..')
        sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize total counts to 10,000 per cell
        
        print('Logarithmizing the data..')
        sc.pp.log1p(adata)  # Logarithmize the data
    
    if scale == True:
        print('Scaling data..')
        sc.pp.scale(adata)  # Scale to zero mean and unit variance
        
    return adata







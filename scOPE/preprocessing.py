# Import modules
import pandas as pd


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
    
    

 


    
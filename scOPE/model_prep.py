import requests

def fetch_ensembl_ids(gene_names):
    server = "https://rest.ensembl.org"
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



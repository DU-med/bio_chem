# 遺伝子情報のキーワード検索
import pandas as pd

df = pd.read_csv('../data/ncbi/homo_gene.txt', sep='\t')
df['Aliases'] = df['Aliases'].fillna('') 
id_df = df.iloc[:,[2,5,7]]
s_df = df.iloc[:,[2,5,6,7]]

# 完全に遺伝子名が一致する遺伝子情報の表示(簡易)
def gene_id(gene):
    return id_df[id_df['Symbol']==gene]

# 完全に遺伝子名が一致する遺伝子情報の表示(詳細)
def gene_id_full(gene):
    return df[df['Symbol']==gene]

# 完全に遺伝子名が一部一致する遺伝子情報の表示(簡易)
def gene_search(gene):
    return s_df[(s_df['Symbol']==gene) | (s_df['Aliases'].str.contains(gene))]

# 完全に遺伝子名が一部一致する遺伝子情報の表示(詳細)
def gene_search_full(gene):
    return df[(df['Symbol']==gene) | (df['Aliases'].str.contains(gene))]
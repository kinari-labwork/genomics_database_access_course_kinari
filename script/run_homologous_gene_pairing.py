import pandas as pd

pd.set_option('display.max_colwidth', None)  
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)


mice_annotated_exons = pd.read_pickle('./data/output/mm39_annotated_exons.pickle')
human_annotated_exons = pd.read_pickle('./data/output/hg38_annotated_exons.pickle')

homologous_df = pd.read_csv('data/homologous_pair/homologous_pair.txt', sep=',', header = 0)
homologous_df.columns = homologous_df.columns.str.strip()
print(homologous_df.columns)

print(len(homologous_df)) #92566
print(homologous_df["Gene stable ID"].duplicated().sum()) #6232
print(homologous_df["Mouse homology type"].value_counts())

"""
Mouse homology type
ortholog_one2one      17178
ortholog_many2many     5426
ortholog_one2many      2987
"""

print(mice_annotated_exons.head())
print(mice_annotated_exons.columns)
print(human_annotated_exons.head())
print(human_annotated_exons.columns)
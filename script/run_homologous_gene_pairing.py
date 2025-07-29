import pandas as pd
from src.homologous_gene_pairing import extract_genename_as_column, make_homologous_gene_tuple, compare_and_pair_exons
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

print(mice_annotated_exons.columns)
print(human_annotated_exons.columns)

mice_annotated_exons = extract_genename_as_column(mice_annotated_exons)
human_annotated_exons = extract_genename_as_column(human_annotated_exons)

homologous_df = make_homologous_gene_tuple(homologous_df)

paired_exons = compare_and_pair_exons(homologous_df, mice_annotated_exons.groupby('gene_id'), human_annotated_exons.groupby('gene_id'))

# 結果を保存
paired_exons_df = pd.DataFrame(paired_exons, columns=['human_gene_id', 'mouse_gene_id', 'human_exon_id', 'mouse_exon_id', 'homology_percent', 'mouse_homology_type', 'overall_homology_percent'])
paired_exons_df.to_pickle("./data/output/paired_exons.pkl")

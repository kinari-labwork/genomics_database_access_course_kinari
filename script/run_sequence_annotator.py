import pandas as pd
from src.sequence_anotator import extract_exons, annotate_exons_with_sequences

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)

mice_gtf_file = pd.read_csv('data/mm39/mm39.gtf', sep='\t', comment='#', header=None, dtype={0: str})
mice_fasta_file = 'data/mm39/mm39.fa'
human_gtf_file = pd.read_csv('data/hg38/hg38.gtf', sep='\t', comment='#', header=None, dtype={0: str})
human_fasta_file = 'data/hg38/hg38.fa'

# エキソン情報を抽出
mice_exon_data = extract_exons(mice_gtf_file)
human_exon_data = extract_exons(human_gtf_file)

# エキソンにゲノム配列をアノテーション
mice_annotated_exons = annotate_exons_with_sequences(mice_exon_data, mice_fasta_file)
human_annotated_exons = annotate_exons_with_sequences(human_exon_data, human_fasta_file)

# 結果を保存
mice_annotated_exons.to_pickle('./data/output/mm39_annotated_exons.pickle')
human_annotated_exons.to_pickle('./data/output/hg38_annotated_exons.pickle')

# 結果を確認
print(f"Mice Annotated Exons:{mice_annotated_exons.head()}")
print(f"Number of Mice Annotated Exons: {len(mice_annotated_exons)}") #689517


print(f"Human Annotated Exons:{human_annotated_exons.head()}")
print(f"Number of Human Annotated Exons: {len(human_annotated_exons)}")
 
print(f"Number of unannotated exons in Mice: {mice_annotated_exons.isna().sum()}") # 0
print(f"Number of unannotated exons in Human: {human_annotated_exons.isna().sum()}")

print(mice_annotated_exons.columns)
import pandas as pd
from src.sequence_anotator import extract_exons, annotate_exons_with_sequences

mice_gtf_file = pd.read_csv('data/mm39/mm39.gtf', sep='\t', comment='#', header=None)
mice_fasta_file = 'data/mm39/mm39.fa'
human_gtf_file = pd.read_csv('data/hg38/hg38.gtf', sep='\t', comment='#', header=None)
human_fasta_file = 'data/hg38/hg38.fa'

# エキソン情報を抽出
mice_exon_data = extract_exons(mice_gtf_file)
human_exon_data = extract_exons(human_gtf_file)

# エキソンにゲノム配列をアノテーション
mice_annotated_exons = annotate_exons_with_sequences(mice_exon_data, mice_fasta_file)
human_annotated_exons = annotate_exons_with_sequences(human_exon_data, human_fasta_file)
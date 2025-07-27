import pandas as pd

def extract_genename_as_column(annotated_exons):
    """
    Pourpose: attribute列に複数の情報が格納されているので、gene nameを抽出して新しい列に格納する
    input: annotated_exons (DataFrame): アノテーションされたエキソンのデータ
    output: annotated_exons (DataFrame): gene nameを含む新しい列を追加したDataFrame
    """
    # Extract gene names from the 9th column (attributes)
    annotated_exons['gene_id'] = annotated_exons[8].str.extract(r'gene_id "([^"]+)"')
    annotated_exons['gene_name'] = annotated_exons[8].str.extract(r'gene_name "([^"]+)"')
    annotated_exons['transcript_id'] = annotated_exons[8].str.extract(r'transcript_id "([^"]+)"')
    annotated_exons['transcript_biotype'] = annotated_exons[8].str.extract(r'transcript_biotype "([^"]+)"')
    annotated_exons['exon_id'] = annotated_exons[8].str.extract(r'exon_id "([^"]+)"')
    annotated_exons['exon_number'] = annotated_exons[8].str.extract(r'exon_number "([^"]+)"')
    annotated_exons.drop(columns=[8], inplace=True)
    return annotated_exons

def make_homologous_gene_dict(homologous_df):
    """
    Pourpose: homologous_dfからhg38のgene_idをキー、mm39のgene_idをvalueとする辞書を作成する
    """
import pandas as pd
import numpy as np
from Bio.Align import PairwiseAligner
from multiprocessing import Pool
from tqdm import tqdm  # 進捗バー用ライブラリ


def extract_genename_as_column(annotated_exons) -> pd.DataFrame:
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

def make_homologous_gene_tuple(homologous_df) -> tuple:
    """
    pourpose: 後で使用するために、homologous_dfからマウスとヒトの遺伝子のペアをタプルとして抽出する
    input: homologous_df (DataFrame): homologous pairのデータ
    output: homologous_tuple (tuple): ヒトとマウスの遺伝子のペアを含むタプルのリスト
    """
    homologous_df = homologous_df.dropna(subset="Mouse gene stable ID")
    homologous_df = homologous_df[homologous_df["Mouse homology type"].isin(["ortholog_one2one", "ortholog_one2many"])]
    homologous_tuple = tuple(zip(homologous_df["Gene stable ID"], 
                                 homologous_df["Mouse gene stable ID"], 
                                 homologous_df["Mouse homology type"], 
                                 homologous_df["%id. query gene identical to target Mouse gene"]))
    return homologous_tuple


def calculate_homology_percent(seq1, seq2) -> float:
    """Biopythonを使って2つの配列の相同性(%)を計算する"""
    if not isinstance(seq1, str) or not isinstance(seq2, str) or not seq1 or not seq2:
        return 0.0 
    aligner = PairwiseAligner()
    aligner.mode = 'global'  # グローバルアライメントを使用
    alignments = aligner.align(seq1, seq2)
    alignment = alignments[0]
    difference = np.array(alignment.aligned)[0][:, 1] - np.array(alignment.aligned)[0][:, 0]
    identity = np.sum(difference)
    return (identity / len(seq1)) * 100 if len(seq1) > 0 else 0.0 # scoreをmatchに変えて実行することもできる

def calculate_homology_for_exon_pair(args):
    h_exon, mouse_exons_df, human_gene_id, mouse_gene_id, mouse_homology_type, overall_homology_percent = args
    best_score = -1
    best_mouse_exon = None

    for _, m_exon in mouse_exons_df.iterrows():
        homology_score = calculate_homology_percent(h_exon['sequence'], m_exon['sequence'])
        if homology_score > best_score:
            best_score = homology_score
            best_mouse_exon = m_exon

    if best_mouse_exon is not None:
        return {
            "human_gene_id": human_gene_id,
            "human_exon_id": h_exon['exon_id'],
            "mouse_gene_id": mouse_gene_id,
            "mouse_exon_id": best_mouse_exon['exon_id'],
            "homology_percent": round(best_score, 2),
            "mouse_homology_type": mouse_homology_type,
            "overall_homology_percent": overall_homology_percent
        }
    return None

def compare_and_pair_exons_parallel(homologous_tuple, mice_grouped, human_grouped) -> list:
    paired_exons = []
    total_genes = len(homologous_tuple)

    tasks = []
    for homologous_gene, (human_gene_id, mouse_gene_id, mouse_homology_type, overall_homology_percent) in enumerate(homologous_tuple):
        print(f"[{homologous_gene + 1}/{total_genes}] Preparing tasks for {human_gene_id} and {mouse_gene_id}...")

        if mouse_gene_id not in mice_grouped.groups or human_gene_id not in human_grouped.groups:
            continue

        mouse_exons_df = mice_grouped.get_group(mouse_gene_id)
        human_exons_df = human_grouped.get_group(human_gene_id)

        for _, h_exon in human_exons_df.iterrows():
            tasks.append((h_exon, mouse_exons_df, human_gene_id, mouse_gene_id, mouse_homology_type, overall_homology_percent))

    # 並列処理でタスクを実行
    with Pool() as pool:
        results = list(tqdm(pool.imap(calculate_homology_for_exon_pair, tasks), total=len(tasks), desc="Processing tasks"))

    # 結果を収集
    paired_exons = [result for result in results if result is not None]
    return paired_exons
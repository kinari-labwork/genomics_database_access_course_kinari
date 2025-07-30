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


def calculate_identity(seq1: str, seq2: str) -> float:
    """
    Biopythonを使って2つの配列の同一性(identity %)を計算する。
    同一性 = (一致した残基数 / アライメント長) * 100
    """
    # 配列が空、または文字列でない場合は0を返す
    if not all(isinstance(s, str) and s for s in [seq1, seq2]):
        return 0.0

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    # スコアリングは結果に影響しないが、アライメントの質のために設定
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    # alignメソッドはアライメントのリストを返す
    alignments = aligner.align(seq1, seq2)

    alignment = alignments[0]
    
    # アライメントされた配列を取得 (例: 'AT-G' と 'ATTG')
    aligned_seq1, aligned_seq2 = alignment
    
    # 一致する残基の数を数える
    matches = sum(c1 == c2 for c1, c2 in zip(aligned_seq1, aligned_seq2))
    
    # アライメントの全長（ギャップ含む）
    alignment_length = len(aligned_seq1)

    # 同一性を計算
    return (matches / alignment_length) * 100 if alignment_length > 0 else 0.0

def process_gene_pair(args):
    """
    【並列処理のワーカー関数】
    1つの遺伝子ペア（ヒトとマウス）を受け取り、その中で最も相同性の高いエクソンペアを見つける。
    """
    human_gene_id, mouse_gene_id, mouse_homology_type, overall_homology_percent, human_exons_df, mouse_exons_df = args
    
    paired_results_for_gene = []

    # このヒト遺伝子に属する各エクソンについてループ
    for _, h_exon in human_exons_df.iterrows():
        best_score = -1
        best_mouse_exon = None

        # 対応するマウス遺伝子の全エクソンと比較
        for _, m_exon in mouse_exons_df.iterrows():
            # 上で修正した identity 計算関数を呼び出す
            identity_score = calculate_identity(h_exon['sequence'], m_exon['sequence'])
            
            if identity_score > best_score:
                best_score = identity_score
                best_mouse_exon = m_exon

        if best_mouse_exon is not None:
            paired_results_for_gene.append({
                "human_gene_id": human_gene_id,
                "human_exon_id": h_exon['exon_id'],
                "mouse_gene_id": mouse_gene_id,
                "mouse_exon_id": best_mouse_exon['exon_id'],
                "homology_percent": round(best_score, 2),
                "mouse_homology_type": mouse_homology_type,
                "overall_homology_percent": overall_homology_percent
            })
            
    return paired_results_for_gene


def compare_and_pair_exons_parallel(homologous_tuple, mice_grouped, human_grouped) -> list:
    """
    【メイン関数】遺伝子ペア単位でタスクを準備し、並列処理を実行する
    """
    tasks = []
    print("Preparing tasks for parallel processing (one task per gene pair)...")
    # 遺伝子ペアごとにタスクを作成する
    for human_gene_id, mouse_gene_id, mouse_homology_type, overall_homology_percent in tqdm(homologous_tuple, desc="Preparing tasks"):
        # 必要なエクソンデータが存在するか確認
        if mouse_gene_id in mice_grouped.groups and human_gene_id in human_grouped.groups:
            mouse_exons_df = mice_grouped.get_group(mouse_gene_id)
            human_exons_df = human_grouped.get_group(human_gene_id)
            
            # 必要な情報をすべてタプルにまとめる
            tasks.append((human_gene_id, mouse_gene_id, mouse_homology_type, overall_homology_percent, human_exons_df, mouse_exons_df))

    print(f"Total {len(tasks)} gene pairs to process.")
    
    # 並列処理でタスクを実行
    with Pool() as pool:
        # pool.imap_unorderedを使うと、完了したタスクから順次結果を受け取れるため、メモリ効率が良い
        all_paired_exons = []
        for result in tqdm(pool.imap_unordered(process_gene_pair, tasks), total=len(tasks), desc="Pairing exons"):
            # 各遺伝子ペアの結果（リスト）を統合していく
            if result:
                all_paired_exons.extend(result)
    
    return all_paired_exons
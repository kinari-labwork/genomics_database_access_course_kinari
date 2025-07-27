import pandas as pd
from Bio import pairwise2

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
    homologous_tuple = tuple(zip(homologous_df["Gene stable ID"], homologous_df["Mouse gene stable ID"]))
    return homologous_tuple


def calculate_homology_percent(seq1, seq2) -> float:
    """Biopythonを使って2つの配列の相同性(%)を計算する"""
    if not isinstance(seq1, str) or not isinstance(seq2, str) or not seq1 or not seq2:
        return 0.0
    # globalxxはシンプルなスコアリング（一致:1, 不一致:0）でアライメント
    score = pairwise2.align.globalxx(seq1, seq2, score_only=True)
    # 短い方の配列長を分母にして同一性を計算
    shorter_len = min(len(seq1), len(seq2))
    return (score / shorter_len) * 100 if shorter_len > 0 else 0.0

def compare_and_pair_exons(homologous_tuple, mice_grouped, human_grouped) -> list:
    """
    Purpose: 相同遺伝子内のエキソンをペアワイズ比較し、
             各ヒトエキソンに対し最も相同性の高いマウスエキソンをペアリングする
    Input:
        homologous_tuple (list of tuple): マウスとヒトの遺伝子ペア
        mice_grouped (DataFrameGroupBy): マウスのエキソン情報（'sequence'列を含む）
        human_grouped (DataFrameGroupBy): ヒトのエキソン情報（'sequence'列を含む）
    Output:
        paired_exons (list): ペアリング結果を格納したリスト
    """
    paired_exons = []
    total_genes = len(homologous_tuple)

    # 相同遺伝子ペアでループ
    for homologous_gene, (human_gene_id, mouse_gene_id) in enumerate(homologous_tuple):
        print(f"[{homologous_gene + 1}/{total_genes}] Comparing {human_gene_id} and {mouse_gene_id}...")

        # グループが存在しない場合はスキップ
        if mouse_gene_id not in mice_grouped.groups or human_gene_id not in human_grouped.groups:
            continue

        mouse_exons_df = mice_grouped.get_group(mouse_gene_id)
        human_exons_df = human_grouped.get_group(human_gene_id)

        # 各ヒトエキソンについてループ
        for _, h_exon in human_exons_df.iterrows():
            best_score = -1
            best_mouse_exon = None

            # 全てのマウスエキソンと総当たりで比較
            for _, m_exon in mouse_exons_df.iterrows():
                # 相同性を計算
                homology_score = calculate_homology_percent(h_exon['sequence'], m_exon['sequence'])
                
                # これまでで最もスコアが高ければ更新
                if homology_score > best_score:
                    best_score = homology_score
                    best_mouse_exon = m_exon
            
            # ベストマッチが見つかった場合、そのペアの情報を保存
            if best_mouse_exon is not None:
                paired_exons.append({
                    "human_gene_id": human_gene_id,
                    "human_exon_id": h_exon['exon_id'],
                    "mouse_gene_id": mouse_gene_id,
                    "mouse_exon_id": best_mouse_exon['exon_id'],
                    "homology_percent": round(best_score, 2)
                })

    return paired_exons
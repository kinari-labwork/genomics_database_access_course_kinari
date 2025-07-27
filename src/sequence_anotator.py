import pandas as pd
from pyfaidx import Fasta

def covert_sequence_to_reverse_complement(sequence: str) -> str:
    """
    purpose: 塩基配列を逆相補鎖に変換する
    input: sequence (str): 塩基配列
    output: reverse_complement (str): 逆相補鎖の塩基配列
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}  # Nは不明な塩基
    reverse_complement = ''.join(complement[base] for base in reversed(sequence))
    return reverse_complement



def extract_exons(gtf_file: pd.DataFrame) -> pd.DataFrame:
    """
    purpose: GTFファイルからエキソン情報を抽出し、重複を除去する
    input: gtf_file (DataFrame): GTFファイルのデータ
    output: exon_data (DataFrame): エキソン情報だけを含むDataFrame
    """
    # フィーチャータイプが "exon" の行を抽出
    exon_data = gtf_file[gtf_file[2] == "exon"]
    # 重複を除去
    exon_data = exon_data.drop_duplicates(subset=[0, 3, 4, 6])
    # 必要なカラムだけを選択
    # カラム名を設定
    return exon_data

def get_exon_sequences(row: pd.Series, fasta_obj: Fasta) -> str:
    """
    purpose: エキソンの開始位置と終了位置を基に、FASTAファイルからエキソンのシーケンスを取得する
    input: row (Series): エキソンの行データ
    output: sequence (str): エキソンの塩基配列
    """
    # エキソンの開始位置と終了位置を取得
    start = row[3] - 1  # 0-indexedに変換
    end = row[4]
    chrom = row[0]
    
    # FASTAファイルからエキソンのシーケンスを取得
    sequence = str(fasta_obj[chrom][start:end])
    if row[6] == '-':
        sequence = covert_sequence_to_reverse_complement(sequence)
    return sequence

def annotate_exons_with_sequences(exon_data: pd.DataFrame, fasta_file: str) -> pd.DataFrame:
    """
    purpose: エキソンにゲノム配列をアノテーションする
    input: exon_data (DataFrame): エキソン情報を含むDataFrame
           fasta_file (str): FASTAファイルのパス
    """
    fasta = Fasta(fasta_file)
    exon_data['sequence'] = exon_data.apply(lambda row: get_exon_sequences(row, fasta), axis=1)
    return exon_data


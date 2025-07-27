import pandas as pd

def extract_genename_as_column(annotated_exons):
    """
    Pourpose: attribute列に複数の情報が格納されているので、gene nameを抽出して新しい列に格納する
    """
    # Extract gene names from the 9th column (attributes)
    annotated_exons['gene_name'] = annotated_exons[8].str.extract(r'gene_name "([^"]+)"')

    return annotated_exons
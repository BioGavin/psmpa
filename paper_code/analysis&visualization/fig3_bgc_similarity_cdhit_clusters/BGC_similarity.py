from collections import defaultdict
import re
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from itertools import combinations
import numpy as np
import logging
from pathlib import Path

def parse_cd_hit_clstr(file_path):
    cluster_dict = defaultdict(list)
    current_cluster = None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = int(line.split()[1])
            else:
                match = re.search(r'>(refseq\d+)', line)
                if match and current_cluster is not None:
                    seq_id = match.group(1)
                    cluster_dict[current_cluster].append(seq_id)
    
    return dict(cluster_dict)


def get_bgc_vector(df):
    bgc_dict = df.to_dict(orient="index")
    bgc_dict = {k: list(v.values()) for k, v in bgc_dict.items()}
    return bgc_dict



def calculate_average_bgc_count(df):
    trait_columns = df.columns[-8:]
    # 根据 16S_ID 列进行分组，并计算每组中 BGC count 列的平均值
    mean_float_bgc_count = df.groupby('16S_ID')[trait_columns].mean()
    # 对平均值进行四舍五入取整
    mean_int_bgc_count = mean_float_bgc_count.round()
    return mean_float_bgc_count, mean_int_bgc_count



def calculate_median_bgc_count(df):
    trait_columns = df.columns[-8:]
    # 根据 16S_ID 列进行分组，并计算每组中最后8列的中位数
    median_float_bgc_count = df.groupby('16S_ID')[trait_columns].median()
    # 对中位数进行四舍五入取整
    median_int_bgc_count = median_float_bgc_count.round()
    return median_float_bgc_count, median_int_bgc_count

def compute_cluster_similarities(cluster_dict, bgc_dict):
    cluster_sim_dict = {}

    for cluster_id, refseq_list in cluster_dict.items():
        if len(refseq_list) <= 1:
            continue  # 跳过孤儿簇

        vectors = [bgc_dict[refseq] for refseq in refseq_list if refseq in bgc_dict]
        if len(vectors) < 2:
            print("可能有refseq在bgc_dict中缺失")
            continue  # 可能有refseq在bgc_dict中缺失

        # 计算 pairwise cosine similarity matrix
        sim_matrix = cosine_similarity(vectors)
        pairwise_sims = [sim_matrix[i][j] for i, j in combinations(range(len(vectors)), 2)]

        sims = np.array(pairwise_sims)
        stats = {
            "mean": sims.mean(),
            "max": sims.max(),
            "min": sims.min(),
            "std": sims.std(),
            "median": np.median(sims),
            "q25": np.percentile(sims, 25),
            "q75": np.percentile(sims, 75),
        }

        cluster_sim_dict[cluster_id] = stats

    return cluster_sim_dict


def replace_zero_rows(df):
    # 需要修改的列（除了 Std）
    stat_cols = ["mean", "max", "min", "median", "q25", "q75"]

    # 识别所有这些列都为 0 的行
    zero_rows = (df[stat_cols + ["std"]] == 0).all(axis=1)

    # 将这些行的除 Std 外的列替换为 1
    df.loc[zero_rows, stat_cols] = 1.0

    return df


def main():
    output_dir = Path("cdhit_results")
    # 设置日志记录
    logging.basicConfig(filename='BGC_similarity.log', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    bgc_table = "TableS1_HQ_filtered16S.csv"
    df = pd.read_csv(bgc_table)
    mean_float_bgc_count, mean_int_bgc_count = calculate_average_bgc_count(df)
    median_float_bgc_count, median_int_bgc_count = calculate_median_bgc_count(df)
    bgc_count_dict = {
        "mean_float_bgc_count": mean_float_bgc_count,
        "mean_int_bgc_count": mean_int_bgc_count,
        "median_float_bgc_count": median_float_bgc_count,
        "median_int_bgc_count": median_int_bgc_count
    }
    
    for label, bgc_count in bgc_count_dict.items():
        bgc_vector = get_bgc_vector(bgc_count)

        threshold_ls = [80, 85] + list(range(90, 100))
        for i in threshold_ls:
            cluster_dict = parse_cd_hit_clstr(f"16SrRNA{i}.clstr")
            count = sum(1 for v in cluster_dict.values() if len(v) > 1)
            logging.info(f"16SrRNA{i} total clusters: {len(cluster_dict)}; with multiple sequences: {count}")
            cluster_sim_dict = compute_cluster_similarities(cluster_dict, bgc_vector)
            # 将结果转换为 DataFrame
            cluster_sim_df = pd.DataFrame.from_dict(cluster_sim_dict, orient="index")
            cluster_sim_df.index.name = "Cluster_ID"
            cluster_sim_df = replace_zero_rows(cluster_sim_df)
            cluster_sim_df.to_csv(output_dir / f"{label}_cdhit{i}.csv")



if __name__ == '__main__':
    main()
        




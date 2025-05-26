import os
from os import listdir
from os.path import dirname, join
from collections import namedtuple
import pandas as pd


class ReadBGCTable():
    def __init__(self, bgc_table_dir):
        self.bgc_table_dir = bgc_table_dir

    @property
    def antismash(self):
        if self.bgc_table_dir.endswith('.gz'):
            table = pd.read_csv(self.bgc_table_dir, compression='gzip', sep='\t', index_col=0)
        elif self.bgc_table_dir.endswith('.tsv'):
            table = pd.read_csv(self.bgc_table_dir, sep='\t', index_col=0)
        else:
            raise Exception("Input file format cannot be read")
        table = table.rename(columns={"Region_Num": "sum"})
        table.rename(columns=lambda c: "sum" if c.lower() == "sum" else c, inplace=True)
        return table

    @property
    def psmpa(self):
        if self.bgc_table_dir.endswith('.gz'):
            table = pd.read_csv(self.bgc_table_dir, compression='gzip', sep='\t', index_col=0)
        elif self.bgc_table_dir.endswith('.tsv'):
            table = pd.read_csv(self.bgc_table_dir, sep='\t', index_col=0)
        else:
            raise Exception("Input file format cannot be read")
        table.dropna(axis=0, how='any', inplace=True)  # 删除没有预测出结果的行
        table.rename(columns=lambda c: "sum" if c.lower() == "sum" else c,
                     inplace=True)
        return table


def drop_no_result_seq(antismash_df, psmpa_df):
    # 统一行：删除antiSMASH结果中PSMPA未预测出结果的条目
    antismash_index = antismash_df.index.to_list()
    psmpa_index = psmpa_df.index.to_list()
    diff = list(set(antismash_index).difference(set(psmpa_index)))
    antismash = antismash_df.drop(index=diff)
    return antismash


def preproduce(antismash_df, psmpa_df):
    """对antiSMASH和PSMPA的结果进行预处理。"""
    antismash_df_out = drop_no_result_seq(antismash_df, psmpa_df)
    psmpa_df_out = psmpa_df

    PreproduceDf = namedtuple('PreproduceDf', ['antismash_df', 'psmpa_df'])
    preproducedf = PreproduceDf(antismash_df_out, psmpa_df_out)
    return preproducedf


def loss(col):
    abs_col = []
    for x in col:
        x = abs(x)
        abs_col.append(x)
    y = sum(abs_col) / len(abs_col)
    return y


def sum2first(df):
    sum = df['sum']
    df.drop(columns=['sum'], inplace=True)
    df.insert(0, 'sum', sum)
    return df


def interval_statistics(df, target_col, bins, labels, counts_name):
    df[target_col] = df[target_col].abs()
    segments = pd.cut(df[target_col], bins, labels=labels, right=True)
    counts = pd.value_counts(segments, sort=False)
    counts.name = counts_name
    return counts


def calc_bias_df(psmpa_df, antismash_df, custom_order):
    # custom_order = ['sum', 'PKSI', 'PKSother', 'NRPS', 'RiPPs',
    #                 'Saccharides', 'Terpene', 'PKS-NRP_Hybrids', 'Others']  # 自定义列顺序

    bias = (psmpa_df - antismash_df).dropna(how="all", axis=0).dropna(how="all", axis=1)
    bias = bias[custom_order]
    return bias


def calc_percentage(df):
    column_sums = df.sum()
    df_percentage = df.div(column_sums, axis=1) * 100
    return df_percentage


def main(save_bias=False, bin_sum_bias=False, stat_pident=False, pident_key="V3_V4", stat_classes=False,
         classes_key="V3_V4"):
    """
    :param save_bias: bool
    :param bin_sum_bias: bool
    :param stat_pident: bool
    :param pident_key: 'full_length' or 'V3_V4'
    :return:
    """
    # 1. 基本变量
    # output_dir = join(dirname(dirname(__file__)), "plot",
    #                   os.path.basename(__file__).strip(".py") + '_output')  # 设置输出文件夹
    output_dir = join("plot", os.path.basename(__file__).strip(".py") + '_output')

    # input_dir = join(dirname(dirname(__file__)), 'results')  # 设置输入文件夹
    input_dir = 'test_results'  # 设置输入文件夹

    input_subdir = ['full_length', 'V3_V4']  # 设置输入的子文件夹
    psmpa_methods = ['emp_prob', 'mp', 'pic', 'scp', 'subtree_average', 'mean_float', 'mean_int', 'median_float',
                     'median_int']
    bgc_classes = ['sum', 'PKSI', 'PKSother', 'NRPS', 'RiPPs',
                   'Saccharides', 'Terpene', 'PKS-NRP_Hybrids', 'Others']  # 自定义列顺序

    antismash_result_file = 'test_data/5000GCA_antiSMASH_result.tsv'  # 设置antismash结果文件

    counts_dict = {}

    # 2. 创建输出文件夹
    if not os.path.exists(output_dir):
        # print(f'Warning: output folder "{output_dir}" already exists!')
        # else:
        os.mkdir(output_dir)

    # 3. 读取antiSMASH分析结果
    antismash_df = ReadBGCTable(antismash_result_file).antismash

    # 4. 读取PSMPA预测结果
    # 4.1 构造数据结构存储数据
    BiasDF = namedtuple("BiasDF", ['in_type', 'method', 'df'])
    bias_dfs = []

    # 为统计不同pident的准确性预存数据，存储不同输入序列类型两种PSMPA方法预测得到的accession，如果是PSMPA2的结果，则存储pident信息
    accession_dict = {}

    # 4.2 循环读取不同方法的PSMPA结果
    for in_type in input_subdir:
        for method in psmpa_methods:
            psmpa_result_folder = listdir(join(input_dir, in_type, method))
            psmpa_result_file_name = [f for f in psmpa_result_folder if f.startswith('psmpa')][0]
            psmpa_result_file = join(join(input_dir, in_type, method), psmpa_result_file_name)
            psmpa_df = ReadBGCTable(psmpa_result_file).psmpa  # 读取PSMPA预测结果
            if 'pident' in psmpa_df.columns:
                accession_dict.setdefault(f'{in_type}.psmpa2', psmpa_df['pident'])
            else:
                accession_dict.setdefault(f'{in_type}.psmpa1', psmpa_df.index)

            bias = calc_bias_df(psmpa_df, antismash_df, bgc_classes)  # 计算bias

            # counts_dict.setdefault(in_type, len(bias))
            counts_dict[f"{in_type}_{method}"] = len(bias)

            if save_bias:
                if not os.path.exists(join(output_dir, in_type)):
                    os.mkdir(os.path.join(output_dir, in_type))
                    # pass
                    # print(f'Warning: output subfolder "{os.path.join(output_dir, in_type)}" already exists!')
                # else:
                #     os.mkdir(os.path.join(output_dir, in_type))
                bias_name = f'{method}_bias.tsv.gz'
                bias.to_csv(os.path.join(output_dir, in_type, bias_name), sep='\t', index=True, compression='gzip')
            bias_df = BiasDF(in_type=in_type, method=method, df=bias)
            bias_dfs.append(bias_df)

    # 5. 对sum偏差进行区间统计
    if bin_sum_bias:
        counts_series_dict = {"full_length": [],
                              "V3_V4": []}
        bins = [-float('inf'), 1, 5, float('inf')]
        labels = ['[0,1]', '(1,5]', '(5,inf]']

        counts_name_dict = {}
        for method in psmpa_methods:
            if method in ['mean_float', 'mean_int', 'median_float', 'median_int']:
                counts_name_dict[method] = f'psmpa2_{method}'
            else:
                counts_name_dict[method] = f'psmpa1_{method}'

        for df in bias_dfs:
            counts_series = interval_statistics(df.df, target_col='sum', bins=bins, labels=labels,
                                                counts_name=counts_name_dict[df.method])
            counts_series_dict[df.in_type].append(counts_series)

        for in_type, v in counts_series_dict.items():
            sum_counts_df = pd.concat(v, axis=1)
            sum_counts_df.index.name = 'sum'
            sum_counts_df_percentage = calc_percentage(sum_counts_df)
            sum_counts_df_percentage.to_csv(join(output_dir, f'BGC_{in_type}.tsv'), sep='\t')

    # 6. 对不同pident的准确性进行统计
    if stat_pident:
        if pident_key:
            psmpa1_accession = accession_dict[f'{pident_key}.psmpa1']
            psmpa2_accession_pident = accession_dict[f'{pident_key}.psmpa2']
            target_accession = psmpa2_accession_pident.index.intersection(psmpa1_accession)

            pident_bins = [0, 80, 90, 95, float('inf')]
            pident_bins_labels = ['0-80', '80-90', '90-95', '95-100']

            # 使用cut函数进行分组
            target_accession_pident = pd.cut(psmpa2_accession_pident[target_accession], bins=pident_bins,
                                             labels=pident_bins_labels, right=False)
            for pident_bin in pident_bins_labels:
                pident_bin_accession = target_accession_pident[target_accession_pident == pident_bin].index

                counts_dict.setdefault(pident_bin, len(pident_bin_accession))
                bins = [-float('inf'), 1, 5, float('inf')]
                labels = ['[0,1]', '(1,5]', '(5,inf]']

                counts_series_ls = []
                counts_name_dict = {}
                for method in psmpa_methods:
                    if method in ['mean_float', 'mean_int', 'median_float', 'median_int']:
                        counts_name_dict[method] = f'psmpa2_{method}'
                    else:
                        counts_name_dict[method] = f'psmpa1_{method}'
                # print(pident_bin_accession)
                for df in bias_dfs:  # 不能对所有bias进行操作，要选择
                    if df.in_type == pident_key:
                        counts_series = interval_statistics(df.df.loc[pident_bin_accession], target_col='sum',
                                                            bins=bins,
                                                            labels=labels, counts_name=counts_name_dict[df.method])
                        counts_series_ls.append(counts_series)
                #
                sum_counts_df = pd.concat(counts_series_ls, axis=1)
                sum_counts_df.index.name = 'sum'
                sum_counts_df_percentage = calc_percentage(sum_counts_df)
                sum_counts_df_percentage.to_csv(join(output_dir, f'pident_{pident_bin}.tsv'), sep='\t')

        else:
            raise Exception(
                "Please choose a value from either 'full_length' or 'V3_V4' to set as the value of pident_key. ")

    # 7. 对每个类别的BGC进行统计
    if stat_classes:
        if classes_key:

            bins = [-float('inf'), 1, float('inf')]
            labels = ['[0,1]', '(1,inf]']

            for bgc_class in bgc_classes[1:]:
                counts_series_ls = []
                counts_name_dict = {}
                for method in psmpa_methods:
                    if method in ['mean_float', 'mean_int', 'median_float', 'median_int']:
                        counts_name_dict[method] = f'psmpa2_{method}'
                    else:
                        counts_name_dict[method] = f'psmpa1_{method}'
                for df in bias_dfs:  # 不能对所有bias进行操作，要选择
                    if df.in_type == classes_key:
                        counts_dict.setdefault(bgc_class, len(df.df))
                        # print(bgc_class, len(df.df))
                        counts_series = interval_statistics(df.df, target_col=bgc_class, bins=bins,
                                                            labels=labels, counts_name=counts_name_dict[df.method])
                        counts_series_ls.append(counts_series)
                #
                sum_counts_df = pd.concat(counts_series_ls, axis=1)
                sum_counts_df.index.name = 'sum'
                sum_counts_df_percentage = calc_percentage(sum_counts_df)
                sum_counts_df_percentage.to_csv(join(output_dir, f'{bgc_class}.tsv'), sep='\t')

                # print(bgc_class)
        else:
            raise Exception(
                "Please choose a value from either 'full_length' or 'V3_V4' to set as the value of classes_key. ")

    # 统计每个分析中预测的结果数量
    counts_res = pd.Series(counts_dict)
    counts_res.name = "counts"
    counts_res.to_csv(os.path.join(output_dir, f'counts.tsv'), sep='\t')


if __name__ == '__main__':
    main(save_bias=False, bin_sum_bias=True,
         stat_pident=True, pident_key="V3_V4",
         stat_classes=True, classes_key="V3_V4")

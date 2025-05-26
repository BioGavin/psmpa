# BGC similarity after cdhit

- 对16S序列按不同阈值聚类

```bash
# 阈值90以上的聚类
for c in $(seq 0.90 0.01 0.99); do cd-hit-est -i filtered_intergrated_16SrRNA.fasta -o 16SrRNA${c#0.} -c $c -n 9; done
# 阈值90以上的聚类
cd-hit-est -i filtered_intergrated_16SrRNA.fasta -o 16SrRNA85 -c 0.85 -n 6
cd-hit-est -i filtered_intergrated_16SrRNA.fasta -o 16SrRNA80 -c 0.80 -n 5
# invalid clstr threshold, should >=0.8
```

- 编写脚本根据cdhit结果计算簇内BGC相似度

```bash
python3 BGC_similarity.py
```

- 绘图详见 `plot.ipynb`

# *psmpa* 基础数据库的构建流程

## pump 数据库内容确定

### **Metadata**

以 2021-07-11 NCBI-Assembly 数据库中 RefSeq 的 Metadata 为依据建立相关数据库，Metadata中共有220524个基因组信息。

### **Taxon**

根据 Metadata，通过 taxonkit (v0.8.0) 获取所有菌株的 Lineage 信息。

### **16S rRNA database**

根据 Metadata，通过 NCBI 直接下载或 prokka (v1.12) 分析获取 16S rRNA。

### **antiSMASH analysis**

根据 Metadata，下载对应的组装基因组数据，通过 antiSMASH (v5.2.0) 分析获取 BGC 分布情况。

### **psmpa basedata**

根据以上信息的完整性进行筛选，以同时具有 Lineage 信息、16S rRNA、antiSMASH 分析结果的菌株构建 psmpa1 数据库。

**最终，共有216408个菌株形成 psmpa 数据库。**



## File Catalog

- metadata_20210711.tsv：2021-07-11 NCBI-Assembly 数据库中 RefSeq 的 metadata
- psmpa_metadata.tsv：psmpa 数据库中216408个菌株的 metadata
- psmpa_antismash_result.tsv：psmpa 数据库中216408个菌株基因组的 antiSMASH 分析结果

- antismash_5.2.0_result.tsv.gz: save all antiSMASH analysis results.

```shell
# antiSMASH command
antismash --genefinding-tool prodigal genome.fna
```


# PSMPA

[![PyPI version](https://img.shields.io/badge/pypi%20package-1.1.1-brightgreen)](https://pypi.org/project/psmpa/) [![Licence](https://img.shields.io/badge/licence-GPLv3-blue)](https://opensource.org/licenses/GPL-3.0/) [![Web](https://img.shields.io/badge/version-web-red)](https://www.psmpa.net/analysis)

PSMPA is a Python pipeline to predict secondary metabolism potential using amplicans  for a single strain or microbial communities. We provide a web version for the convenience of users, you can [click here](https://www.psmpa.net/analysis).

![]( https://cdn.jsdelivr.net/gh/BioGavin/Pic/imgpsmpa_logo2.png)

# Requirements

Specific libraries are required for PSMPA. We provide a requirements file to install everything at once.
Here, We recommende Conda for environment deployment.

```shell
conda env create -f environment.yml
```
If you have successfully created this environment, don't forget to activate it.
```shell
conda activate psmpa
```

# Installation & Help
Install the PSMPA package to the environment.
```shell
pip install psmpa
```
So far, if you installed successfully, you can run this command for more help information.
```shell
psmpa1 -h
```
or
```shell
psmpa2 -h
```

# Sample
## *psmpa1*

Usage:
- for single or multiple bacteria analysis
```shell
psmpa1 -s 16S.fna -o psmpa1_test_out
```

- for microbiological samples analysis
```shell
psmpa1 -s sequences.fasta -i feature-table.biom -o psmpa1_sample_test_out
```


## *psmpa2*
Usage:
- for single or multiple bacteria analysis
```shell
psmpa2 -s 16S.fna -o psmpa2_test_out
```

- for microbiological samples analysis
```shell
psmpa2 -s sequences.fasta -i feature-table.biom -o psmpa2_sample_test_out
```



# Explanation

If empty rows appear in the BGCs predicted result, the likely reason is that the BLAST analysis did not match any sequences. So, if sample analysis is performed, sequences with no results are ignored.


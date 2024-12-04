# PSMPA-Atlas-Exploration

The PSMPA-Atlas data was extracted by the script [`ex_antismash_bgc.py`](scripts/ex_antismash_bgc.py) and stored in [TableS1.tsv](data/TableS1.tsv).

The exploratory analysis of PSMPA-Atlas mentioned in the paper can be found in the [PSMPA-Atlas-Investigation.ipynb](notebooks/PSMPA-Atlas-Investigation.ipynb).

We used [`as2bs.py`](scripts/as2bs.py) to convert antiSMASH BGC classes into BiG-SCAPE BGC classes, for example:

```bash
python3 scripts/as2bs.py -i data/TableS1.tsv -o data/TableS1.8BGC.tsv
```




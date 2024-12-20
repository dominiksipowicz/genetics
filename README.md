# CSV headers

CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,NG173LPFBH

CHROM - chromosome
POS - position
ID - variant ID
REF - reference allele
ALT - alternative allele
QUAL - quality score
FILTER - filter status
INFO - additional information
FORMAT - format

> I use Nebula to sequence my DNA (https://portal.nebula.org/)

# Steps

1. Download the VCF file from Nebula (compressed)
2. Unzip the file VCF.gz

```bash
gunzip VCF.gz
```

3. Run stats for one VCF file and save it to a CSV file

```bash
python gene_stats.py

```

## example table CSV

| CHROM | POS   | ID           | REF                    | ALT | QUAL               | FILTER   | INFO                                       | FORMAT                                                                        | SAMPLE        |
| ----- | ----- | ------------ | ---------------------- | --- | ------------------ | -------- | ------------------------------------------ | ----------------------------------------------------------------------------- | ------------- |
| chr1  | 10616 | rs376342519  | CCGCCGTTGCAAAGGCGCGCCG | C   | 8.40999984741211   |          | <cyvcf2.cyvcf2.INFO object at 0x11f856160> | ['GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'GQ', 'PL', 'GP', 'PRI', 'SB', 'MB'] | [1, 1, False] |
| chr1  | 10816 | rs1266288166 | C                      | CCA | 10.470000267028809 |          | <cyvcf2.cyvcf2.INFO object at 0x11f856190> | ['GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'GQ', 'PL', 'GP', 'PRI', 'SB', 'MB'] | [1, 1, False] |
| chr1  | 11409 | rs9803797    | A                      | G   | 13.770000457763672 |          | <cyvcf2.cyvcf2.INFO object at 0x11f8560a0> | ['GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'GQ', 'PL', 'GP', 'PRI', 'SB', 'MB'] | [1, 1, False] |
| chr1  | 13302 | rs75241669   | C                      | T   | 9.34000015258789   | LowDepth | <cyvcf2.cyvcf2.INFO object at 0x11f856070> | ['GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'GQ', 'PL', 'GP', 'PRI', 'SB', 'MB'] | [1, 1, False] |
| chr1  | 13327 | rs2691329    | G                      | C   | 8.069999694824219  | LowDepth | <cyvcf2.cyvcf2.INFO object at 0x11f856040> | ['GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'GQ', 'PL', 'GP', 'PRI', 'SB', 'MB'] | [0, 1, False] |
| chr1  | 14673 | rs4690       | G                      | C   | 7.690000057220459  |          | <cyvcf2.cyvcf2.INFO object at 0x11f855d10> | ['GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'GQ', 'PL', 'GP', 'PRI', 'SB', 'MB'] | [0, 1, False] |

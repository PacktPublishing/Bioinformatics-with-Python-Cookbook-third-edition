# Chapter 1 – rpy2 and the 1000 Genomes Project

This chapter contains my own reimplementation and notes based on the original *Bioinformatics-with-Python-Cookbook-third-edition*, adapted for current tools and data sources.

---

## Enhancements and additions

- Added an explicit download and verification step in Python for the 1000 Genomes sequence index:

```python
url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.sequence.index"
```
- Demonstrated integration with R via the rpy2 interface.

- Converted the R data frame to pandas format for downstream analysis.

- Applied explicit numeric type coercion and axis scaling to improve plot clarity.

- Ensured ggplot objects are explicitly rendered to an R graphics device for reproducible output.

## Observations

Two main points were addressed while reproducing the original workflow:

1. Column type handling

BASE_COUNT and READ_COUNT represent numeric quantities but may be imported as character vectors, which can affect ggplot2 plotting behaviour if not coerced.

2. Explicit plot rendering

When using ggplot2 via rpy2, plots are not automatically rendered to a file. Explicitly printing the ggplot object to an R graphics device ensures the expected output is produced.

## Key takeaway

This chapter illustrates common interoperability considerations when bridging Python and R workflows and highlights the importance of explicit data validation and visualization control for reproducible bioinformatics analyses.
## Chapter 1 – rpy2 and the 1000 Genomes Project

This chapter contains my own reimplementation and notes based on the
original *Bioinformatics-with-Python-Cookbook-third-edition*, adapted for current tools and data sources.

---

## What was changed from the original cookbook

- Fixed a data loading failure caused by outdated FTP links.
- Added an explicit download and verification step in Python.
### Data download

The 1000 Genomes sequence index is downloaded at runtime:

```python
url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.sequence.index"
```
- Demonstrated safe integration with R via the `rpy2` interface.
- Converted the R data frame to pandas format for downstream analysis.

---

## Issues identified

Two main issues were encountered when reproducing the original figure:

### 1. Column type coercion

Although `BASE_COUNT` and `READ_COUNT` represent numeric quantities, they were
imported as character vectors. This caused ggplot2 to treat the axes as
discrete rather than continuous, leading to incorrect scaling and axis
behaviour.

### 2. Implicit plotting behaviour

When using ggplot2 via `rpy2`, plots are not automatically rendered to file.
The ggplot object must be explicitly printed to an active R graphics device
to produce the expected output.

---

## Resolution

To ensure correct plotting behaviour:

- `BASE_COUNT` and `READ_COUNT` were explicitly coerced to numeric types in R.
- Axis handling was corrected to enforce continuous numeric scales.
- The ggplot object was explicitly printed to a PNG graphics device.

These changes resulted in a plot that correctly reproduces the expected
separation between YRI and CEU samples and the distinct sequencing coverage
groups.

---

## Key takeaway

This chapter illustrates common interoperability pitfalls when bridging
Python and R workflows and highlights the importance of explicit data
validation and visualization control.

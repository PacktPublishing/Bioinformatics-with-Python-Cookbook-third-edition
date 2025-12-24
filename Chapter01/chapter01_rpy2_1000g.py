"""
Chapter 1 – Using rpy2 with the 1000 Genomes Project

This script demonstrates:
- Downloading the 1000 Genomes sequence index
- Interfacing Python with R via rpy2
- Subsetting data in R
- Handling type coercion issues
- Reproducing a ggplot2 scatter plot via rpy2
"""

import os
import requests
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects import r

# -------------------------------------------------------------------
# R environment setup (Windows example)
# -------------------------------------------------------------------

os.environ["R_HOME"] = r"C:\Program Files\R\R-4.4.2"
os.environ["PATH"] += os.pathsep + r"C:\Program Files\R\R-4.4.2\bin\x64"

# -------------------------------------------------------------------
# Download 1000 Genomes sequence index (explicit + reproducible)
# -------------------------------------------------------------------

URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
    "phase3/20130502.phase3.sequence.index"
)
OUTFILE = "20130502.phase3.sequence.index"

if not os.path.exists(OUTFILE):
    response = requests.get(URL)
    response.raise_for_status()
    with open(OUTFILE, "wb") as f:
        f.write(response.content)

# -------------------------------------------------------------------
# Load data into R
# -------------------------------------------------------------------

read_delim = r("read.delim")

r.assign(
    "seq.data",
    read_delim(
        OUTFILE,
        header=True,
        stringsAsFactors=False
    )
)

# -------------------------------------------------------------------
# Explicit type coercion (critical fix)
# -------------------------------------------------------------------

r("""
seq.data$BASE_COUNT <- as.numeric(seq.data$BASE_COUNT)
seq.data$READ_COUNT <- as.numeric(seq.data$READ_COUNT)
""")

# -------------------------------------------------------------------
# Subset YRI and CEU populations
# -------------------------------------------------------------------

r("""
yri_ceu <- seq.data[
    seq.data$POPULATION %in% c("YRI", "CEU") &
    seq.data$BASE_COUNT < 2e9 &
    seq.data$READ_COUNT < 3e7,
]
""")

yri_ceu = r("yri_ceu")

# -------------------------------------------------------------------
# Build ggplot2 scatter plot
# -------------------------------------------------------------------

scatter = (
    ggplot2.ggplot(yri_ceu) +
    ggplot2.aes_string(
        x="BASE_COUNT",
        y="READ_COUNT",
        shape="factor(POPULATION)",
        col="factor(ANALYSIS_GROUP)"
    ) +
    ggplot2.geom_point()
)

# -------------------------------------------------------------------
# Render plot to file (explicit print required in rpy2)
# -------------------------------------------------------------------

r.png(
    "out.png",
    width=1200,
    height=900,
    res=150,
    type="cairo-png"
)

r("print")(scatter)
r("dev.off")()

print("Plot saved to out.png")

#!/usr/bin/env python
# coding: utf-8

# # Chapter 1 – Using rpy2 with the 1000 Genomes Project
# 
# This chapter demonstrates how to interface Python with R using `rpy2`, and highlights common pitfalls when working with external bioinformatics
# resources such as the 1000 Genomes Project.
# 




import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.4.2"
os.environ["PATH"] += os.pathsep + r"C:\Program Files\R\R-4.4.2\bin\x64"
import os
from IPython.display import Image
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.functions import SignatureTranslatedFunction
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects import conversion




read_delim = robjects.r('read.delim')
seq_data = read_delim('sequence.index', header=True,stringsAsFactors=False)




print('This dataframe has %d columns and %d rows' 
% (seq_data.ncol, seq_data.nrow))
print(seq_data.colnames)


# ## Issue: read.delim returns empty data
# 
# When running the original cookbook code, `read.delim()` returned a dataframe with 0 rows and 1 column:
# 
# X404..Not.FoundECHO.is.on.
# 
# This indicates that the remote resource was no longer accessible, and R silently parsed an HTTP 404 response instead of tabular data.
# 

# ## Robust solution: download data explicitly in Python
# 
# To avoid silent failures, the sequence index file was downloaded explicitly using Python before being read into R via `rpy2`.
# 




import requests

url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.sequence.index"
response = requests.get(url)

# Save to file
with open("20130502.phase3.sequence.index", "wb") as f:
    f.write(response.content)

print("File downloaded successfully.")




with open("20130502.phase3.sequence.index", "r") as f:
    for i in range(10):  # Show the first 10 lines
        print(f.readline())





read_delim = robjects.r('read.delim')
seq_data = read_delim('20130502.phase3.sequence.index', header=True,stringsAsFactors=False)





print('This dataframe has %d columns and %d rows' 
% (seq_data.ncol, seq_data.nrow))
print(seq_data.colnames)





as_integer = robjects.r('as.integer')
match = robjects.r.match




my_col = match('READ_COUNT', seq_data.colnames)[0] # vector returned
print('Type of read count before as.integer: %s' % seq_data[my_col - 1].rclass[0])
seq_data[my_col - 1] = as_integer(seq_data[my_col - 1])
print('Type of read count after as.integer: %s' % seq_data[my_col - 1].rclass[0])




robjects.r.assign('seq.data', seq_data)




from rpy2.robjects.functions import SignatureTranslatedFunction




ggplot2.theme = SignatureTranslatedFunction(ggplot2. theme, init_prm_translate = {'axis_text_x': 'axis. text.x'})




bar = ggplot2.ggplot(seq_data) + \
      ggplot2.aes_string(x='CENTER_NAME') + \
      ggplot2.geom_bar() + \
      ggplot2.theme(axis_text_x=ggplot2.element_text(angle=90, hjust=1))

robjects.r.png('out.png', type='cairo-png')
bar.plot()
dev_off = robjects.r('dev.off')
dev_off()




Image(filename='out.png')



robjects.r('yri_ceu <- seq.data[seq.data$POPULATION %in% c("YRI", "CEU") & seq.data$BASE_COUNT < 2E09 & seq.data$READ_COUNT < 3E07, ]')
yri_ceu = robjects.r('yri_ceu')




scatter = (ggplot2.ggplot(yri_ceu) +
           ggplot2.aes_string(x='BASE_COUNT', y='READ_COUNT',
                              shape='factor(POPULATION)', col='factor(ANALYSIS_GROUP)') +
           ggplot2.geom_point())

robjects.r.png('out1.png', type='cairo-png')
scatter.plot()
robjects.r('dev.off')()




Image(filename='out1.png')




robjects.r("""
seq.data$BASE_COUNT <- as.numeric(seq.data$BASE_COUNT)
seq.data$READ_COUNT <- as.numeric(seq.data$READ_COUNT)
""")




robjects.r("""
yri_ceu <- seq.data[
    seq.data$POPULATION %in% c("YRI", "CEU") &
    seq.data$BASE_COUNT < 2e9 &
    seq.data$READ_COUNT < 3e7,
]
""")

yri_ceu = robjects.r("yri_ceu")




scatter = (
    ggplot2.ggplot(yri_ceu) +
    ggplot2.aes_string(
        x="BASE_COUNT",
        y="READ_COUNT",
        shape="factor(POPULATION)",
        col="factor(ANALYSIS_GROUP)"
    ) +
    ggplot2.geom_point() +
    ggplot2.scale_x_continuous(
        breaks=robjects.FloatVector([0, 5e8, 1e9, 1.5e9, 2e9])
    ) +
    ggplot2.scale_y_continuous(
        breaks=robjects.FloatVector([0, 1e7, 2e7, 3e7])
    )
)




robjects.r("str(yri_ceu$BASE_COUNT)")
robjects.r("str(yri_ceu$READ_COUNT)")




from rpy2.robjects import r

r.png("out2.png", width=1200, height=900, res=150, type="cairo-png")
r("print")(scatter)
r("dev.off")()




Image(filename='out2.png')







# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## The cell below will get the data file, you only need to run it once 

# %% [markdown]
# (you do not need to do this if you have done it in the Interfacing_R notebook)

# %%
# !rm sequence.index 2>/dev/null
# !wget -nd http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.sequence.index -O sequence.index

# %%
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2

# %load_ext rpy2.ipython

# %% language="R"
# seq.data <- read.delim('sequence.index', header=TRUE, stringsAsFactors=FALSE)
# seq.data$READ_COUNT <- as.integer(seq.data$READ_COUNT)
# seq.data$BASE_COUNT <- as.integer(seq.data$BASE_COUNT)

# %%
# seq_data = %R seq.data
print(type(seq_data))  #pandas dataframe???

# %%
my_col = list(seq_data.columns).index("CENTER_NAME")
seq_data['CENTER_NAME'] = seq_data['CENTER_NAME'].apply(lambda x: x.upper())

# %%
# %R -i seq_data
# %R print(colnames(seq_data))

# %% language="R"
# seq_data <- seq_data[seq_data$WITHDRAWN==0, ]
# seq_data$POPULATION <- as.factor(seq_data$POPULATION)

# %% language="R"
# bar <- ggplot(seq_data) +  aes(factor(CENTER_NAME)) + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# print(bar)

# %% language="R"
# seq_data$POPULATION <- as.factor(seq_data$POPULATION)
# yri_ceu <- seq_data[seq_data$POPULATION %in% c("YRI", "CEU") & seq_data$BASE_COUNT < 2E9 & seq_data$READ_COUNT < 3E7, ]

# %% language="R"
# scatter <- ggplot(yri_ceu, aes(x=BASE_COUNT, y=READ_COUNT, col=factor(ANALYSIS_GROUP), shape=POPULATION)) + geom_point()
# print(scatter)

# %% language="R"
# library(gridExtra)
# library(grid)
# g <- grid.arrange(bar, scatter, ncol=1)
# g

# %% language="R"
# png('fig.png')
# g
# dev.off()

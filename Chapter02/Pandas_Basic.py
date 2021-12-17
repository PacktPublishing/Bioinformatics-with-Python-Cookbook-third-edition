# # Using Pandas to process vaccine adverse events
#
# ## Data Access
#
# Go to https://vaers.hhs.gov/data/datasets.html and Download 2021 **zip** Data. Please do not download only the CSV File.
#
# Drop it on the directory where this notebook is.


# !unzip 2021VAERSData.zip
# !gzip -9 *csv


import pandas as pd
import matplotlib.pyplot as plt

vdata = pd.read_csv(
    "2021VAERSDATA.csv.gz", encoding="iso-8859-1")

vdata.columns

vdata.dtypes

vdata.shape

vdata.iloc[0]

vdata = vdata.set_index("VAERS_ID")

vdata.loc[916600]

vdata.head(3)

vdata.iloc[:3]

vdata.iloc[:5, 2:4]

vdata["AGE_YRS"].max()

vdata.AGE_YRS.max()

vdata["AGE_YRS"].sort_values().plot(use_index=False)

vdata["AGE_YRS"].sort_values().plot(use_index=False)

fig, ax = plt.subplots(1, 2, sharey=True, dpi=300)
fig.suptitle("Age of adverse events")
vdata["AGE_YRS"].sort_values().plot(
    use_index=False, ax=ax[0],
    xlabel="Obervation", ylabel="Age")
vdata["AGE_YRS"].plot.hist(bins=20, orientation="horizontal")
fig.savefig("adverse.png")

vdata["AGE_YRS"].dropna().apply(lambda x: int(x)).value_counts()
# not documented

vdata.DIED.value_counts(dropna=False)
# NA is a problem, how to be implemented


vdata["is_dead"] = (vdata.DIED == "Y")


dead = vdata[vdata.is_dead]
vax = pd.read_csv("2021VAERSVAX.csv.gz", encoding="iso-8859-1").set_index("VAERS_ID")
print(vax.columns)
print(vax.shape)
print(vax.VAX_TYPE.unique())

vax.groupby("VAX_TYPE").size().sort_values()

vax19 = vax[vax.VAX_TYPE == "COVID19"]
vax19_dead = dead.join(vax19)
# join on id, discuss
vax19_dead.index.value_counts()

baddies = vax19_dead.groupby("VAX_LOT").size().sort_values(ascending=False)
for i, (lot, cnt) in enumerate(baddies.items()):
    print(lot, cnt, len(vax19_dead[vax19_dead.VAX_LOT == lot].groupby("STATE")))
    if i == 10:
        break


# The data above is not totally correct - at least in terms of interpretation, but for that we need to check the next recipe

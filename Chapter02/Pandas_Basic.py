# # Using Pandas to process vaccine adverse events
#
# ## Data Access
#
# Go to https://vaers.hhs.gov/data/datasets.html and Download 2021 **zip** Data. Please do not download only the CSV File.
#
# Drop it on the directory where this notebook is.


# !unzip 2021VAERSData.zip
# !gzip -9 *csv


import gzip
import pandas as pd
import matplotlib.pyplot as plt


vdata = pd.read_csv("2021VAERSDATA.csv.gz", encoding="iso-8859-1")

vdata.columns

vdata.dtypes

vdata.iloc[0]  # vdata.loc[0]

vdata = vdata.set_index("VAERS_ID")
vdata.loc[916600]

vdata.head(5)

vdata.iloc[:5]

vdata.iloc[:5,:3]

vdata["AGE_YRS"].max()

vdata.AGE_YRS.max()

vdata["AGE_YRS"].unique()

vdata["AGE_YRS"].sort_values().plot(use_index=False)

vdata["AGE_YRS"].plot.hist(bins=20)

fig, ax = plt.subplots(1, 2, sharey=True)
fig.suptitle("Age of adverse events")
vdata["AGE_YRS"].sort_values().plot(
    use_index=False, ax=ax[0],
    xlabel="Obervation", ylabel="Age")
vdata["AGE_YRS"].plot.hist(bins=20, orientation="horizontal")

vdata["AGE_YRS"].value_counts()

vdata["AGE_YRS"].value_counts(dropna=False)

vdata["AGE_YRS"].dropna().apply(lambda x: int(x)).value_counts()

vdata.DIED.value_counts(dropna=False)


vdata["is_dead"] = (vdata.DIED == "Y")


dead = vdata[vdata.is_dead]
dead


vax = pd.read_csv("2021VAERSVAX.csv.gz", encoding="iso-8859-1").set_index("VAERS_ID")
print(vax.columns)
print(vax.shape)
print(vax.VAX_TYPE.unique())

vax.iloc[0]

vax.groupby("VAX_TYPE").size().sort_values()

vax19 = vax[vax.VAX_TYPE == "COVID19"]

vax_dead = vax.join(dead)
# join on id, discuss

baddies = vax_dead.groupby("VAX_LOT").size().sort_values(ascending=False)
for i, (lot, cnt) in enumerate(baddies.items()):
    print(lot, cnt, len(vax_dead[vax_dead.VAX_LOT == lot].groupby("STATE")))
    if i == 10:
        break

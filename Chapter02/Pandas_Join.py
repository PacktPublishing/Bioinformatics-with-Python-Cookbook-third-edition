# # Pandas advanced

import numpy as np
import pandas as pd

# # Code to sample original data
#
# ```
# vdata = pd.read_csv("2021VAERSDATA.csv.gz", encoding="iso-8859-1")
# vdata.sample(frac=0.9).to_csv("vdata_sample.csv.gz", index=False)
# vax = pd.read_csv("2021VAERSVAX.csv.gz", encoding="iso-8859-1")
# vax.sample(frac=0.9).to_csv("vax_sample.csv.gz", index=False)
# ```

vdata = pd.read_csv("vdata_sample.csv.gz") # No encoding
vax = pd.read_csv("vax_sample.csv.gz")

vdata_with_vax = vdata.join(
    vax.set_index("VAERS_ID"),
    on="VAERS_ID",
    how="inner")

len(vdata), len(vax), len(vdata_with_vax)

lost_vdata = vdata.loc[~vdata.index.isin(vdata_with_vax.index)]
lost_vdata

lost_vax = vax[~vax["VAERS_ID"].isin(vdata_with_vax["VAERS_ID"])]
lost_vax


# Left, Right and outer caveats


vdata_with_vax_left = vdata.join(
    vax.set_index("VAERS_ID"),
    on="VAERS_ID")

vdata_with_vax_left.groupby("VAERS_ID").size().sort_values()

len(vdata_with_vax_left), len(vdata_with_vax_left.VAERS_ID.unique())

# +
#vdata_all = pd.read_csv("2021VAERSDATA.csv.gz", encoding="iso-8859-1")
#vax_all = pd.read_csv("2021VAERSVAX.csv.gz", encoding="iso-8859-1")
# -

dead = vdata[vdata.DIED == "Y"]
vax19 = vax[vax.VAX_TYPE == "COVID19"]
vax19_dead = vax19.join(dead.set_index("VAERS_ID"), on="VAERS_ID", how="right")
# join on id, discuss

len(vax19), len(dead), len(vax19_dead)

len(vax19_dead[vax19_dead.VAERS_ID.duplicated()])

len(vax19_dead) - len(dead)

vax19_dead["STATE"] = vax19_dead["STATE"].str.upper()
dead_lot = vax19_dead[["VAERS_ID", "VAX_LOT", "STATE"]].set_index(["VAERS_ID", "VAX_LOT"])
dead_lot_clean = dead_lot[~dead_lot.index.duplicated()]
dead_lot_clean = dead_lot_clean.reset_index()
dead_lot_clean[dead_lot_clean.VAERS_ID.isna()]

baddies = dead_lot_clean.groupby("VAX_LOT").size().sort_values(ascending=False)
for i, (lot, cnt) in enumerate(baddies.items()):
    print(lot, cnt, len(dead_lot_clean[dead_lot_clean.VAX_LOT == lot].groupby("STATE")))
    if i == 10:
        break

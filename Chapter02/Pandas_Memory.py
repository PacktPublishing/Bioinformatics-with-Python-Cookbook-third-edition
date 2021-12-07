# # Pandas advanced

import numpy as np
import pandas as pd

vdata = pd.read_csv("2021VAERSDATA.csv.gz", encoding="iso-8859-1")

vdata.info(memory_usage="deep")

for name in vdata.columns:
    col_bytes = vdata[name].memory_usage(index=False, deep=True)
    col_type = vdata[name].dtype
    print(
        name,
        col_type, col_bytes // (1024 ** 2))

vdata.DIED.memory_usage(index=False, deep=True)

vdata.DIED.fillna(False).astype(bool).memory_usage(index=False, deep=True)

vdata.STATE.unique()

vdata["STATE"] = vdata.STATE.str.upper()

states = list(vdata["STATE"].unique())
states

vdata["encoded_state"] = vdata.STATE.apply(lambda state: states.index(state))
vdata["encoded_state"] = vdata["encoded_state"].astype(np.uint8)

vdata[["encoded_state", "STATE"]].head(10)

vdata["STATE"].memory_usage(index=False, deep=True)

vdata["encoded_state"].memory_usage(index=False, deep=True)

vdata.index

states = list(pd.read_csv(
    "vdata_sample.csv.gz",
    converters={
       "STATE": lambda state: state.upper()  # You need to know the states in advance
    },
    usecols=["STATE"]
)["STATE"].unique())

vdata = pd.read_csv(
    "vdata_sample.csv.gz",
    index_col="VAERS_ID",
    converters={
       "DIED": lambda died: died == "Y",
       "STATE": lambda state: states.index(state.upper())
    },
    usecols=lambda name: name != "SYMPTOM_TEXT"
)
vdata["STATE"] = vdata["STATE"].astype(np.uint8)
vdata.info(memory_usage="deep")



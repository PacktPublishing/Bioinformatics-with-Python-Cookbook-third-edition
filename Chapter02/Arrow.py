import gzip
import pandas as pd
from pyarrow import csv
import pyarrow.compute as pc

vdata_pd = pd.read_csv("2021VAERSDATA.csv.gz", encoding="iso-8859-1")
columns = list(vdata_pd.columns)
vdata_pd.info(memory_usage="deep")

vdata_arrow = csv.read_csv("2021VAERSDATA.csv.gz")
tot_bytes = sum([
    vdata_arrow[name].nbytes
    for name in vdata_arrow.column_names])
print(f"Total {tot_bytes // (1024 ** 2)} MB")

for name in vdata_arrow.column_names:
    arr_bytes = vdata_arrow[name].nbytes
    arr_type = vdata_arrow[name].type
    pd_bytes = vdata_pd[name].memory_usage(index=False, deep=True)
    pd_type = vdata_pd[name].dtype
    print(
        name,
        arr_type, arr_bytes // (1024 ** 2),
        pd_type, pd_bytes // (1024 ** 2),)


# %timeit pd.read_csv("2021VAERSDATA.csv.gz", encoding="iso-8859-1")
# %timeit csv.read_csv("2021VAERSDATA.csv.gz")


# REMOVE SYMPTOM_TEXT


vdata_pd = pd.read_csv("2021VAERSDATA.csv.gz", encoding="iso-8859-1", usecols=lambda x: x != "SYMPTOM_TEXT")
data_pd.info(memory_usage="deep")

#columns.remove("SYMPTOM_TEXT")
vdata_arrow = csv.read_csv(
    "2021VAERSDATA.csv.gz",
     convert_options=csv.ConvertOptions(include_columns=columns))
vdata_arrow.nbytes

# %timeit pd.read_csv("2021VAERSDATA.csv.gz", encoding="iso-8859-1", usecols=lambda x: x != "SYMPTOM_TEXT")
# %timeit csv.read_csv("2021VAERSDATA.csv.gz", convert_options=csv.ConvertOptions(include_columns=columns))

vdata = vdata_arrow.to_pandas()
vdata.info(memory_usage="deep")



# Theres more
vdata = vdata_arrow.to_pandas(self_destruct=True)


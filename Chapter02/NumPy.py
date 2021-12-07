import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

vdata = pd.read_csv(
    "2021VAERSDATA.csv.gz", encoding="iso-8859-1")

vdata["STATE"] = vdata["STATE"].str.upper()
top_states = pd.DataFrame({
    "size": vdata.groupby("STATE").size().sort_values(ascending=False).head(5)}).reset_index()
top_states["rank"] = top_states.index
top_states = top_states.set_index("STATE")
top_vdata = vdata[vdata["STATE"].isin(top_states.index)]
top_vdata["state_code"] = top_vdata["STATE"].apply(
    lambda state: top_states["rank"].at[state]
).astype(np.uint8)
top_vdata = top_vdata[top_vdata["AGE_YRS"].notna()]
top_vdata.loc[:,"AGE_YRS"] = top_vdata["AGE_YRS"].astype(int)
top_states

age_state = top_vdata[["state_code", "AGE_YRS"]]
age_state["state_code"]
state_code_arr = age_state["state_code"].values
type(state_code_arr), state_code_arr.shape, state_code_arr.dtype

age_state["AGE_YRS"]
age_arr = age_state["AGE_YRS"].values
type(age_arr), age_arr.shape, age_arr.dtype

age_arr.max()

age_state_mat = np.zeros((5,6), dtype=np.uint64)
for row in age_state.itertuples():
    age_state_mat[row.state_code, row.AGE_YRS//20] += 1
age_state_mat

cal = age_state_mat[0,:]
kids = age_state_mat[:,0]

def compute_frac(arr_1d):
    return arr_1d / arr_1d.sum()

frac_age_stat_mat = np.apply_along_axis(compute_frac, 1, age_state_mat)

perc_age_stat_mat = frac_age_stat_mat * 100
perc_age_stat_mat = perc_age_stat_mat.astype(np.uint8)
perc_age_stat_mat

perc_age_stat_mat = perc_age_stat_mat[:, :5]
perc_age_stat_mat


fig = plt.figure()
ax = fig.add_subplot()
ax.matshow(perc_age_stat_mat, cmap=plt.get_cmap("Greys"))
ax.set_yticks(range(5))
ax.set_yticklabels(top_states.index)
ax.set_xticks(range(6))
ax.set_xticklabels(["0-19", "20-39", "40-59", "60-79", "80-99", "100-119"])
fig.savefig("matrix.png")



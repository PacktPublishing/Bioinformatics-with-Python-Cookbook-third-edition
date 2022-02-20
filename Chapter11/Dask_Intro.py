# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import zarr

mosquito = zarr.open('data/AG1000G-AO/2L/calldata/GT')
mosquito
zarr.array(mosquito, chunks=(48525747 // 4, 81, 2), store='data/rechunk')

# -

a = zarr.open('data/rechunk')
dir(a)
#a.chunks

# +
import dask.array as da

mosquito = da.from_zarr('data/AG1000G-AO/2L/calldata/GT')
# ^^^ load array
# -

mosquito

rechunked_mosquito = da.rechunk(mosquito, chunks=(300000, -1, -1))
rechunked_mosquito

chunk_pos_size = gt_2l.chunks[0]
max_pos = gt_2l.shape[0]

complete_data = 0
more_anc_hom = 0
total_pos = 0
# do a single case
for chunk_pos in range(ceil(max_pos / chunk_pos_size)):
    start_pos = chunk_pos * chunk_pos_size
    end_pos = min(max_pos + 1, (chunk_pos + 1) * chunk_pos_size)
    my_chunk = np.array(
        gt_2l[start_pos:end_pos, :, :])
    #print(start_pos, end_pos, my_chunk.shape)
    num_samples = my_chunk.shape[1]
    num_miss, num_anc_hom, num_het = calc_stats(my_chunk)
    chunk_complete_data = np.sum(np.equal(num_miss, 0))
    #print(end_pos - start_pos, my_chunk.shape, num_anc_hom.shape, num_het.shape)
    chunk_more_anc_hom = np.sum(num_anc_hom > num_het)
    complete_data += chunk_complete_data
    more_anc_hom += chunk_more_anc_hom
    total_pos += (end_pos - start_pos)
print(complete_data, more_anc_hom, total_pos)


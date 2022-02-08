# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
### XXX This is probably to remove
# -

sl_2014 = []
drc_2007 = []
for seq in ebola_seqs.taxon_set:
    if seq.label.startswith('EBOV_2014') and len(sl_2014) < 8:
        sl_2014.append(ebola_seqs[seq])
    elif seq.label.startswith('EBOV_2007'):
        drc_2007.append(ebola_seqs[seq])

print(len(sl_2014))
print(len(drc_2007))

pair_stats = popgenstat.PopulationPairSummaryStatistics(sl_2014, drc_2007)

print('Average number of pairwise differences (total): %s' %
      pair_stats.average_number_of_pairwise_differences)
print('Average number of pairwise differences between populations: %s' %
      pair_stats.average_number_of_pairwise_differences_between)
print('Average number of pairwise differences within populations: %s' %
      pair_stats.average_number_of_pairwise_differences_within)
print('Average number of new pairwise differences : %s' %
      pair_stats.average_number_of_pairwise_differences_net)
print('Number of segregating sites: %s' %
      pair_stats.num_segregating_sites)
print("Watterson's theta: %s" %
      pair_stats.wattersons_theta)
print("Wakeley's Psi: %s" % pair_stats.wakeleys_psi)
print("Tajima's D: %s" % pair_stats.tajimas_d)

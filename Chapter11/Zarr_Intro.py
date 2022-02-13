# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# !mkdir -p data/AG1000G-BF-A/
# !gsutil -m rsync -r \
#         -x '.*/calldata/(AD|GQ|MQ)/.*' \
#         gs://vo_agam_release/v3/snp_genotypes/all/AG1000G-BF-A/ \
#         data/AG1000G-BF-A/

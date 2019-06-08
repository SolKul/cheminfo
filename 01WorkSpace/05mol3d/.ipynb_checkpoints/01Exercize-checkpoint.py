# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from rdkit import rdBase
print('rdkit version: {}'.format(rdBase.rdkitVersion))

import pubchempy as pcp

pcp.download('SDF','./4133.sdf',4133,'cid',record_type='3d',overwrite=True)

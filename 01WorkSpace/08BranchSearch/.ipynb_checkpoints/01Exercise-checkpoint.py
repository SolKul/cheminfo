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

from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw

r32=Chem.AddHs(Chem.MolFromSmiles('C(F)F'))
display(Draw.MolToImage(r32))
print(Chem.MolToMolBlock(r32))



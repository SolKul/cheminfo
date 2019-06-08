# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.2
#   kernelspec:
#     display_name: Python cheminfo
#     language: python
#     name: cheminfo
# ---

from rdkit import rdBase
print('rdkit version: {}'.format(rdBase.rdkitVersion))

# +
import pubchempy as pcp
from rdkit import Chem
 
benzene = pcp.get_compounds('benzene', 'name')
if len(benzene) == 1: benzene = benzene[0]
smiles = benzene.canonical_smiles
print(smiles) # 'C1=CC=CC=C1'
mol_ben = Chem.MolFromSmiles(smiles)
print(type(mol_ben)) # rdkit.Chem.rdchem.Mol
print(Chem.MolToMolBlock(mol_ben))

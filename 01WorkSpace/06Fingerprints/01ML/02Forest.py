# -*- coding: utf-8 -*-
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

# +
import os
import os.path
import sys

p_mod=os.path.join(os.environ['HOME'],'notebooks/99MyModules')
sys.path.append(p_mod)

## 1. ライブラリのインポート
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw, PandasTools, Descriptors
 
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
print(rdBase.rdkitVersion) # 2019.03.2
# -

import descarray as da

df_AMES=pd.read_csv('./ci900161g_si_001/smiles_cas_N6512.smi',sep='\t', header=None)
df_AMES.columns = ['smiles', 'CAS_NO', 'activity']
da.descDf(df_AMES)

PandasTools.AddMoleculeColumnToFrame(frame=df_AMES, smilesCol='smiles')
## 4. 読み込めない分子の削除
df_AMES_d=df_AMES.dropna()
df_AMES_d=df_AMES_d.assign(mw=df_AMES_d['ROMol'].map(Descriptors.MolWt))#.map(Descriptors.MolWt)
da.descDf(df_AMES_d)

print('全部で%d'%len(df_AMES_d)) # 6506
n_neg=np.sum(df_AMES_d.activity == 0)
n_pos=np.sum(df_AMES_d.activity == 1)
print('陰性は{}、陽性は{}'.format(n_neg,n_pos))

maccskeys = []
for mol in df_AMES_d.ROMol:
    maccskey = [i for i in AllChem.GetMACCSKeysFingerprint(mol)]
    maccskeys.append(maccskey)
maccskeys = np.array(maccskeys)
print(maccskeys.shape) # (6506, 167)
## 出力ラベル
target = df_AMES_d.activity
print(target.shape) # (6506,)

# +
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(maccskeys, df_AMES_d.activity, random_state=0)
 
tree = DecisionTreeClassifier(random_state=0)
tree.fit(X_train, y_train)
print('accuracy on train set: {:.3f}'.format(tree.score(X_train, y_train)))
print('accuracy on test set: {:.3f}'.format(tree.score(X_test, y_test)))

# +
accs_train = []
accs_test = []
 
for i in range(1,21):
    tree_i = DecisionTreeClassifier(max_depth=i, random_state=0)
    tree_i.fit(X_train, y_train)
    acc_train = tree_i.score(X_train, y_train)
    acc_test = tree_i.score(X_test, y_test)
    accs_train.append(acc_train)
    accs_test.append(acc_test)
# -

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.plot(accs_train)
ax.plot(accs_test)

tree2 = DecisionTreeClassifier(max_depth=2)
tree2.fit(X_train, y_train)

from sklearn.tree import export_graphviz
export_graphviz(tree2, out_file='tree2.dot', class_names=['negative', 'positive'],
                feature_names=range(1,168), filled=True)

import graphviz
with open('tree2.dot') as f:
    dot_graph = f.read()
graphviz.Source(dot_graph)

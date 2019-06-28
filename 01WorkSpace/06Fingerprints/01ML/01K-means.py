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
## 1. ライブラリのインポート
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw, PandasTools, Descriptors
 
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
print(rdBase.rdkitVersion) # 2019.03.2

# +
import os
import os.path
import sys

p_mod=os.path.join(os.environ['HOME'],'notebooks/99MyModules')
sys.path.append(p_mod)

import descarray as da

DataDir=os.path.join(os.environ['HOME'],'notebooks/50Data/')
# -

AmesF=DataDir+'ci900161g_si_001/smiles_cas_N6512.smi'
df_AMES=pd.read_csv(AmesF,sep='\t', header=None)
df_AMES.columns = ['smiles', 'CAS_NO', 'activity']
da.descDf(df_AMES)

PandasTools.AddMoleculeColumnToFrame(frame=df_AMES, smilesCol='smiles')
da.descDf(df_AMES)

## 4. 読み込めない分子の削除
df_AMES_d=df_AMES.dropna()
df_AMES_d=df_AMES_d.assign(mw=df_AMES_d['ROMol'].map(Descriptors.MolWt))#.map(Descriptors.MolWt)
da.descDf(df_AMES_d)

df_AMES_d.loc[:,['ROMol']].iloc[0]

print('全部で%d'%len(df_AMES_d)) # 6506
n_neg=np.sum(df_AMES_d.activity == 0)
n_pos=np.sum(df_AMES_d.activity == 1)
print('陰性は{}、陽性は{}'.format(n_neg,n_pos))

## 5. 分子量分布をプロット
sns.violinplot(x='activity', y='mw', data=df_AMES_d)

## モルガンフィンガープリント
morgan_fp = []
for mol in df_AMES_d.ROMol:
    fp = [i for i in AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)]
    morgan_fp.append(fp)
morgan_fp = np.array(morgan_fp)
print(morgan_fp.shape) # (6506, 2048)
## MACCS Keys
maccskeys = []
for mol in df_AMES_d.ROMol:
    maccskey = [i for i in AllChem.GetMACCSKeysFingerprint(mol)]
    maccskeys.append(maccskey)
maccskeys = np.array(maccskeys)
print(maccskeys.shape) # (6506, 167)
## 出力ラベル
target = df_AMES_d.activity
print(target.shape) # (6506,)

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(morgan_fp, target, random_state=0)
print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)

from sklearn.neighbors import KNeighborsClassifier
morgan_train_acc = []
morgan_test_acc = []
for i in range(1,11):
    knn = KNeighborsClassifier(n_neighbors=i)
    knn.fit(X_train, y_train)
    train_acc = knn.score(X_train, y_train)
    test_acc = knn.score(X_test, y_test)
    print('test accuracy with k={}: {:.3f}'.format(i, test_acc))
    morgan_train_acc.append(train_acc)
    morgan_test_acc.append(test_acc)

import threading
from sklearn.neighbors import KNeighborsClassifier


def calcknn(i):
    knn = KNeighborsClassifier(n_neighbors=i)
    knn.fit(X_train, y_train)
    train_acc = knn.score(X_train, y_train)
    test_acc = knn.score(X_test, y_test)
    print('test accuracy with k={}: {:.3f}'.format(i, test_acc))
    morgan_train_acc.append(train_acc)
    morgan_test_acc.append(test_acc)


thread_1 = threading.Thread(target=calcknn1)
thread_2 = threading.Thread(target=calcknn2)

thread_1.start()
thread_2.start()

morgan_fp.shape

import concurrent.futures

futures=[]
executor = concurrent.futures.ProcessPoolExecutor(max_workers=8)
for i in range(5,8):
    futures.append(executor.submit(calcknn,i))

futures[7].running()

f=concurrent.futures.Future()

executor.

executor.done()

pow(2,3)

np.sum(df_AMES_d.isnull())

CSO2C = Chem.MolFromSmiles('CS(=O)(=O)C')
Draw.MolToImage(Chem.MolFromSmiles('NNC(=O)CNC(=O)C=N#N'))

df_AMES[df_AMES.ROMol >= CSO2C]

Draw.MolToImage(CSO2C)

type(df_AMES.ROMol[0])

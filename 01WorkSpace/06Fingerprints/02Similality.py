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

# # RDKitでフィンガープリントを使った分子類似性の判定

import pandas as pd
from rdkit import rdBase, Chem, DataStructs
print(rdBase.rdkitVersion) # 2017.09.1
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem,Draw
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs, Torsions

suppl = Chem.SDMolSupplier('../サリチル酸メチル.sdf',removeHs=False)
mols=[]
for x in suppl:
    if x is not None:
        mols.append(x)
len(mols)

# ## MACCS Keys
# ケモインフォマティクスでは非常に有名なMDL社の開発した化学構造データベースに由来するフィンガープリントです．
#
# 全部で166の部分構造についての有無を調べ上げたもので，RDKit内の情報保持のために1ビット使っているため全部で167ビットのフィンガープリントになります．部分構造を有する場合は1が，ない場合は0が格納されています．

maccs_fps = [AllChem.GetMACCSKeysFingerprint(mol) for mol in mols]
maccs = DataStructs.BulkTanimotoSimilarity(maccs_fps[0], maccs_fps[1:])

# ## Topologicalフィンガープリント (RDKitフィンガープリント)  
# ハッシュ化フィンガープリントの中では有名なDaylightフィンガープリントに類似したものと考えてよいです．RDKitフィンガープリントとも呼ばれています．
#
# 一定の結合数に相当する原子と結合種類を格納する方法で，事前に部分構造を用意する必要がないことから，より柔軟に分子構造を表現することが可能です．

rdkit_fps = [Chem.Fingerprints.FingerprintMols.FingerprintMol(mol) for mol in mols]
rdkit = DataStructs.BulkTanimotoSimilarity(rdkit_fps[0], rdkit_fps[1:])

# ## Morganフィンガープリント (Circularフィンガープリント) 
# 原子からある距離にある部分構造を数え上げていく，「circular substructures」と呼ばれるタイプのフィンガープリントです．radiusの値で距離を設定します．いわゆるECFP（Extended Connectivity Fingerprint）フィンガープリントと類似のものですが，探索距離の定義が異なる点に注意が必要です．
#
# よく使われるECFP4はradius=2の設定（デフォルト）と大体同じとのことです．FCFP（Functional Connectivity Fingerprint）様のフィンガープリントを用いたい場合にはuseFeatures=Trueを設定します

# ![グラフ](./01ECFP.png)

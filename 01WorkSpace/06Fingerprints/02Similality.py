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

import sys
print(sys.version)

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

morgan_fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
morgan = DataStructs.BulkTanimotoSimilarity(morgan_fps[0], morgan_fps[1:])

# ## Avalonフィンガープリント  
# RDKitから用いることのできるAvalon Chemoinformatics TookKitにもフィンガープリントが実装されています．

avalon_fps = [pyAvalonTools.GetAvalonFP(mol) for mol in mols]
avalon = DataStructs.BulkTanimotoSimilarity(avalon_fps[0], avalon_fps[1:])

# ## アトムペアと二面角  
# アトムペアを用いるフィンガープリントは全ての原子同士についてその最小結合パスを記録します．各原子には原子の種類，結合する原子のうち水素以外の元素の数，π電子数の情報も記録されます．
#
# 二面角を用いるフィンガープリントは，二面角を形成する4連続原子について，やはり原子の種類，結合する原子のうち水素以外の元素の数，π電子数を記録します．

# ## 各フィンガープリントの比較

df = pd.DataFrame({'RDKit': rdkit,
                  'Avalon': avalon,
                  'MACCS': maccs,
                  'Morgan': morgan})
df.corr().round(2)

morgan_default = AllChem.GetMorganFingerprint(mols[0], 2)
morgan_bitvect = AllChem.GetMorganFingerprintAsBitVect(mols[0], 2, 2048)
print(morgan_default.GetLength(), morgan_bitvect.GetNumBits())

from IPython.display import SVG

### Morganフィンガープリント
bitI_morgan = {}
fp_morgan = AllChem.GetMorganFingerprintAsBitVect(mols[0], 2, bitInfo=bitI_morgan)
### RDKitフィンガープリント
bitI_rdkit = {}
fp_rdkit = Chem.RDKFingerprint(mols[0], bitInfo=bitI_rdkit)

for row in range(10):
    print(row)

bitI_rdkit.keys()

SVG(Draw.DrawRDKitBit(mols[0],8,bitI_rdkit,legend='bit: '+str(bit)))

for bit in bitI_rdkit .keys():
    display(SVG(Draw.DrawRDKitBit(mols[0],bit,bitI_rdkit,legend='bit: '+str(bit))))

morgan_turples = ((mols[0], bit, bitI_morgan) for bit in list(bitI_morgan.keys())[:12])
SVG(Draw.DrawMorganBits(morgan_turples, molsPerRow=4, legends=['bit: '+str(x) for x in list(bitI_morgan.keys())[:12]]))

SVG(Draw.DrawRDKitBits([(mols[0],52,bitI_rdkit)]))

SVG(Draw.DrawRDKitBit(mols[0],52,bitI_rdkit))

list(bitI_rdkit.keys())[8]

rdkit_turples =((mols[0], bit, bitI_rdkit) for bit in list(bitI_rdkit.keys())[7:9])

list(bitI_rdkit.keys())[7:9]

import copy

m=copy.copy(mols[0])

id(mols[0])

id(m)

Chem.Kekulize(m)

rdkit_turples = [(m,52,bitI_rdkit)]

Chem.Draw.MolToImage(mol)

mol = Chem.MolFromSmiles('c1cncn(C)1')
bi = {}
_ = AllChem.RDKFingerprint(mol, minPath=1, maxPath=5, bitInfo=bi)
tuples = [(mol, k, bi) for k in bi.keys()]
SVG(Draw.DrawRDKitBit(mol,1505,bi))

bi[1505]

SVG(Draw.DrawRDKitBits([(mol,1505,bi)]))





tuples[38]

e_cause=1505

SVG(Draw.DrawRDKitBits([(mol,1505,bi)]))

SVG(Draw.DrawRDKitBit(mol,1505,bi))

import inspect

print(inspect.getsource(Draw.DrawRDKitBit))

print(inspect.getsource(Draw.DrawRDKitBits))

bi[e_cause]

bondpath=bi[e_cause][0]

Draw.DrawRDKitEnvs([(mol,bondpath)])

print(inspect.getsource(Draw.DrawRDKitEnvs))


def DrawRDKitEnvs(envs, molsPerRow=3, subImgSize=(150, 150), baseRad=0.3, useSVG=True,
                  aromaticColor=(0.9, 0.9, 0.2), extraColor=(0.9, 0.9,0.9),
                  nonAromaticColor=None, legends=None, **kwargs):
    submols = []
    highlightAtoms = []
    atomColors = []
    highlightBonds = []
    bondColors = []
    highlightRadii = []
    for mol, bondpath in envs:
        menv = _getRDKitEnv(mol, bondpath, baseRad, aromaticColor, extraColor, nonAromaticColor,
                            **kwargs)
        submols.append(menv.submol)
        highlightAtoms.append(menv.highlightAtoms)
        atomColors.append(menv.atomColors)
        highlightBonds.append(menv.highlightBonds)
        bondColors.append(menv.bondColors)
        highlightRadii.append(menv.highlightRadii)

    if legends is None:
        legends = [''] * len(envs)

    nRows = len(envs) // molsPerRow
    if len(envs) % molsPerRow:
        nRows += 1

    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
    # Drawing
    if useSVG:
    drawer = rdMolDraw2D.MolDraw2DSVG(fullSize[0], fullSize[1], subImgSize[0], subImgSize[1])
    else:
    drawer = rdMolDraw2D.MolDraw2DCairo(fullSize[0], fullSize[1], subImgSize[0], subImgSize[1])

    drawopt = drawer.drawOptions()
    drawopt.continuousHighlight = False
    drawer.DrawMolecules(submols, legends=legends, highlightAtoms=highlightAtoms,
                       highlightAtomColors=atomColors, highlightBonds=highlightBonds,
                       highlightBondColors=bondColors, highlightAtomRadii=highlightRadii, **kwargs)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


baseRad=0.3
useSVG=True
aromaticColor=(0.9, 0.9, 0.2)
extraColor=(0.9, 0.9,0.9)
nonAromaticColor=None
menv = Draw._getRDKitEnv(mol, bondPath, baseRad, aromaticColor, extraColor, nonAromaticColor)

menv

print(inspect.getsource(Draw.rdMolDraw2D.MolDraw2DSVG))

Draw.DrawRDKitEnv

inspect.getfile(Draw.rdMolDraw2D)

import os


# + {"language": "bash"}
# cp -f /opt/conda/lib/python3.7/site-packages/rdkit/Chem/Draw/rdMolDraw2D.so ./Draw.rdMolDraw2D.so

# + {"language": "bash"}
# ls /opt/conda/lib/python3.7/site-packages/rdkit/Chem/Draw/
# -

def DrawRDKitEnv(mol, bondPath, molSize=(150, 150), baseRad=0.3, useSVG=True,
                 aromaticColor=(0.9, 0.9, 0.2), extraColor=(0.9, 0.9,0.9),
                 nonAromaticColor=None, **kwargs):
    menv = _getRDKitEnv(mol, bondPath, baseRad, aromaticColor, extraColor, nonAromaticColor, **kwargs)

    # Drawing
    if useSVG:
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    else:
        drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])

    drawopt = drawer.drawOptions()
    drawopt.continuousHighlight = False

    drawer.DrawMolecule(menv.submol, highlightAtoms=menv.highlightAtoms,
                      highlightAtomColors=menv.atomColors, highlightBonds=menv.highlightBonds,
                      highlightBondColors=menv.bondColors, highlightAtomRadii=menv.highlightRadii,
                      **kwargs)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


mol = Chem.MolFromSmiles('c1cncn(C)1')
bondPath = [0, 1, 2, 3, 5]
submol = Chem.PathToSubmol(mol, bondPath)
# submol.UpdatePropertyCache()
# submol.Debug()
Chem.Kekulize(submol)

Chem.MolToSmiles(submol)

# rdkit_turples =((mols[0], bit, bitI_rdkit) for bit in list(bitI_rdkit.keys())[8:9])
SVG(Draw.DrawRDKitBits(rdkit_turples))

Draw.DrawMorganBit(mol[0],0,)

rdkit_fps[0]

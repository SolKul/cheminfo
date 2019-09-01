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

from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw

r32=Chem.AddHs(Chem.MolFromSmiles('C(F)F'))
display(Draw.MolToImage(r32))
print(Chem.MolToMolBlock(r32))

MthRad=Chem.AddHs(Chem.MolFromSmiles('[CH3]'))
display(Draw.MolToImage(MthRad))
print(Chem.MolToMolBlock(MthRad))

import itertools
import numpy as np

l_mb=[]
l_b=[[0,0,1],[0,1,0],[1,0,0]]
l_hb=list(itertools.combinations_with_replacement(l_b, 3))
for j in range(len(l_hb)):
    l_hb_p=[]
    for i in range(3):
        l_1hb=l_hb[j][i].copy()
        l_1hb.insert(i+1,0)
        l_hb_p.append(l_1hb)
    m_bh=np.array(l_hb_p)
    v_b_c=np.insert(m_bh[:,0],0,0)
    m_b=np.concatenate([v_b_c[np.newaxis,:],m_bh],axis=0)
    if(np.all(m_b==m_b.T)):
        m_b[0,0]=4-m_b[0,:].sum()
        l_mb.append(m_b)
l_mb

for i in range(4):
    for j in range(4):
        if i>j:
            if(m_b[i,j]>0):
                print(i,j,m_b[i,j])


class Compound:
    def __init__(self,i,j,l_o_atm,m_o_b):
        self.l_atm=[i,j]
        #i,ｊを新しく化合物に追加
        self.l_o_atm=l_o_atm.copy()
        #もとの原子リストをコピー
        self.d_atm={i:l_o_atm[i],j:l_o_atm[j]}
        #化合物内の原子辞書をつくる
        self.m_o_b=m_o_b
        #もとの隣接行列をコピー

    def addatm(self,i):
        self.l_atm.append(i)
        self.d_atm[i]=self.l_o_atm[i]
        #iを新しく化合物に追加

    def exatms(self,l_t_atm):
        self.l_atm.extend(l_t_atm)
        for i in l_t_atm:
            self.d_atm[i]=self.l_o_atm[i]
        #l_t_atmを新しく化合物に追加

    def calcAd(self):
        #化合物に属する原子を全部まとめ終わってから
        #隣接行列、隣接リストを算出
        l_adm=[]
        self.l_n_atm=[]
        for i in self.l_atm:
            l_adm.append(self.m_o_b[i,self.l_atm])
            self.l_n_atm.append(self.d_atm[i])
            #原子のリストも作る。
            
        self.m_adm=np.array(l_adm)
        
        self.l_adl=[]
        for i in range(self.m_adm.shape[0]):
            for j in range(i+1,self.m_adm.shape[1]):
                if(self.m_adm[i,j]==0):
                    continue
                self.l_adl.append([i+1,j+1,self.m_adm[i,j]])


# +
l_p_atm=['C','C','H','H','H','H']
m_b=np.array([[0,3,1,0,0,0],[3,0,0,1,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,1],[0,0,0,0,0,0]])
for i in range(m_b.shape[0]):
    for j in range(m_b.shape[1]):
        if i>j:
            m_b[i,j]=m_b[j,i]
print(m_b)
                
d_atom = {}
#どの原子がどの原子群に属しているかという辞書
d_comp={}
#原子群ごとにどの原子が属しているかという辞書
d_lt={}
comp_c=0
for i in range(m_b.shape[0]):
    for j in range(i+1,m_b.shape[1]):
        if(m_b[i,j]==0):
            continue
        if i not in d_atom.keys():
            if j not in d_atom.keys():
                d_atom[i]=comp_c
                d_atom[j]=comp_c
                d_comp[comp_c]=Compound(i,j,l_p_atm,m_b)
                comp_c+=1
                #i,ｊを新しく原子群に追加
            else:
                d_atom[i]=d_atom[j]
                d_comp[d_atom[j]].addatm(i)
                #iをｊの属している原子群に
        else:
            if j not in d_atom.keys():
                d_atom[j]=d_atom[i]
                d_comp[d_atom[i]].addatm(j)
                #iをｊの属している原子群に
            else:
                if(d_atom[i]<d_atom[j]):
                    Comp=d_comp.pop(d_atom[j])
                    #つながっているとわかった原子群をpopで取り出す
                    d_comp[d_atom[i]].exatms(Comp.l_atm)
                    #番号の若い方に追加
                    for atm_num in Comp.l_atm:
                        d_atom[atm_num]=d_atom[i]
                        #つながっているわかった原子群の属する原子群を若い方に更新
                elif(d_atom[i]>d_atom[j]):
                    Comp=d_comp.pop(d_atom[i])
                    #つながっているとわかった原子群をpopで取り出す
                    d_comp[d_atom[j]].exatms(Comp.l_atm)
                    #番号の若い方に追加
                    for atm_num in Comp.l_atm:
                        d_atom[atm_num]=d_atom[j]

for k in d_comp:
    d_comp[k].calcAd()
    c=d_comp[k]
    print(k,c.d_atm,c.l_atm,'\n',c.l_adl,c.l_n_atm)
# -

d_comp[0].calcAdL()
d_comp[0].m_b_adl

# +
l_p_atm=['C','C','H','H','H','H']
m_b=np.array([[0,3,1,0,0,0],[3,0,0,1,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
for i in range(m_b.shape[0]):
    for j in range(m_b.shape[1]):
        if i>j:
            m_b[i,j]=m_b[j,i]
print(m_b)
                
d_atom = {}
#どの原子がどの原子群に属しているかという辞書
d_comp={}
#原子群ごとにどの原子が属しているかという辞書
d_lt={}
comp_c=0
for i in range(m_b.shape[0]):
    for j in range(i+1,m_b.shape[1]):
        if(m_b[i,j]==0):
            continue
        if i not in d_atom.keys():
            if j not in d_atom.keys():
                d_atom[i]=comp_c
                d_atom[j]=comp_c
                d_comp[comp_c]=[i,j]
                comp_c+=1
                #i,ｊを新しく原子群に追加
            else:
                d_atom[i]=d_atom[j]
                d_comp[d_atom[j]].append(i)
                #iをｊの属している原子群に
        else:
            if j not in d_atom.keys():
                d_atom[j]=d_atom[i]
                d_comp[d_atom[i]].append(j)
            else:
                if(d_atom[i]<d_atom[j]):
                    l_atm_t=d_comp.pop(d_atom[j])
                    #つながっているとわかった原子群をpopで取り出す
                    d_comp[d_atom[i]].extend(l_atm_t)
                    #番号の若い方に追加
                    for atm_num in l_atm_t:
                        d_atom[atm_num]=d_atom[i]
                        #つながっているわかった原子群の属する原子群を若い方に更新
                elif(d_atom[i]>d_atom[j]):
                    l_atm_t=d_comp.pop(d_atom[i])
                    for atm_num in l_atm_t:
                        d_atom[atm_num]=d_atom[j]
                    d_comp[d_atom[j]].extend(l_atm_t)
#         if j not in compounds.keys():
#             compounds[j]=[comp_c]
print(d_atom)
print(d_comp)
# -

l_c_atm=[]
l_m_c=[]
for k in d_comp:
    m_cb=[]
    for atm in d_comp[k]:
        m_cb.append(m_b[atm,d_comp[k]])
    l_c_atm.append(list(itemgetter(*d_comp[k])(l_p_atm)))
    l_m_c.append(np.array(m_cb))
print(l_c_atm)
print(l_m_c)

l_l_b_map=[]
for m_c in l_m_c:
    l_b_map=[]
    for i in range(m_c.shape[0]):
        for j in range(i+1,m_c.shape[1]):
            if(m_c[i,j]==0):
                continue
            l_b_map.append([i,j,m_c[i,j]])
    l_l_b_map.append(l_b_map)
l_l_b_map

for i in range(len(l_c_atm)):
    print(' '*2+str(len(l_c_atm[i]))+' '*2+str(len(l_l_b_map[i]))+'  0  0  0  0  0  0  0  0999 V2000')
    for atm in l_c_atm[i]:
        print('    0.0   '*3+' '+atm+'   0'+'  0'*11)
    for l_b_map in l_l_b_map[i]:
        print(' '*2+str(l_b_map[0]+1)+' '*2+str(l_b_map[1]+1)+' '*2+str(l_b_map[2])+'  0')

s_test='''
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0       0.0       0.0    C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0       0.0       0.0    C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0       0.0       0.0    H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0       0.0       0.0    H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  3  0
  1  3  1  0
  2  4  1  0
M  END'''

test_m=Chem.MolFromMolBlock(s_test,removeHs=False)
AllChem.Compute2DCoords(test_m)
Draw.MolToImage(test_m)

len(l_l_b_map[i])

# +
#aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
#xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
#111222tttsssxxxrrrccc

# +
#  5  4  0  0  0  0  0  0  0  0999 V2000
#    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#  1  2  1  0
# -



  1  2  1  0



from operator import itemgetter
list(itemgetter(*d_comp[k])(l_p_atm))

m_cs

l_1hb.extend([2,30,3])

l_1hb

m_b.shape[0]



m_b

{'Utsunomiya': ['Mito', 'Maebashi', 'Saitama', 'Chiba', 'Tokyo'], 'Wakayama': ['Kyoto', 'Osaka', 'Kobe', 'Nara', 'Tokushima'], 'Nagasaki': ['Saga', 'Kumamoto'], 'Tokyo': ['Mito', 'Utsunomiya', 'Maebashi', 'Saitama', 'Chiba', 'Yokohama'], 'Okayama': ['Tottori', 'Tokushima', 'Takamatsu'], 'Miyazaki': ['Kagoshima'], 'Nagano': ['Maebashi', 'Toyama'], 'Toyama': ['Kanazawa', 'Nagano'], 'Maebashi': ['Utsunomiya', 'Saitama', 'Tokyo', 'Kofu', 'Nagano'], 'Nagoya': ['Gifu', 'Tsu'], 'Chiba': ['Mito', 'Utsunomiya', 'Saitama', 'Tokyo', 'Yokohama'], 'Sendai': ['Yamagata', 'Fukushima'], 'Tokushima': ['Kobe', 'Wakayama', 'Okayama', 'Takamatsu'], 'Kofu': ['Maebashi', 'Shizuoka'], 'Oita': ['Yamaguchi', 'Kumamoto'], 'Kobe': ['Otsu', 'Kyoto', 'Osaka', 'Nara', 'Wakayama', 'Tokushima'], 'Otsu': ['Gifu', 'Tsu', 'Kyoto', 'Osaka', 'Kobe', 'Nara'], 'Tsu': ['Gifu', 'Nagoya', 'Otsu', 'Kyoto', 'Osaka', 'Nara'], 'Nara': ['Tsu', 'Otsu', 'Kyoto', 'Osaka', 'Kobe', 'Wakayama'], 'Hiroshima': ['Matsuyama'], 'Matsuyama': ['Hiroshima', 'Kochi'], 'Fukui': ['Kanazawa', 'Gifu'], 'Kagoshima': ['Miyazaki'], 'Saitama': ['Mito', 'Utsunomiya', 'Maebashi', 'Chiba', 'Tokyo', 'Yokohama'], 'Fukuoka': ['Saga', 'Kumamoto'], 'Takamatsu': ['Okayama', 'Tokushima', 'Kochi'], 'Saga': ['Fukuoka', 'Nagasaki', 'Kumamoto'], 'Kyoto': ['Tsu', 'Otsu', 'Osaka', 'Kobe', 'Nara', 'Wakayama'], 'Osaka': ['Tsu', 'Otsu', 'Kyoto', 'Kobe', 'Nara', 'Wakayama'], 'Gifu': ['Fukui', 'Nagoya', 'Tsu', 'Otsu'], 'Shizuoka': ['Kofu'], 'Kumamoto': ['Fukuoka', 'Saga', 'Nagasaki', 'Oita'], 'Mito': ['Utsunomiya', 'Saitama', 'Chiba', 'Tokyo'], 'Yokohama': ['Saitama', 'Chiba', 'Tokyo'], 'Fukushima': ['Sendai', 'Yamagata'], 'Yamagata': ['Sendai', 'Fukushima'], 'Yamaguchi': ['Oita'], 'Kanazawa': ['Toyama', 'Fukui'], 'Tottori': ['Okayama'], 'Kochi': ['Takamatsu', 'Matsuyama']}


list(range(2,10,2))

m_b[1][m_b[1]>0]

m_b=np.array(l_hb_p)
v_b_c=np.insert(m_b[:,0],0,0)
np.concatenate([v_b_c[np.newaxis,:],m_b],axis=0)



l_hb[0]

l_hb[0][0]
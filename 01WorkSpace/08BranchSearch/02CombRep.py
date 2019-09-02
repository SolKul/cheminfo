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

import numpy as np
import itertools

# +
#分子内の原子数がmのとき、
#原子価nの結合の仕方は
#ラジカルを自身への結合と考えれば、
#m種類のものから重複して何回でも取り出せるときの
#n個選び取る重複組合せ
n_atm=3
n_valence=1
l_cr=list(itertools.combinations_with_replacement(list(range(n_atm)),n_valence))
print(l_cr)

for i in range(len(l_cr)):
    l_bond=[]
    for j in range(n_atm):
        l_bond.append(l_cr[i].count(j))
        #ある1種の数字が何回出てくるかを数える
    print(l_bond)
# -

l_c_atm=['C','H','H','H']
l_vlc=[4,1,1,1]
n_atm=len(l_c_atm)
#原子のリスト
l_bond_p=[[[0]*n_atm]]
#取りうる隣接行列の行
#最初の行は後から決まるので0で良い
for k in range(1,len(l_c_atm)):
    n_atm_r=k+1
    #注目すればいい原子の数
    #三角形の横方向の数+1
    #詳しくは足してnになる組み合わせの数の応用
    #足してn以下になる組み合わせの数を参照
    n_valence=l_vlc[k]
    l_cr=list(itertools.combinations_with_replacement(list(range(n_atm_r)),n_valence))
    l_bond_r_p=[]
    #隣接行列のある1行が取りうる行のリスト
    for i in range(len(l_cr)):
        l_bond=[0]*n_atm
        for j in range(n_atm_r-1):
            l_bond[j]=l_cr[i].count(j)
        l_bond_r_p.append(l_bond)
    l_bond_p.append(l_bond_r_p)
l_bond_p

ar_bond_p

ar_bond_p[i,:j+1,j]

ar_bond_p

np.array([[4,1,1,1]])


def Base_n_to_10(l_X,n):
    out = 0
    for i in range(1,len(l_X)+1):
        out += int(l_X[-i])*(n**(i-1))
    return out#int out


Base_n_to_10(ar_bond_cons[0,0,:],4)

# +
ar_bond_p=np.array(list(itertools.product(*l_bond_p)))
for i in range(ar_bond_p.shape[0]):
    for j in range(ar_bond_p.shape[1]):
        ar_bond_p[i,:j,j]=ar_bond_p[i,j,:j]
        #対称行列にする
ar_b_cons=ar_bond_p.sum(axis=1)<=np.array([[4,1,1,1]])
#結合数が原子価を以下のものを抽出
ar_b_ind=np.all(ar_b_cons,axis=1)
ar_bond_cons=ar_bond_p[ar_b_ind]

# for i in range(ar_bond_cons.shape[0]):
#     for j in range(ar_bond_cons.shape[1]):
#         ar_bond_cons[i,j,j]=l_vlc[j]-(ar_bond_cons[i,j,:j].sum()+ar_bond_cons[i,j,j+1:].sum())
# -

l_unique=[0,1,1,1]

ar_c_atm=np.array(l_c_atm)
ar_uni=np.unique(ar_c_atm)
print(ar_uni)
l_c_table=[]
l_cons_b=[]
for i in range(len(ar_base)):
    l_sorted=[]
    for atm in ar_uni:
        ar_uni_b=(ar_c_atm==atm)
        l_sorted.extend(np.sort(ar_base[i][ar_uni_b]))
    if l_sorted not in l_c_table:
        l_c_table.append(l_sorted)
        l_cons_b.append(i)

ar_bond_cons[l_cons_b]

l_c_table

l_sorted==[0,0,0,0]

np.sort(ar_base[0,:])

# +
for i in range(ar_bond_cons.shape[0]):
    for j in range(ar_bond_cons.shape[1]):
        ar_bond_cons[i,j,j:]=0
        
l_l_base=[]
for i in range(ar_bond_cons.shape[0]):
    l_base=[]
    for j in range(ar_bond_cons.shape[1]):
        l_base.append(Base_n_to_10(ar_bond_cons[i,j,:],l_vlc[j]+1))
    l_l_base.append(l_base)
ar_base=np.array(l_l_base)
# -

ar_base[0]

ar_bond_cons

ar_bond_cons

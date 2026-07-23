import matplotlib.pyplot as plt
import matplotlib as mpl
import multiprocess as mp
import seaborn as sns
import numpy as np
import pandas as pd
import h5py
import pickle
import scipy
from scipy import stats
import cooler
import cooltools
import bioframe
from cooltools import insulation
import os



files=['old.mm10.mapq_30.1000.mcool::/resolutions/100000']

files2=[]
for f in files:
    c=cooler.Cooler(f)
    p=c.pixels()[:]
    p=p.loc[abs(p.bin1_id-p.bin2_id)!=0]
    bins=c.bins()[:]
    bins=bins[['chrom','start','end']]
    pixels=p
    x=f.split('.')[0]
    name_pattern=x+'_100k.cool'
    cooler.create_cooler(name_pattern, bins, pixels)
    files2.append(name_pattern)

number=[]
for f in files2:
    c=cooler.Cooler(f)
    number.append(c.info['sum'])


numb=min(number)-1

for f in files2:
    cooltools.sample(f,'sampled_'+f,count=numb,exact=True)

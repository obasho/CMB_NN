import pysm3
import pysm3.units as u
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import os
import multiprocessing


def bin2nbit(binary,n):
  while(len(binary)<n):
    binary='0'+binary
  return binary
for i in range(12):
    os.mkdir(str(i))
    os.mkdir(str(i)+'o')



nside=1024
frequency=[85,95,145,155,220,270]
nside=1024
frequency=[85,95,145,155,220,270]
fwhm=[25.5,22.7,25.5,22.7,13,13]*np.pi/(180*60)
sensit=np.array([1.31,1.15,1.78,1.91,4.66,7.99])



def dataprep(b):
    r=-1
    for q in range(100):
        qb=q+100*b
        if r<0 or r>5:
            sky= pysm3.Sky(nside=1024, preset_strings=["s1", "f1", "a1", "d1"])
            sky_cmb=pysm3.Sky(nside=1024,preset_strings="c1")
        cmb=sky_cmb.get_emission(frequency[r]*u.GHz)
        map = sky.get_emission(frequency[r] * u.GHz)
        map = map.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(frequency[r]*u.GHz))
        cmb = cmb.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(frequency[r]*u.GHz))
        mapg=hp.sphtfunc.smoothing(map[0],fwhm[r])
        cmbg=hp.sphtfunc.smoothing(cmb[0],fwhm[r])
        fmap=mapg+cmbg+np.random.normal(scale=sensit[r],size=len(mapg))
        mp=np.zeros((12,1024,1024))
        m=fmap
        mpo=np.zeros((12,1024,1024))
        mo=cmb[0]
        for k in range(12):
            pixar='0'*20
            kb=bin2nbit(bin(k)[2:],4)
            pixar=np.array(list(kb+pixar))
            for i in range(1024):
                ib=np.array(list(bin2nbit(bin(i)[2:],10)))
                for index, value in np.ndenumerate(ib):
                    pixar[4+2*index[0]]=value
                for j in range(1024):
                    jb=np.array(list(bin2nbit(bin(j)[2:],10)))
                    for index, value in np.ndenumerate(jb):
                        pixar[5+2*index[0]]=value
                    PixNested=int(''.join([str(score) for score in pixar]),2)
                    PixRing=hp.nest2ring(1024,PixNested)
                    mp[k,i,j]=m[PixRing]
                    mpo[k,i,j]=mo[PixRing]
            plt.imsave(str(k)+'/'+str(qb))
            plt.imsave(str(k)+'o'+'/'+str(qb))



processes=[]
for b in range(16):
    p=multiprocessing.Process(target=dataprep,kwargs=[b])
    p.start()
    processes.append(p)
for process in processes:
    process.join()
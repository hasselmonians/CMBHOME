# ReadBehavior.py
"""Python module to import behavioral data from AVI"""
import pims
import scipy.ndimage
from joblib import Parallel, delayed
import multiprocessing
import time
import numpy as np
import sys
from scipy.ndimage.morphology import generate_binary_structure
from matplotlib import pyplot as plt


def GetMean(v):
    m = np.zeros_like(v[0],dtype='double')
    
    for i in range(0,v.get_metadata()['nframes'],100):
        m+=v[i]
        
    m=m/(v.get_metadata()['nframes']/100)
    return m
    
    
def ShowFrame(frame, m, pos):
    plt.subplot(2, 2, 1)
    plt.imshow(frame)
    plt.scatter(x=[pos[0]], y=[pos[1]], c='r', s=10)
    plt.scatter(x=[pos[2]], y=[pos[3]], c='b', s=10)
    
    plt.subplot(2,2,2)
    plt.imshow(frame.astype('double')[:,:,0] - m[:,:,0])
    plt.scatter(x=[pos[0]], y=[pos[1]], c='r', s=10)
    plt.scatter(x=[pos[2]], y=[pos[3]], c='b', s=10)
    
    plt.subplot(2,2,3)
    plt.imshow(frame.astype('double')[:,:,1] - m[:,:,1])
    plt.scatter(x=[pos[0]], y=[pos[1]], c='r', s=10)
    plt.scatter(x=[pos[2]], y=[pos[3]], c='b', s=10)
    
    plt.subplot(2,2,4)
    plt.imshow(frame.astype('double')[:,:,2] - m[:,:,2])
    plt.scatter(x=[pos[0]], y=[pos[1]], c='r', s=10)
    plt.scatter(x=[pos[2]], y=[pos[3]], c='b', s=10)
    
    plt.show()


def ProcessFrame(frame, m):
    try:
        dat = frame.astype('double')
        dat = dat - m
        
        
        r = dat[:,:,0]
        g = dat[:,:,1]
        
        
        rb = r > 0.8*r.max()
        gb = g > 0.8*g.max()
    
            
        lbl, nBlob = scipy.ndimage.label(rb, structure = generate_binary_structure(2,2))
        index=np.unique(lbl)
        mvar = np.zeros([index.max(),])
        for i in range(1, index.max()):
            mvar[i] = r[np.where(lbl==i)].max()
        correctCluster = mvar.argmax()
        rpos = scipy.ndimage.measurements.center_of_mass(rb, lbl, correctCluster)
        
        lbl, nBlob = scipy.ndimage.label(gb, structure = generate_binary_structure(2,2))
        index=np.unique(lbl)
        mvar = np.zeros([index.max(),])
        for i in range(1, index.max()):
            mvar[i] = g[np.where(lbl==i)].max()
        correctCluster = mvar.argmax()
        gpos = scipy.ndimage.measurements.center_of_mass(gb, lbl, correctCluster)
        
        pos = np.concatenate([rpos, gpos])
    except:
        pos = np.empty([1,4])
        pos[:] = np.nan
    return pos

def ReadIt(fileName):
    v = pims.Video(fileName)
    m = GetMean(v)
    
    nFrames = v.get_metadata()['nframes']
    pos = np.zeros([nFrames,4])
 
    for frame in range(nFrames):
        pos[frame,:] = ProcessFrame(v[frame], m)
        #ShowFrame(v[frame], m, pos[frame,:])
        if np.mod(frame,np.floor(nFrames/100))==0:
            print(100*frame/nFrames)

    np.savetxt(fileName[0:-3] + 'csv',pos)


def ReadItPar(fileName):
    print("Reading:", fileName)
    t1 = time.time()
    v = pims.Video(fileName)
    m = GetMean(v)
    num_cores = multiprocessing.cpu_count()
    

    results = Parallel(n_jobs=num_cores)(delayed(ProcessFrame)(frame, m) for frame in v[:])

    pos = np.asarray(results)
    t2 = time.time()
    print('Elapsed: %s' % (t2-t1))
    np.savetxt(fileName[0:-3] + 'csv',pos)


if __name__ == '__main__':
    if int(sys.argv[2]) == 1:
        ReadItPar(sys.argv[1])
    else:
        ReadIt(sys.argv[1])

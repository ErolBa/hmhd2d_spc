import h5py
import numpy as np
import matplotlib.pyplot as plt


#================= Routines to Read Arrays ===========================
def read2(dataset,num,xmin=-np.inf,xmax=np.inf, \
          zmin=-np.inf,zmax=np.inf):
# read in 2D data in the region defined by 
# [xmin,xmax]x[zmin,zmax]

  filename='data'+'%04d'%num+'.hdf'
  f=h5py.File(filename,'r')
  ttime=float(np.array(f['ttime']))
  x=np.array(f['x'])
  z=np.array(f['z']) 

  xid=np.where((x-xmin)*(x-xmax)<=0)[0]
  xs,xe=xid[0],xid[-1]+1

  zid=np.where((z-zmin)*(z-zmax)<=0)[0]
  zs,ze=zid[0],zid[-1]+1

  v=np.array(f[dataset][zs:ze,xs:xe])
  dtype=v.dtype
  if (dtype in ('uint8','uint16','int8','int16')):
    # binned data
    scale=np.array(f[dataset+'_scale'])
    offset=np.array(f[dataset+'_offset'])
    v=v*scale+offset
  f.close() 
  return ttime,x[xs:xe],z[zs:ze],v

#--------------------------------------------------------------------
def read2x(dataset,num,k,xmin=-np.inf,xmax=np.inf):
# read 1D profile along x

  filename='data'+'%04d'%num+'.hdf'
  f=h5py.File(filename,'r')
  ttime=float(np.array(f['ttime']))
  x=np.array(f['x'])
  z=np.array(f['z']) 
  xid=np.where((x-xmin)*(x-xmax)<=0)[0]
  xs,xe=xid[0],xid[-1]+1
  v=np.array(f[dataset][k,xs:xe])
  dtype=v.dtype

  if (dtype in ('uint8','uint16','int8','int16')):
    # binned data
    scale=np.array(f[dataset+'_scale'])
    offset=np.array(f[dataset+'_offset'])
    v=v*scale+offset
  return ttime,x[xs:xe],z[k],v

#--------------------------------------------------------------------
def read2z(dataset,num,i,zmin=-np.inf,zmax=np.inf):
# read 1D profile along z

  filename='data'+'%04d'%num+'.hdf'
  f=h5py.File(filename,'r')
  ttime=float(np.array(f['ttime']))
  x=np.array(f['x'])
  z=np.array(f['z']) 
  zid=np.where((z-zmin)*(z-zmax)<=0)[0]
  zs,ze=zid[0],zid[-1]+1
  v=np.array(f[dataset][zs:ze,i])
  dtype=v.dtype

  if (dtype in ('uint8','uint16','int8','int16')):
    # binned data
    scale=np.array(f[dataset+'_scale'])
    offset=np.array(f[dataset+'_offset'])
    v=v*scale+offset
  return ttime,x[i],z[zs:ze],v


#--------------------------------------------------------------------


def rdhst(filename):
# read history file
  filename=filename+'.hdf'
  f=h5py.File(filename,'r')
  hist=np.array(f['thist'])
  f.close()
  return hist
#--------------------------------------------------------------------

def combinehst(filehead,startid,endid,output):
  """
  Combine many history files together
  filehead :  the head of filename
  startid : start id number
  endid : end id number
  output : output file name
  Note : the file id is 2-digit; file name extension .hdf is assumed
  Newer data will overwrite older data if the time windows overlap.
  """
  hist=np.array([])
  for j in range(startid,endid+1):
    filename=filehead+'%02d'%j
    hist1=rdhst(filename)
    if hist.shape==(0,):
      hist=hist1
    else:  
      hist=hist[:,hist[0,:]<hist1[0,0]]  #find the part that doesn't overlap
      hist=np.hstack((hist,hist1))
  # write file
  filename=output+'.hdf'
  f=h5py.File(filename,'w')
  ds=f.create_dataset('thist',hist.shape,hist.dtype,chunks=True) 
  ds[...]=hist.copy()
  f.close()
  return hist    

#=============== Routines for Visualization ===========================
def plot2x(dataset,num,k,xmin=-np.inf,xmax=np.inf):
# plot 1D profile along x

  ttime,x,z,v=read2x(dataset,num,k,xmin,xmax)
  plt.plot(x,v)
  plt.xlabel('x')
  plt.title(dataset+' at t=%.3f along z=%.3f'%(ttime,z))

#--------------------------------------------------------------------
def plot2z(dataset,num,i,zmin=-np.inf,zmax=np.inf):
# plot 1D profile along z

  ttime,x,z,v=read2z(dataset,num,i,zmin,zmax)
  plt.plot(z,v)
  plt.xlabel('z')
  plt.title(dataset+' at t=%.3f along x=%.3f'%(ttime,x))
#--------------------------------------------------------------------
def plot2(dataset,num,xmin=-np.inf,xmax=np.inf, \
          zmin=-np.inf,zmax=np.inf):
# plot 2D profile along x-z plane

  ttime,x,z,v=read2(dataset,num,xmin,xmax,zmin,zmax)
  plt.clf()
  plt.pcolormesh(x,z,v)
  plt.colorbar()
  plt.xlabel('x')
  plt.ylabel('z')
  plt.title(dataset+' at t=%.3f'%ttime)

#--------------------------------------------------------------------
def movie2x(dataset,start,end,k,interval=1,tpause=0.05, \
            xmin=-np.inf,xmax=np.inf):
# plot a 1D time sequence along x from start to end as 
# a movie

  plt.ion() 
  for num in range(start,end+1,interval):
    plt.clf()
    plot2x(dataset,num,k,xmin,xmax)
    plt.draw()
    plt.pause(tpause)

#--------------------------------------------------------------------
def movie2z(dataset,start,end,i,interval=1,tpause=0.05, \
            zmin=-np.inf,zmax=np.inf):
# plot a 1D time sequence along z from start to end as 
# a movie

  plt.ion() 
  for num in range(start,end+1,interval):
    plt.clf()
    plot2z(dataset,num,i,zmin,zmax)
    plt.draw()
    plt.pause(tpause)

#--------------------------------------------------------------------
def movie2(dataset,start,end,interval=1,tpause=0.05, \
           xmin=-np.inf,xmax=np.inf,zmin=-np.inf,zmax=np.inf):
# plot a 2D time sequence as a movie

  plt.ion() 
  for num in range(start,end+1,interval):
    plot2(dataset,num,xmin,xmax,zmin,zmax)
    plt.draw()
    plt.pause(tpause)
#--------------------------------------------------------------------






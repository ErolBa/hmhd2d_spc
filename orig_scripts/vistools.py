from __future__ import print_function 
import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob as glob

def pcolormeshc(x,y,f,**kwargs):
    x1=np.hstack((1.5*x[0]-0.5*x[1],
                 (x[1:]+x[:-1])*0.5,
                 1.5*x[-1]-0.5*x[-2]))
    y1=np.hstack((1.5*y[0]-0.5*y[1],
                 (y[1:]+y[:-1])*0.5,
                  1.5*y[-1]-0.5*y[-2]))
    return plt.pcolormesh(x1,y1,f,**kwargs)

class dataset:
    def __init__(self,filename,xmin=-np.inf,xmax=np.inf, \
                 zmin=-np.inf,zmax=np.inf):
        self.name=filename
        filelist=glob.glob(filename+'*.hdf')
        if filelist==[]:
            print('Data set does not exist!')
            return
        else:
            f=h5py.File(filelist[0],'r')
            x=np.array(f['x'])
            z=np.array(f['z']) 
            self.xmin=max(xmin,x.min())
            self.xmax=min(xmax,x.max())
            self.zmin=max(zmin,z.min())
            self.zmax=min(zmax,z.max())
            xid=np.where((x>=xmin) & (x<=xmax))[0]
            zid=np.where((z>=zmin) & (z<=zmax))[0]
            self.xs,self.xe=xid[0],xid[-1]+1
            self.zs,self.ze=zid[0],zid[-1]+1
            self.x=x[self.xs:self.xe]
            self.z=z[self.zs:self.ze]
            f.close()
    #--------------------------------------------------------------------
    def time(self,num):
        filename=self.name+'%04d'%num+'.hdf'
        f=h5py.File(filename,'r')
        t=float(np.array(f['ttime']))
        f.close()
        return t 
            
    #================= Routines to Read Arrays ===========================
    def __call__(self,var,num):        
        # read in 2D data in the region defined by 
        # [xmin,xmax]x[zmin,zmax]

        filename=self.name+'%04d'%num+'.hdf'
        f=h5py.File(filename,'r')

        v=np.array(f[var][self.zs:self.ze,self.xs:self.xe])
        dtype=v.dtype
        if (dtype in ('uint8','uint16','int8','int16')):
            # binned data
            scale=np.array(f[var+'_scale'])
            offset=np.array(f[var+'_offset'])
            v=v*scale+offset
        f.close() 
        return v
   #--------------------------------------------------------------------
    def read2x(self,var,num,k):
    # read 1D profile along x
 
        filename=self.name+'%04d'%num+'.hdf'
        f=h5py.File(filename,'r')
        v=np.array(f[var][self.zs+k,self.xs:self.xe])
        dtype=v.dtype

        if (dtype in ('uint8','uint16','int8','int16')):
            # binned data
            scale=np.array(f[var+'_scale'])
            offset=np.array(f[var+'_offset'])
            v=v*scale+offset
        return v
    #--------------------------------------------------------------------
    def read2z(self,var,num,i):
    # read 1D profile along z

        filename=self.name+'%04d'%num+'.hdf'
        f=h5py.File(filename,'r')
        v=np.array(f[var][self.zs:self.ze,self.xs+i])
        dtype=v.dtype

        if (dtype in ('uint8','uint16','int8','int16')):
            # binned data
            scale=np.array(f[var+'_scale'])
            offset=np.array(f[var+'_offset'])
            v=v*scale+offset
        return v

    #=============== Routines for Visualization ===========================
    def plot2x(self,var,num,k,s='',output=False,dpi=100,**kwargs):
    # plot 1D profile along x
        xmin=self.xmin
        xmax=self.xmax
 
        v=self.read2x(var,num,k)
        plt.plot(self.x,v,s,**kwargs)
        plt.xlim([xmin,xmax])
        plt.xlabel('x')
        plt.title(var+' at t=%.3f along z=%.3f'%(self.time(num),self.z[k]))

        if output==True:
            zd=np.floor(np.log10(self.z.size))+1
            append='_z%0'+'%d'%zd+'d_%04d.png'
            fname=var+append%(k,num)
            plt.savefig(fname,dpi=dpi)
        

    #--------------------------------------------------------------------
    def plot2z(self,var,num,i,s='',output=False,dpi=100,**kwargs):
    # plot 1D profile along z
        zmin=self.zmin
        zmax=self.zmax
 
        v=self.read2z(var,num,i)
        plt.plot(self.z,v,s,**kwargs)
        plt.xlim([zmin,zmax])
        plt.xlabel('z')
        plt.title(var+' at t=%.3f along x=%.3f'%(self.time(num),self.x[i]))

        if output==True:
            xd=np.floor(np.log10(self.x.size))+1
            append='_x%0'+'%d'%xd+'d_%04d.png'
            fname=var+append%(i,num)
            plt.savefig(fname,dpi=dpi)
        
    #--------------------------------------------------------------------
    def plot2(self,var,num,style='pcolor',nlevels=15,\
              output=False,dpi=100,**kwargs):
    # plot 2D profile along x-z plane
        xmin=self.xmin
        xmax=self.xmax
        zmin=self.zmin
        zmax=self.zmax
        v=self(var,num)
        plt.clf()
        if style=='pcolor':
            pcolormeshc(self.x,self.z,v,**kwargs)
        elif style=='contour':
            plt.contour(self.x,self.z,v,nlevels,**kwargs)
        elif style=='contourf':
            plt.contourf(self.x,self.z,v,nlevels,**kwargs)
        plt.xlim([xmin,xmax])
        plt.ylim([zmin,zmax])
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('z')
        plt.title(var+' at t=%.3f'%self.time(num))

        if output==True:
            append='_%04d.png'
            fname=var+append%(num)
            plt.savefig(fname,dpi=dpi)

    #--------------------------------------------------------------------
    def movie2x(self,var,start,end,k,s='',interval=1,tpause=0.05,**kwargs):
    # plot a 1D time sequence along x from start to end as 
    # a movie
        plt.ion() 
        for num in range(start,end+1,interval):
            plt.clf()
            self.plot2x(var,num,k,s,**kwargs)
            plt.draw()
            plt.pause(tpause)

    #--------------------------------------------------------------------
    def movie2z(self,var,start,end,i,s='',interval=1,tpause=0.05,**kwargs):
    # plot a 1D time sequence along z from start to end as 
    # a movie
        plt.ion() 
        for num in range(start,end+1,interval):
            plt.clf()
            self.plot2z(var,num,i,s,**kwargs)
            plt.draw()
            plt.pause(tpause)

    #--------------------------------------------------------------------
    def movie2(self,var,start,end,interval=1,tpause=0.05,**kwargs):
    # plot a 2D time sequence as a movie

        plt.ion() 
        for num in range(start,end+1,interval):
            plt.clf()
            self.plot2(var,num,**kwargs)
            plt.draw()
            plt.pause(tpause)


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
            #find the part that doesn't overlap
            hist=hist[:,hist[0,:]<hist1[0,0]] 
            hist=np.hstack((hist,hist1))
    # write file
    filename=output+'.hdf'
    f=h5py.File(filename,'w')
    ds=f.create_dataset('thist',hist.shape,hist.dtype,chunks=True) 
    ds[...]=hist.copy()
    f.close()
    return hist    

#--------------------------------------------------------------------

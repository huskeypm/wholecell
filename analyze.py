#
# For processing hdf files generated by fenics 
# 
from dolfin import *
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata

def InterpolateData(mesh,u,dims=2,mode="line",doplot=False, res=100):    
    if dims==3: 
      raise RuntimeError("Only 2D support right now") 

#    mesh = results.mesh
    #dims = np.max(mesh.coordinates(),axis=0) - np.min(mesh.coordinates(),axis=0)
    mmin = np.min(mesh.coordinates(),axis=0)
    mmax = np.max(mesh.coordinates(),axis=0)
    
    #u = results.u_n.split()[0]
    #u = results.u_n
    up = project(u,FunctionSpace(mesh,"CG",1))
    #(gx,gy,gz) = np.mgrid[0:dims[0]:(res*1j),
    #                      dims[1]/2.:dims[1]/2.:1j,
    #                      0:dims[2]:(res*1j)]
    #if dims==3:

    
    if mode=="slice":
      #(gx,gy,gz) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
      #                      mmin[1]:mmax[1]:(res*1j),  
      #                            0:0:1j]
      #img0 = griddata(mesh.coordinates(),up.vector(),(gx,gy,gz))
      (gx,gy) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
                            mmin[1]:mmax[1]:(res*1j)]  
      img0 = griddata(mesh.coordinates(),up.vector(),(gx,gy))

      #print np.shape(img0)
      img0 = np.reshape(img0[:,:,0],[res,res])

      if doplot:
        plt.pcolormesh(img0)
        plt.colorbar()
        plt.gcf().savefig("test.png")
      return img0

    if mode=="line":
    
      yMid = 0.5*(mmax[1]+mmin[1])
      #(gx,gy,gz) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
      #                      yMid:yMid:1j,  
      #                      0:0:1j]
      #line = griddata(mesh.coordinates(),up.vector(),(gx,gy,gz))

      #(gx,gy) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
      #                      yMid:yMid:1j]  
      tol = 0.1
      yMid = 0.1   
      print mmin,mmax
      (gx,gy) = np.mgrid[mmin[0]:mmax[0]:(res*1j),
                            (yMid):(yMid+tol):100j]  
      line = griddata(mesh.coordinates(),up.vector(),(gx,gy))
      
      print np.shape(line)
      line = np.mean(line,axis=0)
      print np.shape(line)
      #img0 = np.reshape(line[:,0,0],[res])
      line = np.ndarray.flatten(line)
      #print np.shape(line)
      if doplot:
        plt.plot(line)
        plt.gcf().savefig("test1.png") 
      return line 


def ReadHdfSlice(hdfFile,i=0,verbose=True):
  hdf = HDF5File(mpi_comm_world(),hdfFile,'r')
  mesh = Mesh()
  hdf.read(mesh,"mesh",False)
  V = FunctionSpace(mesh,"CG",1)
  u = Function(V)

  dataset = "uCa/vector_%d"%i
  attr = hdf.attributes(dataset)
  if verbose:
      print 'Retrieving time step:', attr['timestamp']
  hdf.read(u, dataset)
  return mesh,u 

def ReadHdf(hdfFile,ssl=False,verbose=False):
  hdf = HDF5File(mpi_comm_world(),hdfFile,'r')
  mesh = Mesh()
  hdf.read(mesh,"mesh",False)
  #hdf.read(z,"nEntries") 
  V = FunctionSpace(mesh,"CG",1)
  u = Function(V)
  R = FunctionSpace(mesh,"R",0)
  ur = Function(R)

  # temp hack for grabbing doubles
  x = Function(V) 
  hdf.read(x,"volSSL")
  volSSL = x.vector()[0]
  hdf.read(x,"volCleft")
  volCleft= x.vector()[0]
  hdf.read(x,"volCyto")  
  volCyto= x.vector()[0]
  print volCleft,volSSL,volCyto
  
  attr = hdf.attributes("uCa")
  nsteps = attr['count']

  # arrays
  concsCa = np.zeros(nsteps) 
  concsCaSSL = np.zeros(nsteps) 
  concsCaCleft = np.zeros(nsteps) 

  
  ts = np.arange(nsteps)
  for i in range(nsteps):
    # uCa Field
    dataset = "uCa/vector_%d"%i
    attr = hdf.attributes(dataset)
    if verbose:
      print 'Retrieving time step:', attr['timestamp']
    hdf.read(u, dataset)
    #print "Assemble %d/%f" % (i,assemble(u*dx))
    concsCa[i] = assemble(u*dx) # /volCyto

    # uCaSSL scalar
    if ssl:
      dataset = "uCaSSL/vector_%d"%i
      attr = hdf.attributes(dataset)
      hdf.read(ur, dataset)
      #print "Assemble %d/%f" % (i,assemble(ur*dx))
      concsCaSSL[i]= assemble(ur*dx) # /volSSL

    # uCaCleft scalar
    dataset = "uCaCleft/vector_%d"%i
    attr = hdf.attributes(dataset)
    hdf.read(ur, dataset)
    #print "Assemble %d/%f" % (i,assemble(ur*dx))
    concsCaCleft[i]= assemble(ur*dx) # /volCleft

    # 
    if verbose:
      print concsCaCleft[i], concsCaSSL[i], concsCa[i]

  hdf.close()
  return ts,concsCaCleft, concsCaSSL, concsCa


#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################

#
# ROUTINE  
#


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -read <hdfname>" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

#
# MAIN routine executed when launching this script from command line 
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  #print "arg"

  # Loops over each argument in the command line 
  ssl = False
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-ssl"):
      ssl = True
    if(arg=="-read"):
      arg1=sys.argv[i+1] 
      ReadHdf(arg1,ssl=ssl,verbose=True)
      quit()
  





  raise RuntimeError("Arguments not understood")





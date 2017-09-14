##
## Case defintion for simulations
##
import numpy as np

class Case:
  def __init__(self,
               tag="hashIdentifier",
               label="plotLabel",
               T = 1000, 
               nIntervals=1,
               args=[], # param defn like "[-T 1000, -Vmax 85]"
               odeName = "shannon_2004.ode",
               name="shannon_2004.pickle"):
    self.tag = tag
    self.label = label
    self.args = args
    self.name = name 
    self.odeName = odeName 
    self.prefix=name # store as prefix, since 'name' gets rewritten in daisychain

    # Grab 'T' from command line args. 
    for arg in args:
      if "-T" in arg:
        T = arg.split()[1]
        T = np.float(T)
        print "Using T=%f from command line args"%T

    self.CommandLine(T,nIntervals)


  def CommandLine(self,T,nIntervals):
   
    #
    Ti = T/nIntervals
    if Ti>1.:
      execName = "python daisychain.py " 
      self.args.append("-iters %d"%nIntervals) 
      pickleName = self.name.replace(".pickle","_%d.pickle"%nIntervals)
      
    else: 
      execName = "python runShannonTest.py -jit"
      Ti = T
      pickleName = self.name 

  
        
    self.args.append("-odeName %s"%self.odeName)   
    self.args.append("-T %f"%Ti)
    self.args.append("-name %s"%self.name)   

    cmdLine = execName+" " 
    cmdLine+= " ".join(self.args)+" " 
    cmdLine+= " &"
    self.cmdLine = cmdLine 

    # update name so we read last file 
    self.name = pickleName 


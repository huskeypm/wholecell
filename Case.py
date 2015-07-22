##
## Case defintion for simulations
##

class Case:
  def __init__(self,
               tag="hashIdentifier",
               label="plotLabel",
               T = 1000, 
               nIntervals=1,
               args=[], # param defn like "[-T 1000, -Vmax 85]"
               name="test.pickle"):
    self.tag = tag
    self.label = label
    self.args = args
    self.name = name 

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

  
        
    self.args.append("-T %f"%Ti)
    self.args.append("-name %s"%self.name)   

    cmdLine = execName+" " 
    cmdLine+= " ".join(self.args)+" " 
    cmdLine+= " &"
    self.cmdLine = cmdLine 

    # update name so we read last file 
    self.name = pickleName 


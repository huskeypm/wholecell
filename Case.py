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
    args.append("-T %f"% T)
    self.args = args
    self.name = name 

    self.CommandLine(T,nIntervals)


  def CommandLine(self,T,nIntervals):
    execName = "python runShannonTest.py -jit"
    cmdLine = execName+" "  +" ".join(self.args)+" -name "+self.name+" &"
    self.cmdLine = cmdLine 


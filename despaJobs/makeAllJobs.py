
import runShannonTest as rst


T = 100* 1e3 # 300s --> ms
varDict = dict()  
varDict["G_CaBk"] = [1.0,2.75,1.25] # define 'fold' range over which G_CaBk is varied: (v[0]..v[1] in increments of v[2])
varDict["G_NaBk"] = [1.0,2.75,1.25]

varDictFixed = dict()
varDictFixed["T"] = 298 # K 
names,allKeys,allVars = rst.GenSweptParams(varDict,varDictFixed,T=T)
len(names)

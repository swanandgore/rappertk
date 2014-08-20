import os
from samplers import PhipsiSampler
cbDatapath = os.environ["RTKROOT"] + "/data/"

class PhipsiSamplerProvider :
    def __init__(self) : self.samplers = {}
    def get(self, resn) :
        if not resn in self.samplers.keys() :
            self.samplers[resn] = PhipsiSampler( cbDatapath + "/PhipsiWeightedProp/ps%s" % resn )
        return self.samplers[resn]

psp = PhipsiSamplerProvider() ## make peptide builders
glySampler = psp.get("GLY")
for i in range(100000) :
    glySampler.printSample()

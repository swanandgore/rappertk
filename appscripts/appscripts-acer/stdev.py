import math

def findMeanStdev(vals) :
    mean, stdev = 0., 0.
    for v in vals : mean += v
    mean /= len(vals)
    for v in vals : stdev += (v-mean) * (v-mean)
    stdev = math.sqrt( stdev / len(vals) )
    return mean, stdev



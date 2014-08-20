import string, sys
filenames = ["scASN", "scCYS", "scGLU", "scHIS", "scLEU", "scMET", "scPRO", "scSER", "scTRP", "scVAL", "scARG", "scASP", "scGLN", "scILE", "scLYS", "scPHE", "scTHR", "scTYR", ]


for fn in filenames :
    freq = {}
    for l in open(fn,'r').readlines() :
        flds = l.split()
        k = flds[0] + flds[1]
        if not k in freq.keys() : freq[k] = 0
        freq[k] += 1
    vals = freq.values()
    vals.sort()
    print fn, vals[0], vals[len(vals)-1]

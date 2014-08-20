'''find the freq of pairs input on std input, one pair per line'''

import sys, string

freq = {}
for l in sys.stdin.readlines() :
    flds = l.split()
    pair = ( string.atoi(flds[0]),string.atoi(flds[1]) )
    try :
        freq[pair] = freq[pair] + 1
    except KeyError :
        freq[pair] = 1

for pair,f in freq.items() :
    print pair[0], pair[1], f

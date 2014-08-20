########## some pdb format info for atom and hetatm lines
atmi = (6,11)
atmn = (12,16)
resn = (17,20)
chid = (21,22)
resi = (22,26)
ic = (26,27)
assert(atmn[0] < resn[0] and atmn[1] <= resn[0])
assert(resn[0] < chid[0] and resn[1] <= chid[0])
assert(chid[0] < resi[0] and chid[1] <= resi[0])
xcrd = (30,38)
ycrd = (38,46)
zcrd = (46,54)
occu = (54,60)
bfac = (60,66)
segi = (72,76)

########## CRYST1 line
a     = ( 6,15)   #a (Angstroms).
b     = (15,24)   #b (Angstroms).
c     = (24,33)   #c (Angstroms).
alpha = (33,40)   #alpha (degrees).
beta  = (40,47)   #beta (degrees).
gamma = (47,54)   #gamma (degrees).
sg    = (55,66)   #Space group.
z     = (66,70)   #Z value.

import math, random

class Anneal :
    def __init__(s, tweakers, pts, startTemp, endTemp, numSteps) :
        s.tweakers = tweakers
        s.pts = pts
        s.Tstart, s.Tend, s.tempSteps = startTemp, endTemp, numSteps

    def backupPoints(s, pts, indices) :
        backup = []
        for i in range(len(indices)) :
            index = indices[i]
            bkup = [ pts[index][0], pts[index][1], pts[index][2] ]
            backup.append(bkup)
        return backup

    def restorePoints(s, pts, backup, indices) :
        for i in range(len(indices)) :
            index = indices[i]
            for k in range(3) : pts[index][k] = backup[i][k]

    def anneal(s) :
        energies = [0.] * len(s.tweakers)
        for tw in s.tweakers : energies.append( tw.energy() )
        dT = (endTemp-startTemp) / (numSteps+0.)
        for ti in range(len(s.tempSteps)) :
            currTemp -= dT * ti
            ptsBackup = None
            while 1 : ## find a tweaker which will have nonzero dE
                twi = int (floor( random() * len(s.tweakers) ))
                ptsBackup = backupPoints(pts, s.tweakers[twi].getOP())
                s.tweakers[twi].build()
                en = s.tweakers[twi].energy()
                if en != energies[twi] :
                    dE = en - energies[twi]
                    energies[twi] = en
                    break
                restorePoints(pts, ptsBackup, s.tweakers[twi].getOP())
            r1 = pow(math.e, -1.*dE/currTemp)
            r2 = random.random()
            if r2 >= r1 : restorePoints(pts, ptsBackup, s.tweakers[twi].getOP())

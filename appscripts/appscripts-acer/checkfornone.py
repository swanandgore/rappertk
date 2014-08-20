
def fixdirnames(dir1):
    lendir = len(dir1)
    #if dir1[lendir-1] == '/' and lendir == 1 :
    #   return dir1

    for dx in range(lendir-1,1,-1):
        if dir1[dx] != '/' :
            dir = dir1[0:dx+1]
            return dir
    return dir1

def checkfornone(pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotLib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel, poorOnly,poorThreshold,loopres,start,startcode,stop,stopcode,chainid,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,closure,addsc,userot):

    if  pdbfile == "None" : pdbfile = None
    if  pdbout == "None" : pdbout= None
    if dir_xyzout == "None" :     dir_xyzout= None
    if mapfn == "None" :     mapfn= None
    if mtzfn == "None" :     mtzfn= None
    if     caRes == "None" :     caRes= None
    if scRes == "None" :    scRes= None
    if caRad == "None" :    caRad= None
    if scRad == "None" :    scRad= None
    if     scReduction == "None" :    scReduction= None
    if     popsize == "None" :     popsize= None
    if     verbose == "None" :     verbose= None
    if backtrack == "None" :     backtrack= None
    if     rotLib == "None" :     rotLib= None
    if nmodels == "None" :     nmodels= None
    if     mconly == "None" :     mconly= None
    if     sconly == "None" :    sconly= None
    if     opsax == "None" :     opsax= None
    if     natt == "None" :     natt= None
    if     cacaCutoff == "None" :     cacaCutoff= None
    if     a == "None" :     a= None
    if     b == "None" :     b= None
    if     c == "None" :     c= None
    if     alpha == "None" :     alpha= None
    if     beta == "None" :     beta= None
    if     gamma == "None" :     gamma= None
    if     sg == "None" :      sg= None
    if     resolution == "None" :     resolution= None
    if     f1label == "None" :    f1label= None
    if     sigf1label == "None" :    sigf1label= None
    if     f2label == "None" :     f2label= None
    if     phiclabel == "None" :     phiclabel= None
    if     usefreer == "None" :     usefreer= None
    if     freeRlabel == "None" :     freeRlabel= None

    if     poorOnly == "None" :     poorOnly= None
    if     poorThreshold == "None" :     poorThreshold= None
    if     loopres == "None" :     loopres= None
    if     start == "None" :     start= None
    if     startcode == "None" :     startcode = None
    if     stop == "None" :    stop= None
    if     stopcode == "None" :    stopcode = None
    if     chainid == "None" :     chainid= None

    if     nativeBfac == "None" :     nativeBfac= None
    if     mcBfac == "None" :     mcBfac= None
    if     scBfac == "None" :     scBfac= None
    if     minXSig == "None" :     minXSig= None
    if     maxXSig == "None" :     maxXSig= None
    if     edOpt == "None" :    edOpt= None
    if     allOpt == "None" :     allOpt= None
    if     closure == "None" :     closure = None
    if     addsc == "None" :     addsc = None
    if     userot == "None" :     userot = None



    if mconly not in ["True", "False"]:
        print "--mconly needs to be either True or False. Incorrect value input " , mconly
        import sys ; sys.exit()

    if userot not in ["True", "False"]:
        print "--userot needs to be either True or False. Incorrect value input " , userot
        import sys ; sys.exit()

    if addsc not in ["True", "False"]:
        print "--addsc needs to be either True or False. Incorrect value input " , addsc
        import sys ; sys.exit()



        


    if sconly not in ["True", "False"]:
        print "--sconly needs to be either True or False. Incorrect value input " , sconly
        import sys ; sys.exit()

    if mconly == 'True' and sconly == 'True':
        print "Choose between building mainchain modelling (--mconly True --sconly False) or sidechain only modelling (--mconly False --sconly True) or all atom modelling (--mconly False --sconly False) ";
        import sys ; sys.exit()


    if scRes not in ["True", "False"]:
        print "--use-sc-restraints needs to be either True or False, Incorrect value input " ,scRes
        import sys ; sys.exit()

    if caRes not in ["True", "False"]:
        print "--use-ca-restraints needs to be either True or False , Incorrect value input",caRes
        import sys ; sys.exit()

    if nativeBfac not in ["True","False"]:
        print "--models-get-native-bfactors needs to be either True or False. Incorrect value input",nativeBfac
        import sys ; sys.exit()

    if opsax not in ["True","False"]:
        print "--opsax needs to be either True or False.  Incorrect value input",opsax
        import sys ; sys.exit()

    if usefreer not in ["True","False"]:
        print "--use-freer needs to be either True or False. Incorrect value input",usefreer
        import sys ; sys.exit()
        
    if edOpt  not in ["True","False"]:
        print "--make-ed-optional needs to be either True or False.  Incorrect value input",edOpt
        import sys ; sys.exit()

    if allOpt not in ["True","False"] :
        print "--make-all-restraints-optional needs to be either True or False .  Incorrect value input",allOpt
        import sys ; sys.exit()
        
    if poorOnly not in ["True","False"]:
        print "--rebuild-poor-regions-only needs to be either True or False,  Incorrect value input",poorOnly
        import sys ; sys.exit()

    if closure not in ["True","False"]:
        print "--use-loopclosure-restraints  needs to be either True or False,  Incorrect value input",closure
        import sys ; sys.exit()        



        

    if (scReduction < 0. or scReduction > 1.) and mconly != "True" :
        print "value of --sidechain-vdw-reduction needs to be >= 0.0 and <= 1." ;
        print "Incorrect value given", scReduction
        import sys ; sys.exit()
        
    if rotLib not in ["PRL", "SCL1.0","SCL0.5","SCL0.2"]:
        print "Unrecognised rotamer library, options are PRL, SCL1.0, SCL0.5 ,  SCL0.2"
        print "Incorrect value given", rotLib
        import sys ; sys.exit()


        
    if caRad < 0. and caRes == "True" :
        print "value of --ca-restraint-radius  needs to be positive value " ;
        print "Incorrect value given", caRad
        import sys ; sys.exit()

    if scRad < 0. and scRes == "True" :
        print "value of --sc-centroid-restraint-radius  needs to be positive value " ;
        print "Incorrect value given", scRad
        import sys ; sys.exit()


    if popsize < 100  :
        print "WARNING !! value of --population-size  should be atleast 100"
        if popsize < 0  :
            print "value of --population-size  should be atleast 100"
            print "Incorrect value given", popsize
            import sys ; sys.exit()                

    if verbose < 1  :
        print "WARNING !! value of --verbose  should be atleast 1"
        print "Setting to default value of 7"
        verbose = 7
    
    if natt < 1  :
        print "WARNING !! value of --attempts  should be atleast 1"
        print "Setting to default value of 5"
        natt = 5


    if cacaCutoff < 0. :
        print "WARNING !! value of --cacaCutoff  should be positive"
        print "Setting to default value of 5.0"
        cacaCutoff = 5.0

    if backtrack == None :
        backtrack = None
    else :
        if backtrack != None : 
            xpos = backtrack.find('X')
            if xpos <= 0:
                print "Error in specifying backtrack, it should be of the form of numstepsXstepsize, where numsteps and stepsize are integer values. For e.g. --backtrack 4X5 , in case of failure to build at a residue the program will  backtrack 4 times by 5 residues."
                print "Incorrect value given", backtrack
                import sys ; sys.exit()

            try : 
                int(backtrack[0:xpos]) ; int(backtrack[xpos+1:])

            except ValueError : 
                print "Error in specifying backtrack, it should be of the form of numstepsXstepsize, where numsteps and stepsize are integer values. For e.g. --backtrack 4X5 , in case of failure to build at a residue the program will  backtrack 4 times by 5 residues."
                print "Incorrect value given ", backtrack
                import sys ; sys.exit()



    if opsax == "True" and (mapfn == None and mtzfn == None) :
        print "WARNING!! OPSAX will  only be used when MTZ/MAP file is given"
        opsax = "False"





    print "========================== Now  printing all KEYWORD values ==========================="




    print
    print "                           Directories and Files"
    print
    
    print "	--xyzin ".ljust(30), 
    print "%10s "%pdbfile 
    print "	--xyzout ".ljust(30), 
    print "%10s "%pdbout 
    print "	--dir-xyzout ".ljust(30), 
    print "%10s "%dir_xyzout 
    print "	--mapin ".ljust(30), 
    print "%10s "%mapfn 
    print "	--hklin ".ljust(30), 
    print "%10s "%mtzfn 
    print "	--use-ca-restraints ".ljust(20), 
    print "%10s "%caRes
    
    print
    print "                           Restraint paramaters "
    print


    print "	--use-sc-restraints ".ljust(30), 
    print "%10s "%scRes 
    print "	--ca-restraint-radius ".ljust(30), 
    print "%10s "%caRad 
    print "	--sc-centroid-restraint-radius ".ljust(30), 
    print "%10s "%scRad

    


    print "	--use-loopclosure-restraints ".ljust(30), 
    print "%10s "%closure
    

    print
    print "                           Build paramaters "
    print
    print "	--sidechain-vdw-reduction ".ljust(30), 
    print "%10s "%scReduction 
    print "	--population-size ".ljust(30), 
    print "%10s "%popsize 
    print "	--verbose ".ljust(30),
    print "%10s "%verbose
    
    print "	--attempts ".ljust(30), 
    print "%10s "%natt



    print "	--backtrack ".ljust(30), 
    print "%10s "%backtrack 
    print "	--mconly ".ljust(30), 
    print "%10s "%mconly 
    print "	--sconly ".ljust(30), 
    print "%10s "%sconly 
    print "	--make-all-restraints-optional ".ljust(30), 
    print "%10s "%allOpt 


    print "	--loopseq ".ljust(30), 
    print "%10s "%loopres 
    print "	--start ".ljust(30), 
    print "%10s "%start 
    print "	--stop ".ljust(30), 
    print "%10s "%stop 


    print "	--start-inscode ".ljust(30), 
    print "%10s "%startcode
    print "	--stop-inscode ".ljust(30), 
    print "%10s "%stopcode

    print "	--chainid ".ljust(30), 
    print "%10s "%chainid

    
    
    
    
    print
    print "                           Crystallographic paramaters "
    print
    print "	--opsax ".ljust(30), 
    print "%10s "%opsax 
    print "	--use-freer ".ljust(30), 
    print "%10s "%usefreer 
    print "	--rotamerlib ".ljust(30), 
    print "%10s "%rotLib 
    print "	--num-models ".ljust(30), 
    print "%10s "%nmodels 
    print "	--a ".ljust(30), 
    print "%10s "%a 
    print "	--b ".ljust(30), 
    print "%10s "%b 
    print "	--c ".ljust(30), 
    print "%10s "%c 
    print "	--alpha ".ljust(30), 
    print "%10s "%alpha 
    print "	--beta ".ljust(30), 
    print "%10s "%beta 
    print "	--gamma ".ljust(30), 
    print "%10s "%gamma 
    print "	--sg ".ljust(30), 
    print "%10s "%sg 
    print "	--resolution ".ljust(30), 
    print "%10s "%resolution 
    print "	--FP ".ljust(30), 
    print "%10s "%f1label 
    print "	--SIGFP ".ljust(30), 
    print "%10s "%sigf1label 
    print "	--FC ".ljust(30), 
    print "%10s "%f2label 
    print "	--PHI ".ljust(30), 
    print "%10s "%phiclabel 
    print "	--FREER ".ljust(30), 
    print "%10s "%freeRlabel 



    print "	--rebuild-poor-regions-only ".ljust(30), 
    print "%10s "%poorOnly 
    print "	--poor-fit-threshold ".ljust(30), 
    print "%10s "%poorThreshold 
    print "	--default-mainchain-b-factor ".ljust(30), 
    print "%10s "%mcBfac 
    print "	--default-sidechain-b-factor ".ljust(30), 
    print "%10s "%scBfac 
    print "	--models-get-native-bfactors ".ljust(30), 
    print "%10s "%nativeBfac 
    print "	--minimum-sig ".ljust(30), 
    print "%10s "%minXSig 
    print "	--maximum-sig ".ljust(30), 
    print "%10s "%maxXSig 
    print "	--cacaCutoff ".ljust(30), 
    print "%10s "%cacaCutoff 
    print "	--make-ed-optional ".ljust(30), 
    print "%10s "%edOpt 

    
 
    print
    print
    print "========================== End KEYWORD values =========================================="
    print
    print
    

    
    
    return pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotLib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel, poorOnly,poorThreshold,loopres,start,startcode,stop,stopcode,chainid,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,closure,addsc,userot 














def checkfornone_ppi(pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotLib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel,poorOnly,poorThreshold,loopres,start,startcode,stop,stopcode,chainid,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,closure,addsc,userot):

    if  pdbfile == "None" : pdbfile = None
    if  pdbout == "None" : pdbout= None
    if dir_xyzout == "None" :     dir_xyzout= None
    if mapfn == "None" :     mapfn= None
    if mtzfn == "None" :     mtzfn= None
    if     caRes == "None" :     caRes= None
    if scRes == "None" :    scRes= None
    if caRad == "None" :    caRad= None
    if scRad == "None" :    scRad= None
    if     scReduction == "None" :    scReduction= None
    if     popsize == "None" :     popsize= None
    if     verbose == "None" :     verbose= None
    if backtrack == "None" :     backtrack= None
    if     rotLib == "None" :     rotLib= None
    if nmodels == "None" :     nmodels= None
    if     mconly == "None" :     mconly= None
    if     sconly == "None" :    sconly= None
    if     opsax == "None" :     opsax= None
    if     natt == "None" :     natt= None
    if     cacaCutoff == "None" :     cacaCutoff= None
    if     a == "None" :     a= None
    if     b == "None" :     b= None
    if     c == "None" :     c= None
    if     alpha == "None" :     alpha= None
    if     beta == "None" :     beta= None
    if     gamma == "None" :     gamma= None
    if     sg == "None" :      sg= None
    if     resolution == "None" :     resolution= None
    if     f1label == "None" :    f1label= None
    if     sigf1label == "None" :    sigf1label= None
    if     f2label == "None" :     f2label= None
    if     phiclabel == "None" :     phiclabel= None
    if     usefreer == "None" :     usefreer= None
    if     freeRlabel == "None" :     freeRlabel= None


    if     nativeBfac == "None" :     nativeBfac= None
    if     mcBfac == "None" :     mcBfac= None
    if     scBfac == "None" :     scBfac= None
    if     minXSig == "None" :     minXSig= None
    if     maxXSig == "None" :     maxXSig= None
    if     edOpt == "None" :    edOpt= None
    if     allOpt == "None" :     allOpt= None



    if mconly not in ["True", "False"]:
        print "--mconly needs to be either True or False. Incorrect value input " , mconly
        import sys ; sys.exit()


    if sconly not in ["True", "False"]:
        print "--sconly needs to be either True or False. Incorrect value input " , sconly
        import sys ; sys.exit()

    if mconly == 'True' and sconly == 'True':
        print "Choose between building mainchain modelling (--mconly True --sconly False) or sidechain only modelling (--mconly False --sconly True) or all atom modelling (--mconly False --sconly False) ";
        import sys ; sys.exit()


    if scRes not in ["True", "False"]:
        print "--use-sc-restraints needs to be either True or False, Incorrect value input " ,scRes
        import sys ; sys.exit()

    if caRes not in ["True", "False"]:
        print "--use-ca-restraints needs to be either True or False , Incorrect value input",caRes
        import sys ; sys.exit()

    if nativeBfac not in ["True","False"]:
        print "--models-get-native-bfactors needs to be either True or False. Incorrect value input",nativeBfac
        import sys ; sys.exit()

    if opsax not in ["True","False"]:
        print "--opsax needs to be either True or False.  Incorrect value input",opsax
        import sys ; sys.exit()

    if usefreer not in ["True","False"]:
        print "--use-freer needs to be either True or False. Incorrect value input",usefreer
        import sys ; sys.exit()
        
    if edOpt  not in ["True","False"]:
        print "--make-ed-optional needs to be either True or False.  Incorrect value input",edOpt
        import sys ; sys.exit()

    if allOpt not in ["True","False"] :
        print "--make-all-restraints-optional needs to be either True or False .  Incorrect value input",allOpt
        import sys ; sys.exit()
        

        

    if (scReduction < 0. or scReduction > 1.) and mconly != "True" :
        print "value of --sidechain-vdw-reduction needs to be >= 0.0 and <= 1." ;
        print "Incorrect value given", scReduction
        import sys ; sys.exit()
        
    if rotLib not in ["PRL", "SCL1.0","SCL0.5","SCL0.2"]:
        print "Unrecognised rotamer library, options are PRL, SCL1.0, SCL0.5 ,  SCL0.2"
        print "Incorrect value given", rotLib
        import sys ; sys.exit()


        
    if caRad < 0. and caRes == "True" :
        print "value of --ca-restraint-radius  needs to be positive value " ;
        print "Incorrect value given", caRad
        import sys ; sys.exit()

    if scRad < 0. and scRes == "True" :
        print "value of --sc-centroid-restraint-radius  needs to be positive value " ;
        print "Incorrect value given", scRad
        import sys ; sys.exit()


    if popsize < 100  :
        print "WARNING !! value of --population-size  should be atleast 100"
        if popsize < 0  :
            print "value of --population-size  should be atleast 100"
            print "Incorrect value given", popsize
            import sys ; sys.exit()                

    if verbose < 1  :
        print "WARNING !! value of --verbose  should be atleast 1"
        print "Setting to default value of 7"
        verbose = 7
    
    if natt < 1  :
        print "WARNING !! value of --attempts  should be atleast 1"
        print "Setting to default value of 5"
        natt = 5


    if cacaCutoff < 0. :
        print "WARNING !! value of --cacaCutoff  should be positive"
        print "Setting to default value of 5.0"
        cacaCutoff = 5.0

    if backtrack == None :
        backtrack = None
    else :
        if backtrack != None : 
            xpos = backtrack.find('X')
            if xpos <= 0:
                print "Error in specifying backtrack, it should be of the form of numstepsXstepsize, where numsteps and stepsize are integer values. For e.g. --backtrack 4X5 , in case of failure to build at a residue the program will  backtrack 4 times by 5 residues."
                print "Incorrect value given", backtrack
                import sys ; sys.exit()

            try : 
                int(backtrack[0:xpos]) ; int(backtrack[xpos+1:])

            except ValueError : 
                print "Error in specifying backtrack, it should be of the form of numstepsXstepsize, where numsteps and stepsize are integer values. For e.g. --backtrack 4X5 , in case of failure to build at a residue the program will  backtrack 4 times by 5 residues."
                print "Incorrect value given ", backtrack
                import sys ; sys.exit()



    if opsax == "True" and (mapfn == None and mtzfn == None) :
        print "WARNING!! OPSAX will  only be used when MTZ/MAP file is given"
        opsax = "False"





    print "========================== Now  printing all KEYWORD values ==========================="




    print
    print "                           Directories and Files"
    print
    
    print "	--xyzin ".ljust(30), 
    print "%10s "%pdbfile 
    print "	--xyzout ".ljust(30), 
    print "%10s "%pdbout 
    print "	--dir-xyzout ".ljust(30), 
    print "%10s "%dir_xyzout 
    print "	--mapin ".ljust(30), 
    print "%10s "%mapfn 
    print "	--hklin ".ljust(30), 
    print "%10s "%mtzfn 
    print "	--use-ca-restraints ".ljust(20), 
    print "%10s "%caRes
    
    print
    print "                           Restraint paramaters "
    print


    print "	--use-sc-restraints ".ljust(30), 
    print "%10s "%scRes 
    print "	--ca-restraint-radius ".ljust(30), 
    print "%10s "%caRad 
    print "	--sc-centroid-restraint-radius ".ljust(30), 
    print "%10s "%scRad

    


    

    print
    print "                           Build paramaters "
    print
    print "	--sidechain-vdw-reduction ".ljust(30), 
    print "%10s "%scReduction 
    print "	--population-size ".ljust(30), 
    print "%10s "%popsize 
    print "	--verbose ".ljust(30),
    print "%10s "%verbose
    
    print "	--attempts ".ljust(30), 
    print "%10s "%natt



    print "	--backtrack ".ljust(30), 
    print "%10s "%backtrack 
    print "	--mconly ".ljust(30), 
    print "%10s "%mconly 
    print "	--sconly ".ljust(30), 
    print "%10s "%sconly 
    print "	--make-all-restraints-optional ".ljust(30), 
    print "%10s "%allOpt 


    
    print
    print "                           Crystallographic paramaters "
    print
    print "	--opsax ".ljust(30), 
    print "%10s "%opsax 
    print "	--use-freer ".ljust(30), 
    print "%10s "%usefreer 
    print "	--rotamerlib ".ljust(30), 
    print "%10s "%rotLib 

    print "	--num-models ".ljust(30), 
    print "%10s "%nmodels 


    print "	--a ".ljust(30), 
    print "%10s "%a 
    print "	--b ".ljust(30), 
    print "%10s "%b 
    print "	--c ".ljust(30), 
    print "%10s "%c 
    print "	--alpha ".ljust(30), 
    print "%10s "%alpha 
    print "	--beta ".ljust(30), 
    print "%10s "%beta 
    print "	--gamma ".ljust(30), 
    print "%10s "%gamma 
    print "	--sg ".ljust(30), 
    print "%10s "%sg 
    print "	--resolution ".ljust(30), 
    print "%10s "%resolution 
    print "	--FP ".ljust(30), 
    print "%10s "%f1label 
    print "	--SIGFP ".ljust(30), 
    print "%10s "%sigf1label 
    print "	--FC ".ljust(30), 
    print "%10s "%f2label 
    print "	--PHI ".ljust(30), 
    print "%10s "%phiclabel 
    print "	--FREER ".ljust(30), 
    print "%10s "%freeRlabel 




    print "	--default-mainchain-b-factor ".ljust(30), 
    print "%10s "%mcBfac 
    print "	--default-sidechain-b-factor ".ljust(30), 
    print "%10s "%scBfac 
    print "	--models-get-native-bfactors ".ljust(30), 
    print "%10s "%nativeBfac 
    print "	--minimum-sig ".ljust(30), 
    print "%10s "%minXSig 
    print "	--maximum-sig ".ljust(30), 
    print "%10s "%maxXSig 
    print "	--cacaCutoff ".ljust(30), 
    print "%10s "%cacaCutoff 
    print "	--make-ed-optional ".ljust(30), 
    print "%10s "%edOpt 

    
 
    print
    print
    print "========================== End KEYWORD values =========================================="
    print
    print
    

    
    
    return pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotLib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel,poorOnly,poorThreshold,loopres,start,startcode,stop,stopcode,chainid,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,closure,addsc,userot 


if __name__ == "__main__" :
    print 
    

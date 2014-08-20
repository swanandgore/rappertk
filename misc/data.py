#{{{ number of rotameric states for various sidechains in Shetty libraries



three2one = {}
three2one["CYS"] = 'C'    
three2one["HIS"] = 'H'    
three2one["ILE"] = 'I'     	
three2one["MET"] = 'M'     	
three2one["SER"] = 'S'    
three2one["VAL"] = 'V'     	
three2one["ALA"] = 'A'    
three2one["GLY"] = 'G'     	
three2one["LEU"] = 'L'     	
three2one["PRO"] = 'P'     	
three2one["THR"] = 'T'    
three2one["PHE"] = 'F'    
three2one["ARG"] = 'R'    
three2one["TYR"] = 'Y'    
three2one["TRP"] = 'W'     	
three2one["ASP"] = 'D'     	
three2one["ASN"] = 'N'     	
three2one["GLU"] = 'E'     	
three2one["GLN"] = 'Q'     	
three2one["LYS"] = 'K'
three2one["MSE"] = 'M'




numRotSCL = {}
numRotSCL["0.2"] = {}
numRotSCL["0.2"]["ARG"] = 1240
numRotSCL["0.2"]["ASN"] = 248
numRotSCL["0.2"]["ASP"] = 204
numRotSCL["0.2"]["CYS"] = 14
numRotSCL["0.2"]["GLN"] = 667
numRotSCL["0.2"]["GLU"] = 582
numRotSCL["0.2"]["HIS"] = 344
numRotSCL["0.2"]["ILE"] = 79
numRotSCL["0.2"]["LEU"] = 207
numRotSCL["0.2"]["LYS"] = 670
numRotSCL["0.2"]["MET"] = 393
numRotSCL["0.2"]["MSE"] = 393
numRotSCL["0.2"]["PHE"] = 354
numRotSCL["0.2"]["PRO"] = 14
numRotSCL["0.2"]["SER"] = 26
numRotSCL["0.2"]["THR"] = 22
numRotSCL["0.2"]["TRP"] = 459
numRotSCL["0.2"]["TYR"] = 435
numRotSCL["0.2"]["VAL"] = 30
numRotSCL["0.5"] = {}
numRotSCL["0.5"]["ARG"] = 415
numRotSCL["0.5"]["ASN"] = 48
numRotSCL["0.5"]["ASP"] = 35
numRotSCL["0.5"]["CYS"] = 4
numRotSCL["0.5"]["GLN"] = 148
numRotSCL["0.5"]["GLU"] = 108
numRotSCL["0.5"]["HIS"] = 62
numRotSCL["0.5"]["ILE"] = 18
numRotSCL["0.5"]["LEU"] = 36
numRotSCL["0.5"]["LYS"] = 195
numRotSCL["0.5"]["MET"] = 85
numRotSCL["0.5"]["MSE"] = 85
numRotSCL["0.5"]["PHE"] = 55
numRotSCL["0.5"]["PRO"] = 3
numRotSCL["0.5"]["SER"] = 8
numRotSCL["0.5"]["THR"] = 5
numRotSCL["0.5"]["TRP"] = 105
numRotSCL["0.5"]["TYR"] = 88
numRotSCL["0.5"]["VAL"] = 8
numRotSCL["1.0"] = {}
numRotSCL["1.0"]["ARG"] = 135
numRotSCL["1.0"]["ASN"] = 15
numRotSCL["1.0"]["ASP"] = 11
numRotSCL["1.0"]["CYS"] = 3
numRotSCL["1.0"]["GLN"] = 43
numRotSCL["1.0"]["GLU"] = 30
numRotSCL["1.0"]["HIS"] = 21
numRotSCL["1.0"]["ILE"] = 6
numRotSCL["1.0"]["LEU"] = 16
numRotSCL["1.0"]["LYS"] = 46
numRotSCL["1.0"]["MET"] = 33
numRotSCL["1.0"]["MSE"] = 33
numRotSCL["1.0"]["PHE"] = 20
numRotSCL["1.0"]["PRO"] = 1
numRotSCL["1.0"]["SER"] = 4
numRotSCL["1.0"]["THR"] = 3
numRotSCL["1.0"]["TRP"] = 32
numRotSCL["1.0"]["TYR"] = 23
numRotSCL["1.0"]["VAL"] = 4
#}}}

#{{{ number of rotameric states for various sidechains in Dunbrack and Richardson libraries
numRotDBB = {} # in Dunbrack's backbone dependent rotamer library
numRotDBB["GLY"] = 0
numRotDBB["ALA"] = 0
numRotDBB["PRO"] = 2
numRotDBB["CYS"] = 3
numRotDBB["SER"] = 3
numRotDBB["VAL"] = 3
numRotDBB["THR"] = 3
numRotDBB["PHE"] = 6
numRotDBB["TYR"] = 6
numRotDBB["ASP"] = 9
numRotDBB["ILE"] = 9
numRotDBB["HIS"] = 9
numRotDBB["LEU"] = 9
numRotDBB["TRP"] = 9
numRotDBB["ASN"] = 18
numRotDBB["GLU"] = 27
numRotDBB["MET"] = 27
numRotDBB["MSE"] = 27
numRotDBB["GLN"] = 36
numRotDBB["ARG"] = 81
numRotDBB["LYS"] = 81

numRotPRL = {} # in Richardson's Penultimate rotamer lib
numRotPRL["VAL"] = 3
numRotPRL["PRO"] = 2
numRotPRL["PHE"] = 5
numRotPRL["ASN"] = 9
numRotPRL["THR"] = 2
numRotPRL["LYS"] = 20
numRotPRL["TYR"] = 5
numRotPRL["SER"] = 3
numRotPRL["ILE"] = 6
numRotPRL["ASP"] = 8
numRotPRL["GLN"] = 14
numRotPRL["GLU"] = 14
numRotPRL["MET"] = 12
numRotPRL["MSE"] = 12
numRotPRL["LEU"] = 5
numRotPRL["TRP"] = 6
numRotPRL["ARG"] = 29
numRotPRL["CYS"] = 3
numRotPRL["HIS"] = 6
numRotPRL["GLY"] = 0
numRotPRL["ALA"] = 0

#numRotPRL["GLY"] = 0
#numRotPRL["ALA"] = 0
#numRotPRL["ARG"] = 34
#numRotPRL["ASN"] = 7
#numRotPRL["ASP"] = 5
#numRotPRL["CYS"] = 3
#numRotPRL["GLN"] = 9
#numRotPRL["GLU"] = 8
#numRotPRL["HIS"] = 8
#numRotPRL["ILE"] = 7
#numRotPRL["LEU"] = 5
#numRotPRL["LYS"] = 27
#numRotPRL["MET"] = 13
#numRotPRL["MSE"] = 13
#numRotPRL["PHE"] = 4
#numRotPRL["PRO"] = 2
#numRotPRL["SER"] = 3
#numRotPRL["THR"] = 3
#numRotPRL["TRP"] = 7
#numRotPRL["TYR"] = 4
#numRotPRL["VAL"] = 3
#}}}

# {{{ atoms in AA residues and their covalent connectivity
resAtoms = {
    'GLY':[' N  ', ' CA ', ' C  ', ' O  '],
    'GLU':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' OE1', ' OE2'],
    'ALA':[' N  ', ' CA ', ' C  ', ' O  ', ' CB '],
    'VAL':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG1', ' CG2'],
    'ARG':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' NE ', ' CZ ', ' NH1', ' NH2'],

    'MET':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' SD ', ' CE ',],
    'MSE':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', 'SE  ', ' CE ',],

    'ASP':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' OD1', ' OD2',],
    'ASN':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' OD1', ' ND2',],
    'PRO':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD '],
    'LEU':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2'],
    'HIS':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' ND1', ' CD2', ' CE1', ' NE2'],
    'GLN':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' OE1', ' NE2'],
    'SER':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' OG '],
    'LYS':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' CE ', ' NZ '],
    'PHE':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ '],
    'TYR':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ', ' OH '],
    'THR':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' OG1', ' CG2'],
    'ILE':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG1', ' CG2', ' CD1'],
    'CYS':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' SG '],
    'TRP':[' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' NE1', ' CE2', ' CE3', ' CZ2', ' CZ3', ' CH2'],
}
mcConn = [
    (' N  ',' CA '), (' CA ',' C  '), (' C  ',' O  ')
]
scConn = {
    'GLY' : [],
    'ALA' : [],
    'VAL' : [ (' CB ',' CG1'), (' CB ',' CG2') ],
    'LEU' : [ (' CB ',' CG '), (' CG ',' CD1'), (' CG ',' CD2') ],
    'ILE' : [ (' CB ',' CG1'), (' CB ',' CG2'), (' CG1',' CD1') ],
    'PRO' : [ (' CB ',' CG '), (' CG ',' CD '), (' CD ',' N  ') ],
    'ASP' : [ (' CB ',' CG '), (' CG ',' OD1'), (' CG ',' OD2') ],
    'GLU' : [ (' CB ',' CG '), (' CG ',' CD '), (' CD ',' OE1'), (' CD ',' OE2') ],
    'HIS' : [ (' CB ',' CG '), (' CG ',' ND1'), (' CG ',' CD2'), (' ND1',' CE1'), (' CD2',' NE2'), (' CE1',' NE2') ],
    'LYS' : [ (' CB ',' CG '), (' CG ',' CD '), (' CD ',' CE '), (' CE ',' NZ ') ],
    'ARG' : [ (' CB ',' CG '), (' CG ',' CD '), (' CD ',' NE '), (' NE ',' CZ '), (' CZ ',' NH1'), (' CZ ',' NH2') ],
    'CYS' : [ (' CB ',' SG ') ],
    'MET' : [ (' CB ',' CG '), (' CG ',' SD '), (' SD ',' CE ') ],
    'MSE' : [ (' CB ',' CG '), (' CG ','SE  '), ('SE  ',' CE ') ],
    'SER' : [ (' CB ',' OG ') ],
    'THR' : [ (' CB ',' OG1'), (' CB ',' CG2') ],
    'ASN' : [ (' CB ',' CG '), (' CG ',' OD1'), (' CG ',' ND2') ],
    'GLN' : [ (' CB ',' CG '), (' CG ',' CD '), (' CD ',' OE1'), (' CD ',' NE2') ],
    'PHE' : [ (' CB ',' CG '), (' CG ',' CD1'), (' CG ',' CD2'), (' CD1',' CE1'), (' CD2',' CE2'), (' CE1',' CZ '), (' CE2',' CZ ') ],
    'TYR' : [ (' CB ',' CG '), (' CG ',' CD1'), (' CG ',' CD2'), (' CD1',' CE1'), (' CD2',' CE2'), (' CE1',' CZ '), (' CE2',' CZ '), (' CZ ',' OH ') ],
    'TRP' : [ (' CB ',' CG '), (' CG ',' CD1'), (' CG ',' CD2'), (' CD1',' NE1'), (' NE1',' CE2'), (' CD2',' CE2'), (' CE2',' CZ2'), (' CZ2',' CH2'), (' CH2',' CZ3'), (' CZ3',' CE3'), (' CE3',' CD2') ],
}

chiRes = {
    'GLY' : [],
    'ALA' : [],
    'VAL' : [ (' N  ', ' CA ', ' CB ',' CG1') ],
    'LEU' : [ (' N  ', ' CA ', ' CB ',' CG '), (' CA ', ' CB ',' CG ',' CD1'), ],
    'ILE' : [ (' N  ', ' CA ', ' CB ',' CG1'), (' CA ', ' CB ',' CG1',' CD1'), ],
    'PRO' : [ (' N  ', ' CA ', ' CB ',' CG '), (' CA ', ' CB ',' CG ',' CD '), ],
    'ASP' : [ (' N  ', ' CA ', ' CB ',' CG '), (' CA ', ' CB ',' CG ',' OD1'), ],
    'GLU' : [ (' N  ', ' CA ', ' CB ',' CG '), (' CA ', ' CB ',' CG ',' CD '), (' CB ', ' CG ',' CD ',' OE1'), ],
    'HIS' : [ (' N  ', ' CA ', ' CB ',' CG '), (' CA ', ' CB ',' CG ',' ND1'), ],
    'LYS' : [ (' N  ', ' CA ', ' CB ',' CG '), (' CA ', ' CB ',' CG ',' CD '), (' CB ', ' CG ',' CD ',' CE '), (' CG ', ' CD ',' CE ',' NZ '), ],
    'ARG' : [ (' N  ', ' CA ', ' CB ',' CG '), (' CA ', ' CB ',' CG ',' CD '), (' CB ', ' CG ',' CD ',' NE '), (' CG ', ' CD ',' NE ',' CZ '), ],
    'CYS' : [ (' N  ', ' CA ', ' CB ',' SG ') ],
    'MET' : [ (' N  ', ' CA ', ' CB ',' CG '), (' CA ', ' CB ',' CG ',' SD '), (' CB ', ' CG ',' SD ',' CE '), ],
    'MSE' : [ (' N  ', ' CA ', ' CB ',' CG '), (' CA ', ' CB ',' CG ','SE  '), (' CB ', ' CG ','SE  ',' CE '), ],
    'SER' : [ (' N  ',' CA ',' CB ',' OG ') ],
    'THR' : [ (' N  ',' CA ',' CB ',' OG1'), ],
    'ASN' : [ (' N  ',' CA ',' CB ',' CG '), (' CA ',' CB ',' CG ',' OD1'), ],
    'GLN' : [ (' N  ',' CA ',' CB ',' CG '), (' CA ',' CB ',' CG ',' CD '), (' CB ',' CG ',' CD ',' OE1'), ],
    'PHE' : [ (' N  ',' CA ',' CB ',' CG '), (' CA ',' CB ',' CG ',' CD1'), ],
    'TYR' : [ (' N  ',' CA ',' CB ',' CG '), (' CA ',' CB ',' CG ',' CD1'), ],
    'TRP' : [ (' N  ',' CA ',' CB ',' CG '), (' CA ',' CB ',' CG ',' CD1'), ],
}

# }}}

# {{{ set bondlengths and angles in AA residues, expected by builders
from builders import Constants
consts = Constants()
consts.set("TAU_QUALITY", 20.)
consts.set("N_CA", 1.458)
consts.set("CA_C", 1.525)
consts.set("CA_CB", 1.525)
consts.set("C_N", 1.329)
consts.set("C_O", 1.231)
consts.set("N_CA_C", 111.2)
consts.set("CA_C_N", 116.2)
consts.set("C_N_CA", 121.7)
consts.set("CA_C_O", 120.8)
consts.set("N_C_O", 123.0)
consts.set("N_C_CA_CB", 120.0)
consts.set("C_N_CA_CB", -123.0)
consts.set("C_CA_CB", 123.0)
consts.set("N_CA_CB", 110.5)

consts.set("ALA_CA_CB",1.521)
consts.set("ALA_CB_CA",1.521)

consts.set("VAL_CA_CB",1.540)
consts.set("VAL_CB_CA",1.540)
consts.set("VAL_CB_CG1",1.521)
consts.set("VAL_CG1_CB",1.521)
consts.set("VAL_CB_CG2",1.521)
consts.set("VAL_CG2_CB",1.521)

consts.set("LEU_CA_CB",1.530)
consts.set("LEU_CB_CA",1.530)
consts.set("LEU_CB_CG",1.530)
consts.set("LEU_CG_CB",1.530)
consts.set("LEU_CG_CD1",1.521)
consts.set("LEU_CD1_CG",1.521)
consts.set("LEU_CG_CD2",1.521)
consts.set("LEU_CD2_CG",1.521)

consts.set("ILE_CA_CB",1.540)
consts.set("ILE_CB_CA",1.540)
consts.set("ILE_CB_CG1",1.530)
consts.set("ILE_CG1_CB",1.530)
consts.set("ILE_CB_CG2",1.521)
consts.set("ILE_CG2_CB",1.521)
consts.set("ILE_CG1_CD1",1.513)
consts.set("ILE_CD1_CG1",1.513)

consts.set("PRO_N_CA",1.466)
consts.set("PRO_CA_N",1.466)
consts.set("PRO_CA_CB",1.530)
consts.set("PRO_CB_CA",1.530)
consts.set("PRO_CB_CG",1.492)
consts.set("PRO_CG_CB",1.492)
consts.set("PRO_CG_CD",1.503)
consts.set("PRO_CD_CG",1.503)
consts.set("PRO_N_CD",1.473)
consts.set("PRO_CD_N",1.473)

consts.set("ASP_CA_CB",1.530)
consts.set("ASP_CB_CA",1.530)
consts.set("ASP_CB_CG",1.516)
consts.set("ASP_CG_CB",1.516)
consts.set("ASP_CG_OD1",1.249)
consts.set("ASP_OD1_CG",1.249)
consts.set("ASP_CG_OD2",1.249)
consts.set("ASP_OD2_CG",1.249)

consts.set("GLU_CA_CB",1.530)
consts.set("GLU_CB_CA",1.530)
consts.set("GLU_CB_CG",1.520)
consts.set("GLU_CG_CB",1.520)
consts.set("GLU_CG_CD",1.516)
consts.set("GLU_CD_CG",1.516)
consts.set("GLU_CD_OE1",1.249)
consts.set("GLU_OE1_CD",1.249)
consts.set("GLU_CD_OE2",1.249)
consts.set("GLU_OE2_CD",1.249)

consts.set("ASN_CA_CB",1.530)
consts.set("ASN_CB_CA",1.530)
consts.set("ASN_CB_CG",1.516)
consts.set("ASN_CG_CB",1.516)
consts.set("ASN_CG_OD1",1.231)
consts.set("ASN_OD1_CG",1.231)
consts.set("ASN_CG_ND2",1.328)
consts.set("ASN_ND2_CG",1.328)

consts.set("GLN_CA_CB",1.530)
consts.set("GLN_CB_CA",1.530)
consts.set("GLN_CB_CG",1.520)
consts.set("GLN_CG_CB",1.520)
consts.set("GLN_CG_CD",1.516)
consts.set("GLN_CD_CG",1.516)
consts.set("GLN_CD_OE1",1.231)
consts.set("GLN_OE1_CD",1.231)
consts.set("GLN_CD_NE2",1.328)
consts.set("GLN_NE2_CD",1.328)

consts.set("LYS_CA_CB",1.530)
consts.set("LYS_CB_CA",1.530)
consts.set("LYS_CB_CG",1.520)
consts.set("LYS_CG_CB",1.520)
consts.set("LYS_CG_CD",1.520)
consts.set("LYS_CD_CG",1.520)
consts.set("LYS_CD_CE",1.520)
consts.set("LYS_CE_CD",1.520)
consts.set("LYS_CE_NZ",1.489)
consts.set("LYS_NZ_CE",1.489)

consts.set("ARG_CA_CB",1.530)
consts.set("ARG_CB_CA",1.530)
consts.set("ARG_CB_CG",1.520)
consts.set("ARG_CG_CB",1.520)
consts.set("ARG_CG_CD",1.520)
consts.set("ARG_CD_CG",1.520)
consts.set("ARG_CD_NE",1.461)
consts.set("ARG_NE_CD",1.461)
consts.set("ARG_NE_CZ",1.329)
consts.set("ARG_CZ_NE",1.329)
consts.set("ARG_CZ_NH1",1.326)
consts.set("ARG_NH1_CZ",1.326)
consts.set("ARG_CZ_NH2",1.326)
consts.set("ARG_NH2_CZ",1.326)

consts.set("SER_CA_CB",1.530)
consts.set("SER_CB_CA",1.530)
consts.set("SER_CB_OG",1.417)
consts.set("SER_OG_CB",1.417)

consts.set("THR_CA_CB",1.540)
consts.set("THR_CB_CA",1.540)
consts.set("THR_CB_OG1",1.433)
consts.set("THR_OG1_CB",1.433)
consts.set("THR_CB_CG2",1.521)
consts.set("THR_CG2_CB",1.521)

consts.set("MET_CA_CB",1.530)
consts.set("MET_CB_CA",1.530)
consts.set("MET_CB_CG",1.520)
consts.set("MET_CG_CB",1.520)
consts.set("MET_CG_SD",1.803)
consts.set("MET_SD_CG",1.803)
consts.set("MET_SD_CE",1.791)
consts.set("MET_CE_SD",1.791)

consts.set("CYS_CA_CB",1.530)
consts.set("CYS_CB_CA",1.530)
consts.set("CYS_CB_SG",1.808)
consts.set("CYS_SG_CB",1.808)

consts.set("TYR_CA_CB",1.530)
consts.set("TYR_CB_CA",1.530)
consts.set("TYR_CB_CG",1.512)
consts.set("TYR_CG_CB",1.512)
consts.set("TYR_CG_CD1",1.389)
consts.set("TYR_CD1_CG",1.389)
consts.set("TYR_CG_CD2",1.389)
consts.set("TYR_CD2_CG",1.389)
consts.set("TYR_CD1_CE1",1.382)
consts.set("TYR_CE1_CD1",1.382)
consts.set("TYR_CD2_CE2",1.382)
consts.set("TYR_CE2_CD2",1.382)
consts.set("TYR_CE1_CZ",1.378)
consts.set("TYR_CZ_CE1",1.378)
consts.set("TYR_CE2_CZ",1.378)
consts.set("TYR_CZ_CE2",1.378)
consts.set("TYR_CZ_OH",1.376)
consts.set("TYR_OH_CZ",1.376)

consts.set("TRP_CA_CB",1.530)
consts.set("TRP_CB_CA",1.530)
consts.set("TRP_CB_CG",1.498)
consts.set("TRP_CG_CB",1.498)
consts.set("TRP_CG_CD1",1.365)
consts.set("TRP_CD1_CG",1.365)
consts.set("TRP_CG_CD2",1.433)
consts.set("TRP_CD2_CG",1.433)
consts.set("TRP_CD1_NE1",1.374)
consts.set("TRP_NE1_CD1",1.374)
consts.set("TRP_CD2_CE2",1.409)
consts.set("TRP_CE2_CD2",1.409)
consts.set("TRP_NE1_CE2",1.370)
consts.set("TRP_CE2_NE1",1.370)
consts.set("TRP_CD2_CE3",1.398)
consts.set("TRP_CE3_CD2",1.398)
consts.set("TRP_CE2_CZ2",1.394)
consts.set("TRP_CZ2_CE2",1.394)
consts.set("TRP_CE3_CZ3",1.382)
consts.set("TRP_CZ3_CE3",1.382)
consts.set("TRP_CZ2_CH2",1.368)
consts.set("TRP_CH2_CZ2",1.368)
consts.set("TRP_CZ3_CH2",1.400)
consts.set("TRP_CH2_CZ3",1.400)

consts.set("HIS_CA_CB",1.530)
consts.set("HIS_CB_CA",1.530)
consts.set("HIS_CB_CG",1.497)
consts.set("HIS_CG_CB",1.497)
consts.set("HIS_CG_ND1",1.378)
consts.set("HIS_ND1_CG",1.378)
consts.set("HIS_NE2_CD2",1.319)
consts.set("HIS_CD2_NE2",1.319)
consts.set("HIS_ND1_CE1",1.345)
consts.set("HIS_CE1_ND1",1.345)
consts.set("HIS_CE1_NE2",1.319)
consts.set("HIS_NE2_CE1",1.319)

consts.set("PHE_CA_CB",1.530)
consts.set("PHE_CB_CA",1.530)
consts.set("PHE_CB_CG",1.502)
consts.set("PHE_CG_CB",1.502)
consts.set("PHE_CG_CD1",1.384)
consts.set("PHE_CD1_CG",1.384)
consts.set("PHE_CG_CD2",1.384)
consts.set("PHE_CD2_CG",1.384)
consts.set("PHE_CD1_CE1",1.382)
consts.set("PHE_CE1_CD1",1.382)
consts.set("PHE_CD2_CE2",1.382)
consts.set("PHE_CE2_CD2",1.382)
consts.set("PHE_CE1_CZ",1.382)
consts.set("PHE_CZ_CE1",1.382)
consts.set("PHE_CE2_CZ",1.382)
consts.set("PHE_CZ_CE2",1.382)




consts.set("VAL_CA_CB_CG2",110.5)
consts.set("VAL_CG2_CB_CA",110.5)
consts.set("VAL_CA_CB_CG1",110.5)
consts.set("VAL_CG1_CB_CA",110.5)

consts.set("LEU_CA_CB_CG",116.3)
consts.set("LEU_CG_CB_CA",116.3)
consts.set("LEU_CB_CG_CD1",110.7)
consts.set("LEU_CD1_CG_CB",110.7)
consts.set("LEU_CB_CG_CD2",110.7)
consts.set("LEU_CD2_CG_CB",110.7)

consts.set("ILE_CA_CB_CG1",110.4)
consts.set("ILE_CG1_CB_CA",110.4)
consts.set("ILE_CA_CB_CG2",110.5)
consts.set("ILE_CG2_CB_CA",110.5)
consts.set("ILE_CB_CG1_CD1",113.8)
consts.set("ILE_CD1_CG1_CB",113.8)

consts.set("PRO_N_CA_CB",103.0)
consts.set("PRO_CB_CA_N",103.0)
consts.set("PRO_CA_CB_CG",104.5)
consts.set("PRO_CG_CB_CA",104.5)
consts.set("PRO_CB_CG_CD",106.1)
consts.set("PRO_CD_CG_CB",106.1)
consts.set("PRO_CG_CD_N",103.2)
consts.set("PRO_N_CD_CG",103.2)
consts.set("PRO_CD_N_CA",112.0)
consts.set("PRO_CA_N_CD",112.0)


consts.set("ASP_CA_CB_CG",112.6)
consts.set("ASP_CG_CB_CA",112.6)
consts.set("ASP_CB_CG_OD1",118.4)
consts.set("ASP_OD1_CG_CB",118.4)
consts.set("ASP_CB_CG_OD2",118.4)
consts.set("ASP_OD2_CG_CB",118.4)

consts.set("GLU_CA_CB_CG",114.1)
consts.set("GLU_CG_CB_CA",114.1)
consts.set("GLU_CB_CG_CD",112.6)
consts.set("GLU_CD_CG_CB",112.6)
consts.set("GLU_CG_CD_OE1",118.4)
consts.set("GLU_OE1_CD_CG",118.4)
consts.set("GLU_CG_CD_OE2",118.4)
consts.set("GLU_OE2_CD_CG",118.4)

consts.set("ASN_CA_CB_CG",112.6)
consts.set("ASN_CG_CB_CA",112.6)
consts.set("ASN_CB_CG_OD1",120.8)
consts.set("ASN_OD1_CG_CB",120.8)
consts.set("ASN_CB_CG_ND2",116.4)
consts.set("ASN_ND2_CG_CB",116.4)

consts.set("GLN_CA_CB_CG",114.1)
consts.set("GLN_CG_CB_CA",114.1)
consts.set("GLN_CB_CG_CD",112.6)
consts.set("GLN_CD_CG_CB",112.6)
consts.set("GLN_CG_CD_OE1",120.8)
consts.set("GLN_OE1_CD_CG",120.8)
consts.set("GLN_CG_CD_NE2",116.4)
consts.set("GLN_NE2_CD_CG",116.4)

consts.set("LYS_CA_CB_CG",114.1)
consts.set("LYS_CG_CB_CA",114.1)
consts.set("LYS_CB_CG_CD",111.3)
consts.set("LYS_CD_CG_CB",111.3)
consts.set("LYS_CG_CD_CE",111.3)
consts.set("LYS_CE_CD_CG",111.3)
consts.set("LYS_CD_CE_NZ",111.9)
consts.set("LYS_NZ_CE_CD",111.9)

consts.set("ARG_CA_CB_CG",114.1)
consts.set("ARG_CG_CB_CA",114.1)
consts.set("ARG_CB_CG_CD",111.3)
consts.set("ARG_CD_CG_CB",111.3)
consts.set("ARG_CG_CD_NE",112.0)
consts.set("ARG_NE_CD_CG",112.0)
consts.set("ARG_CD_NE_CZ",124.2)
consts.set("ARG_CZ_NE_CD",124.2)
consts.set("ARG_NE_CZ_NH1",120.0)
consts.set("ARG_NH1_CZ_NE",120.0)
consts.set("ARG_NE_CZ_NH2",120.0)
consts.set("ARG_NH2_CZ_NE",120.0)

consts.set("SER_CA_CB_OG",111.1)
consts.set("SER_OG_CB_CA",111.1)

consts.set("THR_CA_CB_OG1",109.6)
consts.set("THR_OG1_CB_CA",109.6)
consts.set("THR_CA_CB_CG2",110.5)
consts.set("THR_CG2_CB_CA",110.5)

consts.set("MET_CA_CB_CG",114.1)
consts.set("MET_CG_CB_CA",114.1)
consts.set("MET_CB_CG_SD",112.7)
consts.set("MET_SD_CG_CB",112.7)
consts.set("MET_CG_SD_CE",100.9)
consts.set("MET_CE_SD_CG",100.9)

consts.set("CYS_CA_CB_SG",114.4)
consts.set("CYS_SG_CB_CA",114.4)

consts.set("TYR_CA_CB_CG",113.9)
consts.set("TYR_CG_CB_CA",113.9)
consts.set("TYR_CB_CG_CD1",120.8)
consts.set("TYR_CD1_CG_CB",120.8)
consts.set("TYR_CB_CG_CD2",120.8)
consts.set("TYR_CD2_CG_CB",120.8)
consts.set("TYR_CG_CD1_CE1",121.2)
consts.set("TYR_CE1_CD1_CG",121.2)
consts.set("TYR_CG_CD2_CE2",121.2)
consts.set("TYR_CE2_CD2_CG",121.2)
consts.set("TYR_CD1_CE1_CZ",119.6)
consts.set("TYR_CZ_CE1_CD1",119.6)
consts.set("TYR_CD2_CE2_CZ",119.6)
consts.set("TYR_CZ_CE2_CD2",119.6)
consts.set("TYR_CE1_CZ_OH",119.9)
consts.set("TYR_OH_CZ_CE1",119.9)
consts.set("TYR_CE2_CZ_OH",119.9)
consts.set("TYR_OH_CZ_CE2",119.9)

consts.set("TRP_CA_CB_CG",113.6)
consts.set("TRP_CG_CB_CA",113.6)
consts.set("TRP_CB_CG_CD1",126.9)
consts.set("TRP_CD1_CG_CB",126.9)
consts.set("TRP_CB_CG_CD2",126.8)
consts.set("TRP_CD2_CG_CB",126.8)
consts.set("TRP_CG_CD1_NE1",110.2)
consts.set("TRP_NE1_CD1_CG",110.2)
consts.set("TRP_CG_CD2_CE2",107.2)
consts.set("TRP_CE2_CD2_CG",107.2)
consts.set("TRP_CE2_CD2_CE3",118.8)
consts.set("TRP_CE3_CD2_CE2",118.8)
consts.set("TRP_CD2_CE2_CZ2",122.4)
consts.set("TRP_CZ2_CE2_CD2",122.4)
consts.set("TRP_CD2_CE3_CZ3",118.6)
consts.set("TRP_CZ3_CE3_CD2",118.6)
consts.set("TRP_CE2_CZ2_CH2",117.5)
consts.set("TRP_CH2_CZ2_CE2",117.5)

consts.set("HIS_CA_CB_CG",113.8)
consts.set("HIS_CG_CB_CA",113.8)
consts.set("HIS_CB_CG_ND1",122.7)
consts.set("HIS_ND1_CG_CB",122.7)
consts.set("HIS_CE1_NE2_CD2",106.9)
consts.set("HIS_CD2_NE2_CE1",106.9)
consts.set("HIS_CG_ND1_CE1",109.0)
consts.set("HIS_CE1_ND1_CG",109.0)
consts.set("HIS_ND1_CE1_NE2",111.7)
consts.set("HIS_NE2_CE1_ND1",111.7)

consts.set("PHE_CA_CB_CG",113.8)
consts.set("PHE_CG_CB_CA",113.8)
consts.set("PHE_CB_CG_CD1",120.7)
consts.set("PHE_CD1_CG_CB",120.7)
consts.set("PHE_CB_CG_CD2",120.7)
consts.set("PHE_CD2_CG_CB",120.7)
consts.set("PHE_CG_CD1_CE1",120.7)
consts.set("PHE_CE1_CD1_CG",120.7)
consts.set("PHE_CG_CD2_CE2",120.7)
consts.set("PHE_CE2_CD2_CG",120.7)
consts.set("PHE_CD1_CE1_CZ",120.0)
consts.set("PHE_CZ_CE1_CD1",120.0)
consts.set("PHE_CD2_CE2_CZ",120.0)
consts.set("PHE_CZ_CE2_CD2",120.0)


# }}}

# {{{ vdw radii and reductions for AA residues

## copied from RAPPER's kip.h
## Original PROBE values
## Ref: Word et al. J Mol Biol. 1999 Jan 29;285(4):1711-33
## 0.80 is too large to generate CIS peptides, CA - CA distance ideally
## is 2.771 and minimum separation of two CA - CA atoms is 2.775!
## Thus, 0.79 is the largest possible while allowing for cis peptides
PROBE_RADIUS_SCALE_FACTOR       = 0.79
ProbeRads = {}
ProbeRads['co']     = 1.65 * PROBE_RADIUS_SCALE_FACTOR
ProbeRads['C']      = 1.75 * PROBE_RADIUS_SCALE_FACTOR
ProbeRads['N']      = 1.55 * PROBE_RADIUS_SCALE_FACTOR
ProbeRads['O']      = 1.40 * PROBE_RADIUS_SCALE_FACTOR
ProbeRads['P']      = 1.80 * PROBE_RADIUS_SCALE_FACTOR
ProbeRads['S']      = 1.80 * PROBE_RADIUS_SCALE_FACTOR
ProbeRads['SE']      = 1.80 * PROBE_RADIUS_SCALE_FACTOR
ProbeRads['H']      = 1.00 * PROBE_RADIUS_SCALE_FACTOR
PROBE_HBOND_OVERLAP_MARGIN      = 0.60 * PROBE_RADIUS_SCALE_FACTOR
PROBE_DISULFIDE_OVERLAP_MARGIN  = 2.40 * PROBE_RADIUS_SCALE_FACTOR

scFactor = 1

vdwr, vdwCommon = {}, {}
vdwCommon[' N  '] = ProbeRads['N']
vdwCommon[' CA '] = ProbeRads['C']
vdwCommon[' C  '] = ProbeRads['co']
vdwCommon[' O  '] = ProbeRads['O']
vdwCommon[' CB '] = ProbeRads['C'] * scFactor

protRes = ['GLY', 'ALA', 'ARG', 'VAL', 'LYS', 'ILE', 'HIS', 'THR', 'LEU',
            'CYS', 'MET', 'MSE', 'GLN', 'ASP', 'GLU', 'SER', 'ASN', 'PRO', 'PHE', 'TYR', 'TRP']
for pr in protRes :
    vdwr[pr] = {}
    for k,v in vdwCommon.items() : vdwr[pr][k] = v
    for an in resAtoms[pr] :
        if an in (' N  ', ' CA ', ' C  ', ' O  ', ' CB ') : continue
        if ProbeRads.has_key(an[1]) : vdwr[pr][an] = ProbeRads[an[1]] * scFactor
        elif ProbeRads.has_key(an[0:2]) : vdwr[pr][an] = ProbeRads[an[0:2]] * scFactor
        else : print "VDW radius unknown for --%s--" % an; assert 1==0

vdwr['XXX'] = {}
vdwr['XXX']['C'] = ProbeRads['C']
vdwr['XXX']['N'] = ProbeRads['N']
vdwr['XXX']['O'] = ProbeRads['O']
vdwr['XXX']['P'] = ProbeRads['P']
vdwr['XXX']['S'] = ProbeRads['S']
# }}}

# {{{ atoms and connectivity in deoxy/ribonucleotides

# phosphate and sugar : names and covconn
phosphate_atomnames = [' O1P', ' O2P', ' P  ', ' O5*']
phosphate_covconn = [(' P  ',' O1P'), (' P  ',' O2P'), (' P  ',' O5*')]

deoxyribose_sugar_atomnames = [ ' C1*', ' C2*', ' C3*', ' O3*', ' C4*', ' O4*', ' C5*',]
deoxyribose_sugar_covconn = [ (' C4*', ' C5*'), (' C4*', ' O4*'), (' C1*', ' O4*'),
            (' C3*',' O3*'), (' C4*', ' C3*'), (' C2*', ' C3*'), (' C2*', ' C1*'), ]

ribose_sugar_atomnames = deoxyribose_sugar_atomnames + [' O2*']
ribose_sugar_covconn = deoxyribose_sugar_covconn + [(' C2*', ' O2*')]

P_sugar_covconn = [(' O5*', ' C5*')]

NU5 = ['A', 'T', 'C', 'G', 'U']
nuConn = {}
for nu in NU5 :
    resAtoms[' D'+nu] = phosphate_atomnames + deoxyribose_sugar_atomnames
    resAtoms['  '+nu] = phosphate_atomnames + ribose_sugar_atomnames
    nuConn[' D'+nu] = phosphate_covconn + P_sugar_covconn + deoxyribose_sugar_covconn
    nuConn['  '+nu] = phosphate_covconn + P_sugar_covconn + ribose_sugar_covconn

# purines A,G bases : names and covconn
AG_atomnames = [ ' N1 ', ' C2 ', ' N3 ', ' C4 ', ' C5 ', ' C6 ', ' N7 ', ' C8 ', ' N9 ', ]
AG_covconn = [ (' N1 ', ' C2 '), (' N3 ', ' C2 '), (' N3 ', ' C4 '),
    (' C5 ', ' C4 '), (' C5 ', ' C6 '), (' N1 ', ' C6 '),
    (' C5 ', ' N7 '), (' C8 ', ' N7 '), (' C8 ', ' N9 '), (' C4 ', ' N9 '),
]
sugarConn = [(' N9 ',' C1*')]

resAtoms['  A'] = resAtoms['  A'] + AG_atomnames + [' N6 ']
resAtoms[' DA'] = resAtoms[' DA'] + AG_atomnames + [' N6 ']
resAtoms['  G'] = resAtoms['  G'] + AG_atomnames + [' O6 ', ' N2 ']
resAtoms[' DG'] = resAtoms[' DG'] + AG_atomnames + [' O6 ', ' N2 ']
nuConn['  A'] = nuConn['  A'] + sugarConn + AG_covconn + [(' C6 ', ' N6 ')]
nuConn[' DA'] = nuConn[' DA'] + sugarConn + AG_covconn + [(' C6 ', ' N6 ')]
nuConn['  G'] = nuConn['  G'] + sugarConn + AG_covconn + [(' C2 ', ' N2 '), (' C6 ', ' O6 ')]
nuConn[' DG'] = nuConn[' DG'] + sugarConn + AG_covconn + [(' C2 ', ' N2 '), (' C6 ', ' O6 ')]

# pyrimidines T,C,U bases : names and covconn
TCU_atomnames = [ ' N1 ', ' C2 ', ' O2 ', ' N3 ', ' C4 ', ' C5 ', ' C6 ', ]
TCU_covconn = [ (' N1 ', ' C2 '), (' N3 ', ' C2 '), (' N3 ', ' C4 '),
    (' C4 ', ' C5 '), (' C6 ', ' C5 '), (' C6 ', ' N1 '), (' C2 ',' O2 ')
]
sugarConn = [(' N1 ',' C1*')]

resAtoms['  T'] = resAtoms['  T'] + TCU_atomnames + [' O4 ', ' C5M']
resAtoms[' DT'] = resAtoms[' DT'] + TCU_atomnames + [' O4 ', ' C5M']
resAtoms['  C'] = resAtoms['  C'] + TCU_atomnames + [' N4 ']
resAtoms[' DC'] = resAtoms[' DC'] + TCU_atomnames + [' N4 ']
resAtoms['  U'] = resAtoms['  U'] + TCU_atomnames + [' O4 ']
resAtoms[' DU'] = resAtoms[' DU'] + TCU_atomnames + [' O4 ']
nuConn['  T'] = nuConn['  T'] + sugarConn + TCU_covconn + [(' C5 ', ' C5M'), (' C4 ', ' O4 ')]
nuConn[' DT'] = nuConn[' DT'] + sugarConn + TCU_covconn + [(' C5 ', ' C5M'), (' C4 ', ' O4 ')]
nuConn['  C'] = nuConn['  C'] + sugarConn + TCU_covconn + [(' C4 ', ' N4 ')]
nuConn[' DC'] = nuConn[' DC'] + sugarConn + TCU_covconn + [(' C4 ', ' N4 ')]
nuConn['  U'] = nuConn['  U'] + sugarConn + TCU_covconn + [(' C4 ', ' O4 ')]
nuConn[' DU'] = nuConn[' DU'] + sugarConn + TCU_covconn + [(' C4 ', ' O4 ')]

# }}}

# {{{ bondlength and angle and torsion values in NU residues taken from dna-rna.param in xplore-nih-2.13/toppar

# distances

consts.set("NU_P_O3*", 1.607) # phosphate
consts.set("NU_P_O5*", 1.593)
consts.set("NU_P_O1P", 1.485)
consts.set("NU_P_O2P", 1.485)

consts.set("NU_O5*_C5*", 1.426) # phosphate - sugar

consts.set("NU_C5*_C4*", 1.510) # sugar
consts.set("NU_C4*_C3*", 1.524)
consts.set("NU_C3*_C2*", 1.525)
consts.set("NU_C2*_C1*", 1.528)
consts.set("NU_C1*_O4*", 1.414)
consts.set("NU_C4*_O4*", 1.453)
consts.set("NU_C3*_O3*", 1.423)
consts.set("NU_C2*_O2*", 1.413)

consts.set("NU_T_N1_C1*", 1.473) # sugar - base T
consts.set("NU_U_N1_C1*", 1.469) # sugar - base U
consts.set("NU_C_N1_C1*", 1.470) # sugar - base C
consts.set("NU_G_N9_C1*", 1.459) # sugar - base G
consts.set("NU_A_N9_C1*", 1.462) # sugar - base A

consts.set("NU_C_C2_O2", 1.240) # base C
consts.set("NU_C_C4_N4", 1.335)
consts.set("NU_C_C2_N3", 1.353)
consts.set("NU_C_N3_C4", 1.335)
consts.set("NU_C_C4_C5", 1.425)
consts.set("NU_C_C5_C6", 1.339)
consts.set("NU_C_N1_C2", 1.397)
consts.set("NU_C_N1_C6", 1.367)

consts.set("NU_T_N1_C2", 1.376) # base T
consts.set("NU_T_C2_N3", 1.373)
consts.set("NU_T_N3_C4", 1.382)
consts.set("NU_T_C4_C5", 1.445)
consts.set("NU_T_C5_C6", 1.339)
consts.set("NU_T_C6_N1", 1.378)
consts.set("NU_T_C2_O2", 1.220)
consts.set("NU_T_C4_O4", 1.228)
consts.set("NU_T_C5_C5M", 1.496)

consts.set("NU_U_C2_O2", 1.219) # base U
consts.set("NU_U_C4_O2", 1.232)
consts.set("NU_U_N1_C2", 1.381)
consts.set("NU_U_N1_C6", 1.375)
consts.set("NU_U_C2_N3", 1.373)
consts.set("NU_U_N3_C4", 1.380)
consts.set("NU_U_C4_C5", 1.431)
consts.set("NU_U_C5_C6", 1.337)
consts.set("NU_U_C4_O4", 1.228)

consts.set("NU_A_N1_C2", 1.339) # base A
consts.set("NU_A_C2_N3", 1.331)
consts.set("NU_A_N3_C4", 1.344)
consts.set("NU_A_C4_C5", 1.383)
consts.set("NU_A_C5_C6", 1.406)
consts.set("NU_A_C6_N1", 1.351)
consts.set("NU_A_C5_N7", 1.388)
consts.set("NU_A_N7_C8", 1.311)
consts.set("NU_A_C8_N9", 1.373)
consts.set("NU_A_N9_C4", 1.374)
consts.set("NU_A_C6_N6", 1.335)

consts.set("NU_G_N1_C2", 1.373) # base G
consts.set("NU_G_C2_N3", 1.323)
consts.set("NU_G_N3_C4", 1.350)
consts.set("NU_G_C4_C5", 1.379)
consts.set("NU_G_C5_C6", 1.419)
consts.set("NU_G_C6_N1", 1.391)
consts.set("NU_G_C5_N7", 1.388)
consts.set("NU_G_N7_C8", 1.305)
consts.set("NU_G_C8_N9", 1.374)
consts.set("NU_G_N9_C4", 1.375)
consts.set("NU_G_C2_N2", 1.341)
consts.set("NU_G_C6_O6", 1.237)

# angles

consts.set("NU_O1P_P_O2P", 119.600) # phosphate
consts.set("NU_O5*_P_O1P", 108.100)
consts.set("NU_O5*_P_O2P", 108.300)
consts.set("NU_O3*_P_O5*", 104.000)
consts.set("NU_O2P_P_O3*", 108.300)
consts.set("NU_O1P_P_O3*", 107.400)

consts.set("NU_O5*_C5*_C4*", 110.200) # phosphate sugar link
consts.set("NU_P_O5*_C5*", 120.900)
consts.set("NU_P_O3*_C3*", 119.700)

consts.set("NU_O4*_C4*_C3*", 105.550) # sugar
consts.set("NU_C5*_C4*_C3*", 115.100)
consts.set("NU_C5*_C4*_O4*", 109.300)
consts.set("NU_C1*_O4*_C4*", 109.650)
consts.set("NU_C4*_C3*_C2*", 102.950)
consts.set("NU_C3*_C2*_C1*", 102.100)
consts.set("NU_O4*_C1*_C2*", 106.250)
consts.set("NU_C4*_C3*_O3*", 110.400)
consts.set("NU_C2*_C3*_O3*", 111.300)
consts.set("NU_C1*_C2*_O2*", 110.600)
consts.set("NU_C3*_C2*_O2*", 113.300)

consts.set("NU_T_N1_C1*_C2*", 113.800) # sugar - base links
consts.set("NU_C_N1_C1*_C2*", 113.800)
consts.set("NU_U_N1_C1*_C2*", 113.800)
consts.set("NU_G_N9_C1*_C2*", 113.800)
consts.set("NU_A_N9_C1*_C2*", 113.800)
consts.set("NU_T_O4*_C1*_N1", 108.000)
consts.set("NU_C_O4*_C1*_N1", 108.000)
consts.set("NU_U_O4*_C1*_N1", 108.000)
consts.set("NU_A_O4*_C1*_N9", 108.000)
consts.set("NU_G_O4*_C1*_N9", 108.000)

consts.set("NU_T_C1*_N1_C2", 120) ## added by swanand, bcz 120 looks like the 'natural' value here
consts.set("NU_C_C1*_N1_C2", 120)
consts.set("NU_U_C1*_N1_C2", 120)
consts.set("NU_A_C1*_N9_C4", 120)
consts.set("NU_G_C1*_N9_C4", 120)

consts.set("NU_C_C6_N1_C2", 120.300) # base C
consts.set("NU_C_N1_C2_N3", 119.200)
consts.set("NU_C_C2_N3_C4", 119.900)
consts.set("NU_C_N3_C4_C5", 121.900)
consts.set("NU_C_C4_C5_C6", 117.400)
consts.set("NU_C_C5_C6_N1", 121.000)
consts.set("NU_C_N1_C2_O2", 118.900)
consts.set("NU_C_N3_C2_O2", 121.900)
consts.set("NU_C_N3_C4_N4", 118.000)
consts.set("NU_C_C5_C4_N4", 120.200)
consts.set("NU_C_C6_N1_C1*", 120.800)
consts.set("NU_C_C2_N1_C1*", 118.800)

consts.set("NU_T_C6_N1_C2", 121.300) # base T
consts.set("NU_T_N1_C2_N3", 114.600)
consts.set("NU_T_C2_N3_C4", 127.200)
consts.set("NU_T_N3_C4_C5", 115.200)
consts.set("NU_T_C4_C5_C6", 118.000)
consts.set("NU_T_C5_C6_N1", 123.700)
consts.set("NU_T_N1_C2_O2", 123.100)
consts.set("NU_T_N3_C2_O2", 122.300)
consts.set("NU_T_N3_C4_O4", 119.900)
consts.set("NU_T_C5_C4_O4", 124.900)
consts.set("NU_T_C4_C5_C5M", 119.000)
consts.set("NU_T_C6_C5_C5M", 122.900)
consts.set("NU_T_C6_N1_C1*", 120.400)
consts.set("NU_T_C2_N1_C1*", 118.200)

consts.set("NU_U_C6_N1_C2", 121.000) # base U
consts.set("NU_U_N1_C2_N3", 114.900)
consts.set("NU_U_C2_N3_C4", 127.000)
consts.set("NU_U_N3_C4_C5", 114.600)
consts.set("NU_U_C4_C5_C6", 119.700)
consts.set("NU_U_C5_C6_N1", 122.700)
consts.set("NU_U_N1_C2_O2", 122.800)
consts.set("NU_U_N3_C2_O2", 122.200)
consts.set("NU_U_N3_C4_O4", 119.400)
consts.set("NU_U_C5_C4_O4", 125.900)
consts.set("NU_U_C6_N1_C1*", 121.200)
consts.set("NU_U_C2_N1_C1*", 117.700)

consts.set("NU_A_C6_N1_C2", 118.600) # base A
consts.set("NU_A_N1_C2_N3", 129.300)
consts.set("NU_A_C2_N3_C4", 110.600)
consts.set("NU_A_N3_C4_C5", 126.800)
consts.set("NU_A_C4_C5_C6", 117.000)
consts.set("NU_A_C5_C6_N1", 117.700)
consts.set("NU_A_C4_C5_N7", 110.700)
consts.set("NU_A_C5_N7_C8", 103.900)
consts.set("NU_A_N7_C8_N9", 113.800)
consts.set("NU_A_C8_N9_C4", 105.800)
consts.set("NU_A_N9_C4_C5", 105.800)
consts.set("NU_A_N3_C4_N9", 127.400)
consts.set("NU_A_C6_C5_N7", 132.300)
consts.set("NU_A_N1_C6_N6", 118.600)
consts.set("NU_A_C5_C6_N6", 123.700)
consts.set("NU_A_C8_N9_C1*", 127.700)
consts.set("NU_A_C4_N9_C1*", 126.300)

consts.set("NU_G_C6_N1_C2", 125.100) # base G
consts.set("NU_G_N1_C2_N3", 123.900)
consts.set("NU_G_C2_N3_C4", 111.900)
consts.set("NU_G_N3_C4_C5", 128.600)
consts.set("NU_G_C4_C5_C6", 118.800)
consts.set("NU_G_C5_C6_N1", 111.500)
consts.set("NU_G_C4_C5_N7", 110.800)
consts.set("NU_G_C5_N7_C8", 104.300)
consts.set("NU_G_N7_C8_N9", 113.100)
consts.set("NU_G_C8_N9_C4", 106.400)
consts.set("NU_G_N9_C4_C5", 105.400)
consts.set("NU_G_N3_C4_N9", 126.000)
consts.set("NU_G_C6_C5_N7", 130.400)
consts.set("NU_G_N1_C2_N2", 116.200)
consts.set("NU_G_N3_C2_N2", 119.900)
consts.set("NU_G_N1_C6_O6", 119.900)
consts.set("NU_G_C5_C6_O6", 128.600)
consts.set("NU_G_C8_N9_C1*", 127.000)
consts.set("NU_G_C4_N9_C1*", 126.500)

# dihedrals

consts.set("NU_O3*_O5*_P_O2P", 114.500) # phosphate
consts.set("NU_O3*_O5*_P_O1P", -114.500)

consts.set("NU_C3E_C5*_C4*_C3*_O3*",  81.1) # C3'-endo within sugar ring
consts.set("NU_C3E_O4*_C4*_C3*_O3*", 201.8)
consts.set("NU_C3E_O4*_C1*_C2*_C3*", 335.4)
consts.set("NU_C3E_C1*_C2*_C3*_C4*",  35.9)
consts.set("NU_C3E_C2*_C3*_C4*_O4*", 324.7)
consts.set("NU_C3E_C3*_C4*_O4*_C1*",  20.5)
consts.set("NU_C3E_C4*_O4*_C1*_C2*",   2.8)
consts.set("NU_C3E_C5*_C4*_C3*_C2*", 204.0)
consts.set("NU_C3E_O3*_C3*_C2*_O2*",  44.3)

consts.set("NU_C3E_T_C4*_O4*_C1*_N1", 241.4) # C3'-endo sugar-base
consts.set("NU_C3E_C_C4*_O4*_C1*_N1", 241.4)
consts.set("NU_C3E_U_C4*_O4*_C1*_N1", 241.4)
consts.set("NU_C3E_G_C4*_O4*_C1*_N9", 241.4)
consts.set("NU_C3E_A_C4*_O4*_C1*_N9", 241.4)

consts.set("NU_C3E_T_O4*_C1*_N1_C2", 195.7) # C3'-endo chi dihedrals
consts.set("NU_C3E_C_O4*_C1*_N1_C2", 195.7)
consts.set("NU_C3E_U_O4*_C1*_N1_C2", 195.7)
consts.set("NU_C3E_A_O4*_C1*_N9_C4", 193.3)
consts.set("NU_C3E_G_O4*_C1*_N9_C4", 193.3)


consts.set("NU_C2E_C5*_C4*_C3*_O3*", 146.25) # C2'-endo within sugar ring
consts.set("NU_C2E_O4*_C1*_C2*_C3*",  34.0)
consts.set("NU_C2E_O4*_C4*_C3*_O3*", 267.0)
consts.set("NU_C2E_C1*_C2*_C3*_C4*", 325.75)
consts.set("NU_C2E_C2*_C3*_C4*_O4*",  23.4)
consts.set("NU_C2E_C3*_C4*_O4*_C1*", 357.7)
consts.set("NU_C2E_C4*_O4*_C1*_C2*", 340.0)
consts.set("NU_C2E_C5*_C4*_C3*_C2*", 262.7)
consts.set("NU_C2E_O3*_C3*_C2*_O2*", 319.7)

consts.set("NU_C2E_T_C4*_O4*_C1*_N1", 217.15) # C2'-endo sugar-base
consts.set("NU_C2E_C_C4*_O4*_C1*_N1", 217.15)
consts.set("NU_C2E_U_C4*_O4*_C1*_N1", 217.15)
consts.set("NU_C2E_G_C4*_O4*_C1*_N9", 217.15)
consts.set("NU_C2E_A_C4*_O4*_C1*_N9", 217.15)

consts.set("NU_C2E_T_O4*_C1*_N1_C2", 229.8) # C2'-endo chi dihedrals
consts.set("NU_C2E_C_O4*_C1*_N1_C2", 229.8)
consts.set("NU_C2E_U_O4*_C1*_N1_C2", 229.8)
consts.set("NU_C2E_G_O4*_C1*_N9_C4", 237.0)
consts.set("NU_C2E_A_O4*_C1*_N9_C4", 237.0)

# }}}

#{{{ spacegroup notations table

long2shortHM = {}
long2shortHM["P121"] = "P2"
long2shortHM["P1211"] = "P21"
long2shortHM["C121"] = "C2"
long2shortHM["A121"] = "A2"
long2shortHM["I121"] = "I2"
long2shortHM["P1m1"] = "PM"
long2shortHM["P1c1"] = "PC"
long2shortHM["P1n1"] = "PN"
long2shortHM["P1a1"] = "PA"
long2shortHM["C1m1"] = "CM"
long2shortHM["A1m1"] = "AM"
long2shortHM["I1m1"] = "IM"
long2shortHM["C1c1"] = "CC"
long2shortHM["A1n1"] = "AN"
long2shortHM["I1a1"] = "IA"
long2shortHM["A1a1"] = "AA"
long2shortHM["C1n1"] = "CN"
long2shortHM["I1c1"] = "IC"
long2shortHM["P12/m1"] = "P2/M"
long2shortHM["P121/m1"] = "P21/M"
long2shortHM["C12/m1"] = "C2/M"
long2shortHM["A12/m1"] = "A2/M"
long2shortHM["I12/m1"] = "I2/M"
long2shortHM["P12/c1"] = "P2/C"
long2shortHM["P12/n1"] = "P2/N"
long2shortHM["P12/a1"] = "P2/A"
long2shortHM["P121/c1"] = "P21/C"
long2shortHM["P121/n1"] = "P21/N"
long2shortHM["P121/a1"] = "P21/A"
long2shortHM["C12/c1"] = "C2/C"
long2shortHM["A12/n1"] = "A2/N"
long2shortHM["I12/a1"] = "I2/A"
long2shortHM["A12/a1"] = "A2/A"
long2shortHM["C12/n1"] = "C2/N"
long2shortHM["I12/c1"] = "I2/C"




sgtable = {}
sgtable["P1"] = ["P 1","P1","1"]
sgtable["P-1"] = ["P -1","P-1","2"]
sgtable["P2"] = ["P 2","P2","2"]
sgtable["P21"] = ["P 21","P2(1)","2"]
sgtable["C2"] = ["C 2","C2","4"]
sgtable["PM"] = ["P M","Pm","2"]
sgtable["PC"] = ["P C","Pc","2"]
sgtable["CM"] = ["C M","Cm","4"]
sgtable["CC"] = ["C C","Cc","4"]
sgtable["P2/M"] = ["P 2/M","P2/m","4"]
sgtable["P21/M"] = ["P 21/M","P2(1)/m","4"]
sgtable["C2/M"] = ["C 2/M","C2/m","8"]
sgtable["P2/C"] = ["P 2/C","P2/c","4"]
sgtable["P21/C"] = ["P 21/C","P2(1)/c","4"]
sgtable["C2/C"] = ["C 2/C","C2/c","8"]
sgtable["P222"] = ["P 2 2 2","P222","4"]
sgtable["P2221"] = ["P 2 2 21","P222(1)","4"]
sgtable["P21212"] = ["P 21 21 2","P2(1)2(1)2","4"]
sgtable["P212121"] = ["P 21 21 21","P2(1)2(1)2(1)","4"]
sgtable["C2221"] = ["C 2 2 21","C222(1)","8"]
sgtable["C222"] = ["C 2 2 2","C222","8"]
sgtable["F222"] = ["F 2 2 2","F222","16"]
sgtable["I222"] = ["I 2 2 2","I222","8"]
sgtable["I212121"] = ["I 21 21 21","I2(1)2(1)2(1)","8"]
sgtable["PMM2"] = ["P M M 2","Pmm2","4"]
sgtable["PMC21"] = ["P M C 21","Pmc2(1)","4"]
sgtable["PCC2"] = ["P C C 2","Pcc2","4"]
sgtable["PMA2"] = ["P M A 2","Pma2","4"]
sgtable["PCA21"] = ["P C A 21","Pca2(1)","4"]
sgtable["PNC2"] = ["P N C 2","Pnc2","4"]
sgtable["PMN21"] = ["P M N 21","Pmn2(1)","4"]
sgtable["PBA2"] = ["P B A 2","Pba2","4"]
sgtable["PNA21"] = ["P N A 21","Pna2(1)","4"]
sgtable["PNN2"] = ["P N N 2","Pnn2","4"]
sgtable["CMM2"] = ["C M M 2","Cmm2","8"]
sgtable["CMC21"] = ["C M C 21","Cmc2(1)","8"]
sgtable["CCC2"] = ["C C C 2","Ccc2","8"]
sgtable["AMM2"] = ["A M M 2","Amm2","8"]
sgtable["ABM2"] = ["A B M 2","Abm2","8"]
sgtable["AMA2"] = ["A M A 2","Ama2","8"]
sgtable["ABA2"] = ["A B A 2","Aba2","8"]
sgtable["FMM2"] = ["F M M 2","Fmm2","16"]
sgtable["FDD2"] = ["F D D 2","Fdd2","16"]
sgtable["IMM2"] = ["I M M 2","Imm2","8"]
sgtable["IBA2"] = ["I B A 2","Iba2","8"]
sgtable["IMA2"] = ["I M A 2","Ima2","8"]
sgtable["PMMM"] = ["P M M M","Pmmm","8"]
sgtable["PNNN"] = ["P N N N","Pnnn","8"]
sgtable["PCCM"] = ["P C C M","Pccm","8"]
sgtable["PBAN"] = ["P B A N","Pban","8"]
sgtable["PMMA"] = ["P M M A","Pmma","8"]
sgtable["PNNA"] = ["P N N A","Pnna","8"]
sgtable["PMNA"] = ["P M N A","Pmna","8"]
sgtable["PCCA"] = ["P C C A","Pcca","8"]
sgtable["PBAM"] = ["P B A M","Pbam","8"]
sgtable["PCCN"] = ["P C C N","Pccn","8"]
sgtable["PBCM"] = ["P B C M","Pbcm","8"]
sgtable["PNNM"] = ["P N N M","Pnnm","8"]
sgtable["PMMN"] = ["P M M N","Pmmn","8"]
sgtable["PBCN"] = ["P B C N","Pbcn","8"]
sgtable["PBCA"] = ["P B C A","Pbca","8"]
sgtable["PNMA"] = ["P N M A","Pnma","8"]
sgtable["CMCM"] = ["C M C M","Cmcm","16"]
sgtable["CMCA"] = ["C M C A","Cmca","16"]
sgtable["CMMM"] = ["C M M M","Cmmm","16"]
sgtable["CCCM"] = ["C C C M","Cccm","16"]
sgtable["CMMA"] = ["C M M A","Cmma","16"]
sgtable["CCCA"] = ["C C C A","Ccca","16"]
sgtable["FMMM"] = ["F M M M","Fmmm","32"]
sgtable["FDDD"] = ["F D D D","Fddd","32"]
sgtable["IMMM"] = ["I M M M","Immm","16"]
sgtable["IBAM"] = ["I B A M","Ibam","16"]
sgtable["IBCA"] = ["I B C A","Ibca","16"]
sgtable["IMMA"] = ["I M M A","Imma","16"]
sgtable["P4"] = ["P 4","P4","4"]
sgtable["P41"] = ["P 41","P4(1)","4"]
sgtable["P42"] = ["P 42","P4(2)","4"]
sgtable["P43"] = ["P 43","P4(3)","4"]
sgtable["I4"] = ["I 4","I4","8"]
sgtable["I41"] = ["I 41","I4(1)","8"]
sgtable["P-4"] = ["P -4","P-4","4"]
sgtable["I-4"] = ["I -4","I-4","8"]
sgtable["P4/M"] = ["P 4/M","P4/m","8"]
sgtable["P42/M"] = ["P 42/M","P4(2)/m","8"]
sgtable["P4/N"] = ["P 4/N","P4/n","8"]
sgtable["P42/N"] = ["P 42/N","P4(2)/n","8"]
sgtable["I4/M"] = ["I 4/M","I4/m","16"]
sgtable["I41/A"] = ["I 41/A","I4(1)/a","16"]
sgtable["P422"] = ["P 4 2 2","P422","8"]
sgtable["P4212"] = ["P 4 21 2","P42(1)2","8"]
sgtable["P4122"] = ["P 41 2 2","P4(1)22","8"]
sgtable["P41212"] = ["P 41 21 2","P4(1)2(1)2","8"]
sgtable["P4222"] = ["P 42 2 2","P4(2)22","8"]
sgtable["P42212"] = ["P 42 21 2","P4(2)2(1)2","8"]
sgtable["P4322"] = ["P 43 2 2","P4(3)22","8"]
sgtable["P43212"] = ["P 43 21 2","P4(3)2(1)2","8"]
sgtable["I422"] = ["I 4 2 2","I422","16"]
sgtable["I4122"] = ["I 41 2 2","I4(1)22","16"]
sgtable["P4MM"] = ["P 4 M M","P4mm","8"]
sgtable["P4BM"] = ["P 4 B M","P4bm","8"]
sgtable["P42CM"] = ["P 42 C M","P4(2)cm","8"]
sgtable["P42NM"] = ["P 42 N M","P4(2)nm","8"]
sgtable["P4CC"] = ["P 4 C C","P4cc","8"]
sgtable["P4NC"] = ["P 4 N C","P4nc","8"]
sgtable["P42MC"] = ["P 42 M C","P4(2)mc","8"]
sgtable["P42BC"] = ["P 42 B C","P4(2)bc","8"]
sgtable["I4MM"] = ["I 4 M M","I4mm","16"]
sgtable["I4CM"] = ["I 4 C M","I4cm","16"]
sgtable["I41MD"] = ["I 41 M D","I4(1)md","16"]
sgtable["I41CD"] = ["I 41 C D","I4(1)cd","16"]
sgtable["P-42M"] = ["P -4 2 M","P-42m","8"]
sgtable["P-42C"] = ["P -4 2 C","P-42c","8"]
sgtable["P-421M"] = ["P -4 21 M","P-42(1)m","8"]
sgtable["P-421C"] = ["P -4 21 C","P-42(1)c","8"]
sgtable["P-4M2"] = ["P -4 M 2","P-4m2","8"]
sgtable["P-4C2"] = ["P -4 C 2","P-4c2","8"]
sgtable["P-4B2"] = ["P -4 B 2","P-4b2","8"]
sgtable["P-4N2"] = ["P -4 N 2","P-4n2","8"]
sgtable["I-4M2"] = ["I -4 M 2","I-4m2","16"]
sgtable["I-4C2"] = ["I -4 C 2","I-4c2","16"]
sgtable["I-42M"] = ["I -4 2 M","I-42m","16"]
sgtable["I-42D"] = ["I -4 2 D","I-42d","16"]
sgtable["P4/MMM"] = ["P 4/M M M","P4/mmm","16"]
sgtable["P4/MCC"] = ["P 4/M C C","P4/mcc","16"]
sgtable["P4/NBM"] = ["P 4/N B M","P4/nbm","16"]
sgtable["P4/NNC"] = ["P 4/N N C","P4/nnc","16"]
sgtable["P4/MBM"] = ["P 4/M B M","P4/mbm","16"]
sgtable["P4/MNC"] = ["P 4/M N C","P4/mnc","16"]
sgtable["P4/NMM"] = ["P 4/N M M","P4/nmm","16"]
sgtable["P4/NCC"] = ["P 4/N C C","P4/ncc","16"]
sgtable["P42/MMC"] = ["P 42/M M C","P4(2)/mmc","16"]
sgtable["P42/MCM"] = ["P 42/M C M","P4(2)/mcm","16"]
sgtable["P42/NBC"] = ["P 42/N B C","P4(2)/nbc","16"]
sgtable["P42/NNM"] = ["P 42/N N M","P4(2)/nnm","16"]
sgtable["P42/MBC"] = ["P 42/M B C","P4(2)/mbc","16"]
sgtable["P42/MNM"] = ["P 42/M N M","P4(2)/mnm","16"]
sgtable["P42/NMC"] = ["P 42/N M C","P4(2)/nmc","16"]
sgtable["P42/NCM"] = ["P 42/N C M","P4(2)/ncm","16"]
sgtable["I4/MMM"] = ["I 4/M M M","I4/mmm","32"]
sgtable["I4/MCM"] = ["I 4/M C M","I4/mcm","32"]
sgtable["I41/AMD"] = ["I 41/A M D","I4(1)/amd","32"]
sgtable["I41/ACD"] = ["I 41/A C D","I4(1)/acd","32"]
sgtable["P3"] = ["P 3","P3","3"]
sgtable["P31"] = ["P 31","P3(1)","3"]
sgtable["P32"] = ["P 32","P3(2)","3"]
sgtable["R3"] = ["R 3","R3","9"]
sgtable["P-3"] = ["P -3","P-3","6"]
sgtable["R-3"] = ["R -3","R-3","18"]
sgtable["P312"] = ["P 3 1 2","P312","6"]
sgtable["P321"] = ["P 3 2 1","P321","6"]
sgtable["P3112"] = ["P 31 1 2","P3(1)12","6"]
sgtable["P3121"] = ["P 31 2 1","P3(1)21","6"]
sgtable["P3212"] = ["P 32 1 2","P3(2)12","6"]
sgtable["P3221"] = ["P 32 2 1","P3(2)21","6"]
sgtable["R32"] = ["R 3 2","R32","18"]
sgtable["P3M1"] = ["P 3 M 1","P3m1","6"]
sgtable["P31M"] = ["P 3 1 M","P31m","6"]
sgtable["P3C1"] = ["P 3 C 1","P3c1","6"]
sgtable["P31C"] = ["P 3 1 C","P31c","6"]
sgtable["R3M"] = ["R 3 M","R3m","18"]
sgtable["R3C"] = ["R 3 C","R3c","18"]
sgtable["P-31M"] = ["P -3 1 M","P-31m","12"]
sgtable["P-31C"] = ["P -3 1 C","P-31c","12"]
sgtable["P-3M1"] = ["P -3 M 1","P-3m1","12"]
sgtable["P-3C1"] = ["P -3 C 1","P-3c1","12"]
sgtable["R-3M"] = ["R -3 M","R-3m","36"]
sgtable["R-3C"] = ["R -3 C","R-3c","36"]
sgtable["P6"] = ["P 6","P6","6"]
sgtable["P61"] = ["P 61","P6(1)","6"]
sgtable["P65"] = ["P 65","P6(5)","6"]
sgtable["P62"] = ["P 62","P6(2)","6"]
sgtable["P64"] = ["P 64","P6(4)","6"]
sgtable["P63"] = ["P 63","P6(3)","6"]
sgtable["P-6"] = ["P -6","P-6","6"]
sgtable["P6/M"] = ["P 6/M","P6/m","12"]
sgtable["P63/M"] = ["P 63/M","P6(3)/m","12"]
sgtable["P622"] = ["P 6 2 2","P622","12"]
sgtable["P6122"] = ["P 61 2 2","P6(1)22","12"]
sgtable["P6522"] = ["P 65 2 2","P6(5)22","12"]
sgtable["P6222"] = ["P 62 2 2","P6(2)22","12"]
sgtable["P6422"] = ["P 64 2 2","P6(4)22","12"]
sgtable["P6322"] = ["P 63 2 2","P6(3)22","12"]
sgtable["P6MM"] = ["P 6 M M","P6mm","12"]
sgtable["P6CC"] = ["P 6 C C","P6cc","12"]
sgtable["P63CM"] = ["P 63 C M","P6(3)cm","12"]
sgtable["P63MC"] = ["P 63 M C","P6(3)mc","12"]
sgtable["P-6M2"] = ["P -6 M 2","P-6m2","12"]
sgtable["P-6C2"] = ["P -6 C 2","P-6c2","12"]
sgtable["P-62M"] = ["P -6 2 M","P-62m","12"]
sgtable["P-62C"] = ["P -6 2 C","P-62c","12"]
sgtable["P6/MMM"] = ["P 6/M M M","P6/mmm","24"]
sgtable["P6/MCC"] = ["P 6/M C C","P6/mcc","24"]
sgtable["P63/MCM"] = ["P 63/M C M","P6(3)/mcm","24"]
sgtable["P63/MMC"] = ["P 63/M M C","P6(3)/mmc","24"]
sgtable["P23"] = ["P 2 3","P23","12"]
sgtable["F23"] = ["F 2 3","F23","48"]
sgtable["I23"] = ["I 2 3","I23","24"]
sgtable["P213"] = ["P 21 3","P2(1)3","12"]
sgtable["I213"] = ["I 21 3","I2(1)3","24"]
sgtable["PM-3"] = ["P M -3","Pm-3","24"]
sgtable["PN-3"] = ["P N -3","Pn-3","24"]
sgtable["FM-3"] = ["F M -3","Fm-3","96"]
sgtable["FD-3"] = ["F D -3","Fd-3","96"]
sgtable["IM-3"] = ["I M -3","Im-3","48"]
sgtable["PA-3"] = ["P A -3","Pa-3","24"]
sgtable["IA-3"] = ["I A -3","Ia-3","48"]
sgtable["P432"] = ["P 4 3 2","P432","24"]
sgtable["P4232"] = ["P 42 3 2","P4(2)32","24"]
sgtable["F432"] = ["F 4 3 2","F432","96"]
sgtable["F4132"] = ["F 41 3 2","F4(1)32","96"]
sgtable["I432"] = ["I 4 3 2","I432","48"]
sgtable["P4332"] = ["P 43 3 2","P4(3)32","24"]
sgtable["P4132"] = ["P 41 3 2","P4(1)32","24"]
sgtable["I4132"] = ["I 41 3 2","I4(1)32","48"]
sgtable["P-43M"] = ["P -4 3 M","P-43m","24"]
sgtable["F4-3M"] = ["F 4 -3 M","F4-3m","96"]
sgtable["I-43M"] = ["I -4 3 M","I-43m","48"]
sgtable["P-43N"] = ["P -4 3 N","P-43n","24"]
sgtable["F-43C"] = ["F -4 3 C","F-43c","96"]
sgtable["I-43D"] = ["I -4 3 D","I-43d","48"]
sgtable["PM-3M"] = ["P M -3 M","Pm-3m","48"]
sgtable["PN-3N"] = ["P N -3 N","Pn-3n","48"]
sgtable["PM-3N"] = ["P M -3 N","Pm-3n","48"]
sgtable["PN-3M"] = ["P N -3 M","Pn-3m","48"]
sgtable["FM-3M"] = ["F M -3 M","Fm-3m","192"]
sgtable["FM-3C"] = ["F M -3 C","Fm-3c","192"]
sgtable["FD-3M"] = ["F D -3 M","Fd-3m","192"]
sgtable["FD-3C"] = ["F D -3 C","Fd-3c","192"]
sgtable["IM-3M"] = ["I M -3 M","Im-3m","96"]
sgtable["IA-3D"] = ["I A -3 D","Ia-3d","96"]
#}}}

import numpy as np
import pyblock.pd_utils
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import pyblock

pd.set_option("display.precision", 15)

def RemoveEmptyTokens(tokens):
    newTokens = []

    for token in tokens:
        if token != '' and token != '\n' and token != '\r' and token != '\r\n':
            newTokens.append(token)

    return newTokens

colNameToOrdinal = {
    "Cycle No.":0,
    "E(Proj_CCSD)":1,
    "E(Proj_F12)":2,
    "E(Proj_Tot)":3,
    "n(E)":4,
    "n(E_F12)":5,
    "E(S)":6,
    "n(Excips)":7,
    "n(Excips_D0)":8
}
ordinalToColName = dict((reversed(item) for item in colNameToOrdinal.items()))


def ReadSCESSOutFile(filePath):
    with open(filePath, 'r', encoding='utf-8') as scessFile:
        fileLine = scessFile.readline()
        while (fileLine.rstrip().find("| Cycle No.") == -1) and fileLine != '':
            fileLine = scessFile.readline();

        if fileLine == "":
            raise Exception("Cannot locate the start of the data block in the SCESS file")
        
        dataSpec = {}

        for key in colNameToOrdinal.keys():
            dataSpec[key] = []

        fileLine = scessFile.readline()
        while fileLine != '\n':

            tokens = RemoveEmptyTokens(fileLine.split(' '));

            for i, tkn in enumerate(tokens):
                if i == 1 or i ==2 or i==3 or i==6:
                    dataSpec[ordinalToColName[i]].append(float(tkn))
                else:
                    dataSpec[ordinalToColName[i]].append(int(tkn))

            fileLine = scessFile.readline()

        return pd.DataFrame(dataSpec)
    
def BlockAnalyseData(data):
    (data_length, reblock_data, covariance) = pyblock.pd_utils.reblock(data)
    print(reblock_data)

def PlotProjectedEnergy(df):
    plt.plot(h2oData["Cycle No."], h2oData["E(Proj_Tot)"])
    plt.show()


h2Data = ReadSCESSOutFile('H2 ccsd-qmc-mp2-f12 sto-3g cc-pVDZ.txt')
h2oData = ReadSCESSOutFile('H2O ccsd-qmc-mp2-f12 sto-3g cc-pVDZ.txt')
liHData = ReadSCESSOutFile('LiH ccsd-qmc-mp2-f12 sto-3g cc-pVDZ.txt')
nh3Data = ReadSCESSOutFile('NH3 ccsd-qmc-mp2-f12 sto-3g cc-pVDZ.txt')


BlockAnalyseData(h2Data["E(Proj_Tot)"])
BlockAnalyseData(h2oData["E(Proj_Tot)"])
BlockAnalyseData(liHData["E(Proj_Tot)"])
BlockAnalyseData(nh3Data["E(Proj_Tot)"])


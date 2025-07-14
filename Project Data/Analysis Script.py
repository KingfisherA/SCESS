import numpy as np
import pyblock.pd_utils
import pyblock.plot
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyblock

plt.rcParams["font.family"] = "Clear Sans"

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
    return reblock_data

def PlotProjectedEnergy(df, ax, colour, title):
    ax.title.set_text(title)
    ax.set_ylabel("$E_{Proj}$")
    ax.set_xlabel("Cycle No.")
    ax.plot(df["Cycle No."], df["E(Proj_Tot)"], color=colour)

#Plot a log log graph of n(Excips)/n(D0) against n(Excips)
def PlotShoulderPlot(df, ax, colour, title):
    ratioArray = np.array(df["n(Excips)"]) / np.array(df["n(Excips_D0)"]) 
    ax.title.set_text(title)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r"$N_{excip}$")
    ax.set_ylabel(r"$N_{excip} / N_0$")
    ax.plot(df["n(Excips)"], ratioArray, color=colour)

def PlotReblocking(block_info, ax, colour, title):
    ax.title.set_text(title)
    data_sets = block_info.columns.get_level_values(0).unique()

    for (i, col) in enumerate(data_sets):
        # There should only be (at most) one non-null value for optimal block.
        opt = block_info[block_info[(col,'optimal block')] != ''].index.values
        if opt:
            opt = opt[0]

        #Standard Error
        std_err = block_info[(col, 'standard error')]
        if 'standard error error' in block_info[col]:
            std_err_err = block_info[(col, 'standard error error')]
        else:
            std_err_err = 0*std_err
        line, caps, bars = ax.errorbar(block_info.index, std_err, yerr=std_err_err, 
                                       marker='o', label=col, color=colour, capsize=5)
        caps[0].set_marker('_')
        caps[1].set_marker('_')

        if opt:
            ax.annotate('', (opt, std_err[opt]-std_err_err[opt]),
                         xytext=(0, -20), textcoords='offset points',
                         arrowprops=dict(
                             arrowstyle="->", #color=line[0].get_color()
                             linewidth=1.2*line.get_linewidth(),
                       ),)

        #ax.legend(loc=2)
        ax.set_ylabel('S.E.')
        ax.set_xlabel('Block Size')

h2Data = ReadSCESSOutFile('H2 ccsd-qmc-mp2-f12 sto-3g cc-pVDZ.txt')
h2oData = ReadSCESSOutFile('H2O ccsd-qmc-mp2-f12 sto-3g cc-pVDZ.txt')
liHData = ReadSCESSOutFile('LiH ccsd-qmc-mp2-f12 sto-3g cc-pVDZ.txt')
nh3Data = ReadSCESSOutFile('NH3 ccsd-qmc-mp2-f12 sto-3g cc-pVDZ.txt')


h2BlockData = BlockAnalyseData(h2Data["E(Proj_Tot)"])
h2oBlockData = BlockAnalyseData(h2oData["E(Proj_Tot)"])
liHBlockData = BlockAnalyseData(liHData["E(Proj_Tot)"])
nh3BlockData = BlockAnalyseData(nh3Data["E(Proj_Tot)"])

molToStyle = {"H_2O":{"data":h2oData, "blockData":h2oBlockData, "colour":"red"},
              "H_2":{"data":h2Data,"blockData":h2BlockData, "colour":"blue"},
              "LiH":{"data":liHData,"blockData":liHBlockData, "colour":"green"},
              "NH_3":{"data":nh3Data,"blockData":nh3BlockData, "colour":"teal"}}


shoulderFig = plt.figure(figsize=(10,6))
projEFig = plt.figure(figsize=(10,6))
blockingFig = plt.figure(figsize=(10,6))
gs = gridspec.GridSpec(2,2)
gs.update(hspace=0.4)

for i in range(0,2):
    for j in range(0,2):
        idx = i * 2 + j
        key = list(molToStyle.keys())[idx]
        value = molToStyle[key]
        PlotShoulderPlot(value["data"], shoulderFig.add_subplot(gs[i,j]), value["colour"], f"$\mathrm{{{key}}}$")
        PlotProjectedEnergy(value["data"], projEFig.add_subplot(gs[i,j]), value["colour"], f"$\mathrm{{{key}}}$")
        PlotReblocking(value["blockData"], blockingFig.add_subplot(gs[i,j]), value["colour"], f"$\mathrm{{{key}}}$")

shoulderFig.show()
projEFig.show()
blockingFig.show()
plt.show()

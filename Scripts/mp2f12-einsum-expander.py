import numpy as np
import sys

nOcc = "nOcc"
nOBS = "nOBS"
nCABS = "nCABS"

indices ={
    "a":["nOcc", "nVirt"],
    "b":["nOcc", "nVirt"],
    "c":["nOcc", "nVirt"],
    "d":["nOcc", "nVirt"],

    "i":["0", "nOcc"],
    "j":["0", "nOcc"],
    "m":["0", "nOcc"],
    "n":["0", "nOcc"],
    "k":["0", "nOcc"],
    "l":["0", "nOcc"],

    "p":["0", "nOBS"],
    "q":["0", "nOBS"],
    "r":["0", "nOBS"],
    "s":["0", "nOBS"],

    "x":["nOBS", "nCABS"],
    "y":["nOBS", "nCABS"],
    "w":["nOBS", "nCABS"],
    "z":["nOBS", "nCABS"]
}

LINE_B=0
LINE_TMP=1
LINE_UNKNOWN=-1

def ParseEinsum(einsumStr, tensor1, tensor2, lineType):
    tokens = []

    tokens.append(einsumStr[:einsumStr.find(',')].strip())
    tokens.append(einsumStr[einsumStr.find(',')+1:einsumStr.find('->')].strip())
    tokens.append(einsumStr[einsumStr.find('->')+2:].strip())

    #Slice 

    def ParseSlice(token):
        slice=[[],[]]
        for c in token:
            idx = indices[c]
            slice[0].append(idx[0])
            slice[1].append(idx[1])
        return slice

    slice1=ParseSlice(tokens[0])
    slice2=ParseSlice(tokens[1])

    #Contract
    def ParseContract(token1, token2):
        contractions=[]
        for i, c in enumerate(token1):
            j=token2.find(c)
            if j != -1:
                contractions.append([i,j])

        return contractions

    contractions = ParseContract(tokens[0], tokens[1])

    #Shuffle
    def ParseShuffle(token1, token2, contractions, finalToken):
        shuffle=[]

        for c in contractions:
            token1 = token1[:c[0]] + token1[c[0]+1:]
            token2 = token2[:c[1]] + token2[c[1]+1:]

        token=token1+token2
        for i, c in enumerate(finalToken):
            j=token.find(c)
            shuffle.append(j)
        return shuffle

    shuffle = ParseShuffle(tokens[0], tokens[1], contractions, tokens[2])

    def ArrToStr(arr):
        return str(arr).replace('[', '{').replace(']', '}').replace('\'', '')

    outputStr = f"//{einsumStr}\n"

    if lineType == LINE_TMP:
        outputStr += (
                    f"tmp = {tensor1}.Slice{len(slice1[0])}({ArrToStr(slice1[0])}, {ArrToStr(slice1[1])})\n"
                    f"     .Contract{len(contractions)}(\n"
                    f"          {tensor2}.Slice{len(slice2[0])}({ArrToStr(slice2[0])}, {ArrToStr(slice2[1])}),\n"
                    f"     ({ArrToStr(contractions)}))\n"
                    f" .Shuffle{len(shuffle)}({ArrToStr(shuffle)});\n")
    elif lineType == LINE_B:
        outputStr += (
                    f"B -= {tensor2}.Slice{len(slice2[0])}({ArrToStr(slice2[0])}, {ArrToStr(slice2[1])})\n"
                    f"     .Contract{len(contractions)}({tensor1},\n"
                    f"          ({ArrToStr(contractions)}))\n"
                    f" .Shuffle{len(shuffle)}({ArrToStr(shuffle)});\n")
    elif lineType == LINE_UNKNOWN:
        outputStr += (
                    f"{tensor1}.Slice{len(slice1[0])}({ArrToStr(slice1[0])}, {ArrToStr(slice1[1])})\n"
                    f"     .Contract{len(contractions)}(\n"
                    f"          {tensor2}.Slice{len(slice2[0])}({ArrToStr(slice2[0])}, {ArrToStr(slice2[1])}),\n"
                    f"     ({ArrToStr(contractions)}))\n"
                    f" .Shuffle{len(shuffle)}({ArrToStr(shuffle)})\n")
        
        

    return outputStr


if(sys.argv[1] == "file"):
    with open(sys.argv[2], 'r') as file:
        for line in file:
            line = line.strip()

            lineType = LINE_UNKNOWN
            if line.find('tmp =') != -1:
                lineType = LINE_TMP
            elif line.find('B -=') != -1:
                lineType = LINE_B

            startIdx = line.find('np.einsum(')
            einsumArgs = line[startIdx:]

            startIdx = einsumArgs.find('\'')+1
            endIdx = einsumArgs.find('\'', startIdx)
            einsumStr = einsumArgs[startIdx:endIdx]
            
            tensor1 = ''
            startIdx = einsumArgs.find(',', endIdx)+1
            endIdx = einsumArgs.find('[', startIdx)
            if endIdx == -1:
                endIdx = line.find(',', startIdx)
                tensor1 = einsumArgs[startIdx:endIdx]
                startIdx = endIdx+1
            else:
                tensor1 = einsumArgs[startIdx:endIdx]
                startIdx = einsumArgs.find('],', endIdx)+2

            endIdx = einsumArgs.find('(', startIdx)

            tensor2 = einsumArgs[startIdx:endIdx]

            tensor1= tensor1.strip()
            tensor2 = tensor2.replace('mp2f12.', '')
            tensor2= tensor2.strip()
            
            if(tensor1 == 'k'):
                tensor1='K'
            
            print(ParseEinsum(einsumStr, tensor1, tensor2, lineType))
            #print(f"T1={tensor1}, T2={tensor2}, Einsum=\"{einsumStr}\", Type={lineType}")
else:
    tensor1 = sys.argv[1]
    tensor2 = sys.argv[2]
    einsumStr = sys.argv[3]
    print(ParseEinsum(einsumStr, tensor1, tensor2, LINE_UNKNOWN))


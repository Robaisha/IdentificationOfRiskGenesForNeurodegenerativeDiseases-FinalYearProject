import sys
import numpy as np
import pandas as pd

matrixFile = sys.argv[1]
familyFile = sys.argv[2]
outputFile = sys.argv[3]
traitFile = sys.argv[4]
print(sys.argv)

def readMatrixFile(matrixFile, matrixBegin, sampleID, trait):
    idExprDict = {}
    with open(matrixFile, "r", encoding="utf-8") as f_in:
        lines = f_in.readlines()
        for line in lines[matrixBegin:-1]:
            l = line.strip("\n").split("\t")
            ID = l[0].strip("\"")
            exprList = l[1:]
            if ID not in idExprDict:
                idExprDict[ID] = []
                idExprDict[ID].append(exprList)
            else:
                print("ID matched with multiExprValue")
                print(l)
    headline = lines[sampleID]
    traitInfo = lines[trait]

    return [idExprDict, headline, traitInfo]

def readFamilyFile(familyFile, tableBegin, tableEnd, geneSymbolCol):
    idGeneDict = {}
    with open(familyFile, "r", encoding="utf-8") as f_in:
        lines = f_in.readlines()
    with open(familyFile, "r", encoding="utf-8") as f_in:
        lines = f_in.readlines()
        for line in lines[tableBegin:tableEnd]:
            l = line.strip("\n").split("\t")
            geneSymbol = l[geneSymbolCol] if len(l) > geneSymbolCol else None
            ID = l[0]
            if ID.startswith("ILMN_"):
                ID = ID.replace("ILMN_", "")
            geneName = l[geneSymbolCol] if len(l) > geneSymbolCol else None
            if ID not in idGeneDict:
                if geneName != "":
                    idGeneDict[ID] = geneName

    print("idGeneDict done!")
    return idGeneDict

def removeDulProbes(idExprDict, idGeneDict):
    geneExprDict = {}

    for symbolID in idExprDict:
        symbolID_without_ILMN = symbolID.replace("ILMN_", "")
        # print(symbolID_without_ILMN, 'symbolid in array 1')
        if symbolID_without_ILMN in idGeneDict:
            # print(symbolID_without_ILMN, 'symbolid in array 2')
            geneName = idGeneDict[symbolID_without_ILMN]
            if geneName not in geneExprDict:
                geneExprDict[geneName] = []
            geneExprDict[geneName].extend(idExprDict[symbolID])

    for geneID in geneExprDict:
        i = len(geneExprDict[geneID])
        exp = np.array(geneExprDict[geneID])
        exp = exp.astype(float)
        a = np.sum(exp, axis=0) / i
        geneExprDict[geneID] = a.tolist()

    return geneExprDict


def outputMatrixFile(outputFile, headline, geneExprDict):
    with open(outputFile, "w") as f_out:
        f_out.write(headline)
        for key in geneExprDict:
            content = key + "\t" + "\t".join(str(i) for i in geneExprDict[key]) + "\n"
            f_out.write(content)

def outputTraitFile(traitFile, headline, traitInfo):
    with open(traitFile, "w") as f_out:
        f_out.write(headline)
        f_out.write(traitInfo)

# 01 read the matrix file
print("Loading matrix file")
try:
    l = readMatrixFile(matrixFile, 37, 36, 8)
    idExprDict = l[0]
    headline = l[1]
    traitInfo = l[2]
    print("idExprDict done!")
except Exception as e:
    print("An error occurred while reading the matrix file:", e)

# 02 read the family file
print("Loading family file")
try:
    idGeneDict = readFamilyFile(familyFile, 511, 48618, 0)
except Exception as e:
    print("An error occurred while reading the family file:", e)

# 03 remove the duplicate gene symbols
print("Removing the duplicate gene symbols......")
try:
    geneExprDict = removeDulProbes(idExprDict, idGeneDict)
    print("Gene symbols are unique")
    print("Length of idExprDict:", len(idExprDict))
    print("Length of idGeneDict:", len(idGeneDict))
    print("Length of geneExprDict:", len(geneExprDict))
except Exception as e:
    print("An error occurred while removing duplicate gene symbols:", e)

# 04 output matrixFile and traitFile
try:
    outputMatrixFile(outputFile, headline, geneExprDict)
    outputTraitFile(traitFile, headline, traitInfo)
except Exception as e:
    print("An error occurred while outputting matrixFile and traitFile:", e)

# Read the text file
data_Matrix = pd.read_table("C:\\Users\\Masoo\\Desktop\\FYP_Material\\ProbeToGene\\output_matrix.txt", delimiter="\t")  # Assuming tab-separated, adjust delimiter as needed
# Write the data to a CSV file
data_Matrix.to_csv("C:\\Users\\Masoo\\Desktop\\FYP_Material\\ProbeToGene\\output_matrix.csv", index=False)  # Specify index=False if you don't want row numbers in the CSV file
# Read the text file
data_Matrix = pd.read_table("C:\\Users\\Masoo\\Desktop\\FYP_Material\\ProbeToGene\\output_trait.txt", delimiter="\t")  # Assuming tab-separated, adjust delimiter as needed
# Write the data to a CSV file
data_Matrix.to_csv("C:\\Users\\Masoo\\Desktop\\FYP_Material\\ProbeToGene\\output_trait.csv", index=False)  # Specify index=False if you don't want row numbers in the CSV file



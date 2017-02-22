__author__ = 'surya'

from GeneralMethods import Generalmethods


def tagPeptideCount(pairDic,prot2peptideFile,outName,header):
    peptideDic={}
    prot2peptideDic=Generalmethods.createDic(prot2peptideFile,0,[1,2],header=True)
    header=header+"\tCommonCount\tpep1Count\tpep2Count"
    output=open(outName,"w")
    output.write(header+"\n")
    for eachpair in pairDic:
        p1=pairDic[eachpair][0]
        p2=pairDic[eachpair][1]
        if p1 in prot2peptideDic and p2 in prot2peptideDic:
            pep1=prot2peptideDic[p1][0].split(";")
            pep2=prot2peptideDic[p2][0].split(";")
            commonPeptide=list(set(pep1).intersection(pep2))
            wline=pairDic[eachpair]+[str(len(commonPeptide)),str(len(pep1)),str(len(pep2))]
        else:
            print " one of the protein is not found ",p1, p2
        output.write("\t".join(wline)+"\n")
        peptideDic[eachpair]=wline
    output.close()
    return peptideDic,header

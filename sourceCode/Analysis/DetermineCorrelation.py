__author__ = 'surya'

## find of the common peptide between both proteins in the exp is same

def OnlyReturnNotCommonPepCount(assay1D,assay2D,commonD,minComAssay,minpep):
    ass1={}
    ass2={}
    commonGood=False
    for assName1 in assay1D:
        if assName1 in commonD:
            commonpep=list(set(assay2D[assName1]).intersection(assay1D[assName1]))
            remain1=len(assay1D[assName1])-len(commonpep)
            remain2=len(assay2D[assName1])-len(commonpep)
            if remain1>=minpep:
                ass1[assName1] = remain1
            if remain2>=minpep:
                ass2[assName1] = len(assay2D[assName1])-len(commonpep)
    newcommon=list(set(ass2.keys()).intersection(ass1.keys()))
    if len(newcommon)>minComAssay:
        commonGood=True
        ass1.update({a1:len(assay1D[a1]) for a1 in assay1D if a1 not in commonD})
        ass2.update({a2:len(assay2D[a2]) for a2 in assay2D if a2 not in commonD})

    return commonGood,ass1,ass2,newcommon






## calculate the p value, k and n for the bionomial

## take two proteins, determine the assay they are involve in, if assay is not present in other protein than assign it zero

def createArrays(AssayDic1,AssayDic2,CommAssayList):#,p1,p2): ## val is assay:pepcount
    allAssay1=[]
    allAssay2=[]
    AllAssayList=list(set(AssayDic1.keys()+AssayDic2.keys()))
    for eachAss in AllAssayList:
        if eachAss in AssayDic1:
            allAssay1+=[AssayDic1[eachAss]]
        else:
            allAssay1+=[0]
        if eachAss in AssayDic2:
            allAssay2+=[AssayDic2[eachAss]]
        else:
            allAssay2+=[0]
    return allAssay1,allAssay2#,total1,total2,commonSum1,commonSum2


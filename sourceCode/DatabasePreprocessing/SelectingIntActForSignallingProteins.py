__author__ = 'surya'

def selectPairsIntAct(intActFile,proteinList):
    pairDic={}
    with open(intActFile) as intact:
        next(intact)
        for line in intact:
            splits=line.split("\t")
            if splits[0].strip() in proteinList and splits[1].strip() in proteinList:
                pair1=splits[0].strip()+"_"+splits[1].strip()
                pair2=splits[1].strip()+"_"+splits[0].strip()
                if pair1 not in pairDic and pair2 not in pairDic:
                    pairDic[pair1]=splits[5].strip() ## IntActScore
    print "with total ", len(proteinList), " protein total pair found in intact are ", len(pairDic)
    return pairDic
__author__ = 'surya'

def createProtein2Path(proteinList,file):
    prot2path={}
    for line in open(file):
        splits=line.strip().split("\t")
        reactomeId=splits[0].strip()
        proteins=splits[1:len(splits)]
        for each in proteins:
            if each in proteinList:
                if each not in prot2path:
                    prot2path[each]=[reactomeId]
                else:
                    prot2path[each].append(reactomeId)
    return prot2path

## make a dic for protein2pathway
def dicFromProtein2Pathway(prot2pathwayFile):#,leafPathways,leafCheck=True):
    ReactomeDic={}
    prot2path={}
    for line in open(prot2pathwayFile):
        splits=line.strip().split("\t")
        proteinId=splits[0].strip().split("-")[0] ## as most of the protein have uniprot id "-" number
        if splits[5].strip()=="Homo sapiens":
            reactomeId=splits[1].strip()
            if proteinId not in prot2path:
                prot2path[proteinId]=[reactomeId]
            else:
                prot2path[proteinId].append(reactomeId)

            if reactomeId not in ReactomeDic:
                ReactomeDic[reactomeId]=[proteinId]
            else:
                ReactomeDic[reactomeId].append(proteinId)

    return prot2path,ReactomeDic

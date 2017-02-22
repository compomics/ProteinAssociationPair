__author__ = 'surya'

from GeneralMethods import GeneralStat


def createAssayDic(line,check=False,checkList=[],MinPep=1):
    dic={}
    for each in line.strip("|").split("|"):
        eachList=each.split(":")
        if check:
            if eachList[0] in checkList:
                if int(eachList[1])>=MinPep:
                    dic[eachList[0]]=int(eachList[1])
        else:
            if int(eachList[1])>=MinPep:
                dic[eachList[0]]=int(eachList[1])
        # if len(dic)>50
    return dic

def createAssayDicWithPepName(line,MinPep=1,AssayCheckList=[],AssCheck=False):
    dic={}
    for each in line.strip("|").split("|"):
        eachList=each.split(":")
        assay=eachList[0]
        if AssCheck:
            if assay in AssayCheckList:
                peptideList=list(set(eachList[1].strip(";").split(";")))
                if len(peptideList)>=MinPep:
                    dic[assay]=peptideList
        else:
            peptideList = list(set(eachList[1].strip(";").split(";")))
            if len(peptideList) >= MinPep:
                dic[assay] = peptideList
                # if len(dic)>50
    return dic

def compare2Dic(dic1,dic2):
    assList=[]
    com=0
    assay1=dic1.keys()
    assay2=dic2.keys()
    totalAss=assay1+assay2
    for each in totalAss:
        if each not in assList:
            assList.append(each)
            if each in assay1 and each in assay2:
                com+=1
    return com,assList

def createAssayDicWithPresentOrNot(line,check=False,checkList=[],MinPep=0):
    dic={}
    for each in line.strip("|").split("|"):
        eachList=each.split(":")
        if check:
            if eachList[0] in checkList:
                if int(eachList[1])>=MinPep:
                    dic[eachList[0]]=1
        else:
            if int(eachList[1])>=MinPep:
                dic[eachList[0]]=1
        # if len(dic)>50
    return dic

def createAssayDicWithNormPSM(line,check=False,checkList=[],MinPep=0,protLength=1):
    dic={}
    for each in line.strip("|").split("|"):
        eachList=each.split(":")
        if check:
            if eachList[0] in checkList:
                if int(eachList[1])>=MinPep:
                    dic[eachList[0]]=float(eachList[1])/protLength
        else:
            if int(eachList[1])>=MinPep:
                dic[eachList[0]]=float(eachList[1])/protLength
        # if len(dic)>50
    return dic


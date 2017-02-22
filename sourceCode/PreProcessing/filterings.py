__author__ = 'surya'


from Analysis import GetAssayCount

def preFiltrationWithHighlyOCcuringProteins(MSFile,uniprotProteinList,out):
    assCountList=[]
    with open(MSFile) as OpenedMSFile:
        next(OpenedMSFile)
        for line in OpenedMSFile:
            colList=line.split("\t")
            assCount = int(colList[1].strip())
            if colList[0].strip() in uniprotProteinList:
                out.write(line)
                assCountList.append(assCount)
    assCountList.sort()
    maxCutoff=assCountList[-10]
    print "the maximum cutoff is ",maxCutoff
    out.close()
    return maxCutoff



def CheckAssayWithPeptideIds(MSFile,AssayLen,assayCheck=False,PhosList=[],proteincheck=False,proteinList=[],
                             MinPep=0,SelectedAssayList=[],selectedProteinFile=""):#,crap=[]):
    maxAssayLen=preFiltrationWithHighlyOCcuringProteins(MSFile,proteinList,selectedProteinFile)
    proteinSkiped=0
    protein2Assay={}
    protein2AssayPepString={}
    count=0
    assayList=[]
    check=False
    totalProt=0
    with open(MSFile) as OpenedMSFile:
        next(OpenedMSFile)
        for line in OpenedMSFile:
            colList=line.split("\t")
            totalProt+=1
            assCount = int(colList[1].strip())
            if assCount>=maxAssayLen:
                print colList[0].strip(),assCount
            if assCount > AssayLen and assCount<maxAssayLen:
                ##first select the threshold
                if proteincheck:
                    if colList[0].strip() in proteinList:# and colList[0].strip() not in crap:
                        check=True
                    else:
                        proteinSkiped+=1
                        check=False
                else:
                    check=True
                if check:
                    AssayDic=GetAssayCount.createAssayDicWithPepName(colList[2].strip(),MinPep=MinPep,
                                                                     AssayCheckList=SelectedAssayList,AssCheck=assayCheck)
                    if len(AssayDic)<=AssayLen: #and AssayDic[AssayDic.keys()[0]]==1:
                        count+=1
                        # print(AssayDic)
                    elif len(AssayDic)>AssayLen:
                        # if len(protein2Assay)==50:
                        #     break
                        protein2Assay[colList[0].strip()]=AssayDic
                        protein2AssayPepString[colList[0].strip()] = colList[2].strip()
                        # print colList[2].strip()
                        for each in AssayDic:
                            if each not in assayList:
                                assayList.append(each)
    # for pro in protein2Assay:
    #     selectedProteinFile.write(pro+"\n")
    # selectedProteinFile.close()
    ln="total removed assays are "+str(count)+"\n"+"total assay found are "+str(len(assayList))+"\n"+\
       "total proteins remaining are "+str(len(protein2Assay))+" out of total proteins "+str(totalProt)+"\n"+\
       "total protein skipped because they are not present in pathway/uniprot are "+str(proteinSkiped)+"\n"
    print ln
    return protein2Assay,assayList,protein2AssayPepString,ln





from GeneralMethods import GeneralStat
from Analysis import DetermineCorrelation
import pdb
## create the pairs from the file
# select the pair which have atleast one com assay with more than 2 peptide
## determine the correaltion btw the proteins


def createPairsWithPepCheck(FilProt2Ass2PepCountDic,minAssayCount,IntActPairDic,jaccCutoff,outWrite1,
                        minPepCount,protein2AssayPepString,outWrite5):

    outWrite1.write("protein1\tprotein2\tjaccSim\tCommonAssayCount\tAssay1\tAssay2\tIntActFound\n")

    pairDic={}
    c,oc,tp=0,0,0
    protAfterselection={}
    uniqueProteinList=FilProt2Ass2PepCountDic.keys()
    for index1 in range(len(uniqueProteinList)-1):
        protein1=uniqueProteinList[index1]
        for index2 in range(index1+1,len(uniqueProteinList)):
            protein2=uniqueProteinList[index2]
            pair1=protein1+"_"+protein2
            pair2=protein2+"_"+protein1
            if pair1 not in pairDic and pair2 not in pairDic and protein2AssayPepString[protein1]!=protein2AssayPepString[protein2]:
                ## get the assay dic for each proteins
                assayDic1=FilProt2Ass2PepCountDic[protein1]
                assayDic2=FilProt2Ass2PepCountDic[protein2]
                     ## compare the assays to see common entries between them
                CommonAssayList=list(set(assayDic1.keys()).intersection(assayDic2.keys()))

                if len(CommonAssayList) >minAssayCount:
                    CommonGood,NewAss2PepCount1,NewAss2PepCount2,newassCommon=DetermineCorrelation.OnlyReturnNotCommonPepCount(assayDic1,assayDic2,CommonAssayList,minAssayCount,minPepCount)
                    # CommonAssayCount,AllAssayList=GetAssayCount.compare2Dic(assayDic1,assayDic2)
                ## make a list of the values in the assays if there is no entry than add 0
                    if CommonGood==True:
                        allAssay1,allAssay2=DetermineCorrelation.createArrays(NewAss2PepCount1,NewAss2PepCount2,CommonAssayList)
                        JaccardSim=GeneralStat.jaccardSimilarity(allAssay1,allAssay2) ## this takes two list of two arrays and check the pearson correlations
                        if JaccardSim>=0.1:
                            if pair1 in IntActPairDic or pair2 in IntActPairDic:
                                intactC=1
                            else:
                                intactC=0
                            outWrite1.write(protein1+"\t"+protein2+"\t"+str(JaccardSim)+"\t"+str(len(newassCommon))+"\t"+str(len(NewAss2PepCount1))+"\t"+str(len(NewAss2PepCount2))+"\t"+str(intactC)+"\n")
                            outWrite5.write(protein1+"\t"+protein2+"\t"+str(JaccardSim)+"\t".join(newassCommon)+"\n")
                            tp+=1
                            if JaccardSim>=jaccCutoff:
                                pairDic[pair1]=[JaccardSim,intactC,len(NewAss2PepCount1),len(NewAss2PepCount2),len(newassCommon)] ## pearson correlation,intactCount (if presnt in intact
                                if protein1 not in protAfterselection:
                                    protAfterselection[protein1]=0
                                if protein2 not in protAfterselection:
                                    protAfterselection[protein2]=0
                                c+=1
    ln="total pair found are "+str(len(pairDic))+ " from total pairs "+str(tp)+" from total unique proteins are "+str(len(protAfterselection))+"\n"
    print ln
    outWrite1.close()
    outWrite5.close()
    return pairDic,protAfterselection,ln

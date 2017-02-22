__author__ = 'surya'

from DatabasePreprocessing import protein2pathway
## check for each pathway how many proteins retrieve back

## check for each pathway how many proteins retrieve back

def checkTheCoverageForLeafPathway(prot2pathwayFile,AllProteinPairsDic,outWrite4,outWrite2,outWrite3,PathId2NameDic,prot2PathDicAvailable=False,ReactDic={}):
    header="protein1\tprotein2\tscore\tassay1\tassay2\tCommAssay\tintAct\tpath1\tpath2\tpathC1\tpathC2\tComPath\tComPathCount"
    outWrite2.write(header+"\n")
    outWrite3.write(header+"\n")

    pairFound={}
    if prot2PathDicAvailable:
        prot2path=prot2pathwayFile
        ReactomeDic=ReactDic
    else:
        prot2path,ReactomeDic=protein2pathway.dicFromProtein2Pathway(prot2pathwayFile)#,leafPathways)
    path2PPI={}
    path2uniqueProtein={}
    pair2pathTag={}
    for pairP in AllProteinPairsDic:
        p1=pairP.split("_")[0]
        p2=pairP.split("_")[1]
        if p1 in prot2path and p2 in prot2path:
            for pathP1 in prot2path[p1]:
                if pathP1 in prot2path[p2]:
                    if pairP not in pairFound:
                        pairFound[pairP]=[pathP1]
                    else:
                        pairFound[pairP].append(pathP1)
                    if pathP1 in path2PPI:
                        if p1 not in path2uniqueProtein[pathP1]:
                            path2uniqueProtein[pathP1].append(p1)
                        if p2 not in path2uniqueProtein[pathP1]:
                            path2uniqueProtein[pathP1].append(p2)
                        if pairP not in path2PPI[pathP1][0]:
                            path2PPI[pathP1][0].append(pairP)
                            path2PPI[pathP1][1]+=1
                            path2PPI[pathP1][2].append(str(AllProteinPairsDic[pairP][0])) # score
                            path2PPI[pathP1][3]+=AllProteinPairsDic[pairP][1] # intact
                    else:
                        path2uniqueProtein[pathP1]=[p1,p2]
                        path2PPI[pathP1]=[[pairP],1,[str(AllProteinPairsDic[pairP][0])],AllProteinPairsDic[pairP][1]]
    print "total leaf pathway remaining are ",len(path2PPI), " out of "
    print "total pairs use to map to pathways are ", len(pairFound)
    for eachPair in AllProteinPairsDic:
        prot1=eachPair.split("_")[0].strip()
        prot2=eachPair.split("_")[1].strip()
        if eachPair in pairFound:
            pair2pathTag[eachPair]=[prot1,prot2,str(AllProteinPairsDic[eachPair][0])
                            ,str(AllProteinPairsDic[eachPair][2]),str(AllProteinPairsDic[eachPair][3])
                            ,str(AllProteinPairsDic[eachPair][4]),str(AllProteinPairsDic[eachPair][1]),
                            ";".join(prot2path[prot1]),";".join(prot2path[prot2]),
                            str(len(prot2path[prot1])),str(len(prot2path[prot2])),
                            ";".join(pairFound[eachPair]),str(len(pairFound[eachPair]))]
            outWrite2.write("\t".join(pair2pathTag[eachPair])+"\n")
        elif eachPair not in pairFound:
            if prot1 not in prot2path:
                ppath1="NF"
                p1pathcount=0
            else:
                ppath1=";".join(prot2path[prot1])
                p1pathcount=len(prot2path[prot1])
            if prot2 not in prot2path:
                ppath2="NF"
                p2pathcount=0
            else:
                ppath2=";".join(prot2path[prot2])
                p2pathcount=len(prot2path[prot2])

            pair2pathTag[eachPair] = [eachPair.split("_")[0].strip(), eachPair.split("_")[1].strip(),
                                      str(AllProteinPairsDic[eachPair][0]),str(AllProteinPairsDic[eachPair][2]),
                                      str(AllProteinPairsDic[eachPair][3]),str(AllProteinPairsDic[eachPair][4]),
                                      str(AllProteinPairsDic[eachPair][1]),ppath1,ppath2,str(p1pathcount),
                                      str(p2pathcount),"NF\t0"]
            outWrite3.write("\t".join(pair2pathTag[eachPair])+"\n")

    outWrite2.close()
    outWrite3.close()
    outWrite4.write("pathway\tpathwayName\tproteins\tTotalProtCount\tPPI\tPPICount\tAllPPIScores\tProtFound\tProtFoundCount\tpercentage\tInIntActFound\tminScore\tMaxScore\n")
    for eachpath in path2PPI:
        limit=len(path2PPI[eachpath][0])
        # if len(path2PPI[eachpath][0])>100:
        #     limit=100
        # else:
        #     limit=len(path2PPI[eachpath][0])
        #print path2PPI[eachpath][2][0:limit]
        AllScore=path2PPI[eachpath][2]
        outWrite4.write(eachpath+"\t"+PathId2NameDic[eachpath]+"\t"+"_".join(ReactomeDic[eachpath])+"\t"+str(len(ReactomeDic[eachpath]))
                        +"\t"+";".join(path2PPI[eachpath][0][0:limit])+"\t"+str(path2PPI[eachpath][1])+"\t"+";".join(AllScore[0:limit])
                        +"\t"+";".join(path2uniqueProtein[eachpath])+"\t"+str(len(path2uniqueProtein[eachpath]))
                        +"\t"+str(round((float(len(path2uniqueProtein[eachpath]))/len(ReactomeDic[eachpath]))*100,2))+"\t"+str(path2PPI[eachpath][3])
                        +"\t"+str(min(AllScore))+"\t"+str(max(AllScore))+"\n")

    outWrite4.close()
    return pair2pathTag,header



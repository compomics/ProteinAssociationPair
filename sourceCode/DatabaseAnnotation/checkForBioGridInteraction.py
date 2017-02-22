__author__ = 'surya'



def ParseBioGrid(BiogridPairFile):
    genePairDic={}
    with open(BiogridPairFile) as BPFile:
        next(BPFile)
        for line in BPFile:
            splits=line.strip().split("\t")
            p1=splits[1]
            p2=splits[2]
            pair1=p1+"_"+p2
            pair2=p2+"_"+p1
            if pair1 not in genePairDic and pair2 not in genePairDic:
                genePairDic[pair1]=line.strip()
    print "total geneID interacting pairs found are ",len(genePairDic)," for human only.."
    # outfile.close()
    return genePairDic


def checkForProteinList(pairdic,outputfilename,protein2GeneID,BiogridPairsFile,header):
    protDic = {}
    geneIDPairs = ParseBioGrid(BiogridPairsFile)
    totalFound=0
    new=0
    outfile=open(outputfilename,"w")
    headerList=header.split("\t")
    newheader="\t".join(headerList[0:9] + ["geneIds1", "geneIds2","ComGene","ComGeneCount","BiogridInt"] +
                            headerList[9:])
    outfile.write( newheader+ "\n")
    for eachPair in pairdic:
        p1 = pairdic[eachPair][0]
        p2 = pairdic[eachPair][1]
        found=0
        if p1 in protein2GeneID and protein2GeneID[p1][3].strip()!="":
            geneidList1=protein2GeneID[p1][3].strip(";").split(";")
            geneIds1="_".join(geneidList1)
        else:
            geneidList1=[]
            geneIds1="NF"
        if p2 in protein2GeneID and protein2GeneID[p2][3].strip()!="":
            geneidList2 = protein2GeneID[p2][3].strip(";").split(";")
            geneIds2 = "_".join(geneidList2)
        else:
            geneidList2=[]
            geneIds2="NF"
        if len(geneidList1)>0 and len(geneidList2)>0:
            commonGene=list(set(geneidList1).intersection(geneidList2))
            comGeneCount=str(len(commonGene))
            if len(commonGene)==0:
                comGenName="NF"
            else:
                comGenName="_".join(commonGene)
        else:
            comGenName="NA"
            comGeneCount="0"
        if p1 in protein2GeneID and p2 in protein2GeneID:
            for eachp1 in protein2GeneID[p1][3].strip(";").split(";"):
                for eachp2 in protein2GeneID[p2][3].strip(";").split(";"):
                    genePair1=eachp1+"_"+eachp2
                    genePair2=eachp2+"_"+eachp1
                    if genePair1 in geneIDPairs or genePair2 in geneIDPairs:
                        found+=1
            if found>=1:
                interact="yes"
                totalFound += 1
                if int(pairdic[eachPair][8])==0:
                    new+=1
            else:
                interact="no"
        else:
            interact="NA"
        wline=pairdic[eachPair][0:9]+[geneIds1,geneIds2,comGenName,comGeneCount,interact]+pairdic[eachPair][9:]
        protDic[eachPair] = wline
        outfile.write("\t".join(wline)+"\n")
    outfile.close()
    print "total pairs found to be interacting ",totalFound," but new pairs found interacting are ",new
    return protDic,newheader

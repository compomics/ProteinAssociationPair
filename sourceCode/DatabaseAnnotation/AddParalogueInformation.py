__author__ = 'surya'

""" add the paralog information about the two genes both proteins comes from and than check if they come from paralog proteins or same genes """


from GeneralMethods import Generalmethods

def findParalogue(paralogFile,Uniprot2EnsemblFile,pairDic,outname,header,geneParalogDic=False):#,old2newUniName):
    paralogDic={}
    out = open(outname, "w")
    header=header+"\tGene1\tGene2\tComGenesCount\tComGenes\tParalog"
    out.write(header+"\n")

    ## get all the genes ids for the specific proteins
    protein2genes=Generalmethods.createDic(Uniprot2EnsemblFile,0,1,header=True,multiple=False,addDublicate=True)
    print "total prot2gene entries are " ,len(protein2genes)
    ## get the gene paralogue
    if not geneParalogDic:
        geneParalogDic=Generalmethods.createDic(paralogFile,0,1,header=True,multiple=False,addDublicate=True)
        print "total number of gene pralog are: ",len(geneParalogDic)
    else:
        geneParalogDic=paralogFile
    # now look for the genes either which are common to proteins or have paralogs
    paracount=0
    for eachPair in pairDic:
        newp1=pairDic[eachPair][0]
        newp2=pairDic[eachPair][1]
        if newp1 in protein2genes and newp2 in protein2genes:
            # if protein2genes[p1]!="NF" and protein2genes[p2]!="NF":
            gene2=protein2genes[newp2]
            gene1=protein2genes[newp1]
            geneNAme1=";".join(protein2genes[newp1])
            geneNAme2=";".join(protein2genes[newp2])
            comGene=list(set(gene1).intersection(gene2))
            paralogList=[]
            for eachGene2 in gene2:
                if eachGene2 in geneParalogDic:
                    paralogList+=geneParalogDic[eachGene2]
            paralogList=list(set(paralogList))
            paralogPresent=list(set(gene1).intersection(paralogList))
            if len(paralogPresent)!=0:
                paralog="yes"
                paracount+=1
            else:
                paralog="No"
            if len(comGene)!=0:
                common=";".join(protein2genes[newp2])
                comNum=str(len(comGene))
            else:
                common="NF"
                comNum="0"
        else:
            if newp1 in protein2genes:
                geneNAme1=";".join(protein2genes[newp1])
                geneNAme2="NF"
            elif newp2 in protein2genes:
                geneNAme2=";".join(protein2genes[newp2])
                geneNAme1="NF"
            else:
                geneNAme1="NF"
                geneNAme2="NF"
            paralog="NA"
            common="NA"
            comNum="NA"
        wline=pairDic[eachPair]+[geneNAme1,geneNAme2,comNum,common,paralog]
        paralogDic[eachPair]=pairDic[eachPair]+[geneNAme1,geneNAme2,comNum,common,paralog]
        out.write("\t".join(wline)+"\n")
    out.close()
    print "total paralog found are ",paracount
    return paralogDic,header



def paralogueAnnotationForRandomSample(geneParalogDic,protein2genes,pairDic):
    paracount = 0
    for eachPair in pairDic:
        pairSplit=eachPair.split("_")
        newp1 = pairSplit[0]
        newp2 = pairSplit[1]
        if newp1 in protein2genes and newp2 in protein2genes:
            # if protein2genes[p1]!="NF" and protein2genes[p2]!="NF":
            gene2 = protein2genes[newp2]
            gene1 = protein2genes[newp1]
            geneNAme1 = ";".join(protein2genes[newp1])
            geneNAme2 = ";".join(protein2genes[newp2])
            comGene = list(set(gene1).intersection(gene2))
            paralogList = []
            for eachGene2 in gene2:
                if eachGene2 in geneParalogDic:
                    paralogList += geneParalogDic[eachGene2]
            paralogList = list(set(paralogList))
            paralogPresent = list(set(gene1).intersection(paralogList))
            if len(paralogPresent) != 0:
                paralog = "yes"
                paracount += 1
    return paracount

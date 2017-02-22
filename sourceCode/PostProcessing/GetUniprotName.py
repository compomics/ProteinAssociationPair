__author__ = 'surya'


## get all the pairs and parse
def AddProteinName(pairdic,UniprotDic,outputfilename,header):
    protDic={}
    outfile = open(outputfilename, "w")
    newheader = "protein1\tprotein2\tpName1\tpName2\t" + "\t".join(header.split("\t")[2:])
    outfile.write(newheader+"\n")
    found=0
    for eachPair in pairdic:
        p1=pairdic[eachPair][0]
        p2=pairdic[eachPair][1]
        if p1 in UniprotDic and p2 in UniprotDic:
            wline = [pairdic[eachPair][0],pairdic[eachPair][1],UniprotDic[p1][1],UniprotDic[p2][1]]+pairdic[eachPair][2:]
            protDic[eachPair]=wline
            outfile.write("\t".join(wline) + "\n")
    outfile.close()
    print 'total pairs reamining after naming are ',len(protDic)
    return protDic,newheader


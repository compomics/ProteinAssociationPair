__author__ = 'surya'

def labelPairWithFilteredProteins(pairDic,outputfilename,header):
    count=0
    newPairDic={}
    outfile=open(outputfilename,"w")
    newheader = header + "\tkeratin\tIg"
    outfile.write(newheader + "\n")
    for eachPair in pairDic:
        p1Name=pairDic[eachPair][2]
        p2Name=pairDic[eachPair][3]
        keraPresent,IgPresent=0,0
        if "Ig" in p1Name or "Immunoglobulin" in p1Name:
            IgPresent += 1
        elif "Keratin," in p1Name:
            keraPresent += 1
        if "Ig" in p2Name or "Immunoglobulin" in p2Name:
            IgPresent += 1
        elif "Keratin," in p2Name:
            keraPresent += 1
        ln=pairDic[eachPair]+[str(keraPresent),str(IgPresent)]
        newPairDic[eachPair]=ln
        outfile.write("\t".join(ln)+"\n")
    print "total pair skipped are: ",count
    outfile.close()
    return newPairDic,newheader
__author__ = 'surya'


def tagCraps(crapdic,pairdic,outputName,header):
    crapPairDic={}
    output=open(outputName,"w")
    header=header+"\tCrap"
    output.write(header+"\n")
    found=0
    for eachPair in pairdic:
        p1=pairdic[eachPair][0]
        p2=pairdic[eachPair][1]
        if p1 in crapdic or p2 in crapdic:
            crap="Yes"
            found+=1
        else:
            crap="No"
        wline=pairdic[eachPair]+[crap]
        crapPairDic[eachPair]=wline
        output.write("\t".join(wline)+"\n")
    print "total crap proteins found as pair are ",found
    output.close()
    return crapPairDic,header



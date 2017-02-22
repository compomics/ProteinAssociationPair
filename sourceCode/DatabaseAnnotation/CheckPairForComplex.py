__author__ = 'surya'


def ParseCorumComplexFile(file1,CheckcolNumValue=[0,[0]]):
    Complex2Protdic={}
    prot2complex={}
    with open(file1) as filename:
        next(filename)
        for line in filename:
            splits=line.split("\t")
            complex=splits[0].strip()+"_"+splits[1].strip() #id+name
            proteinsList=splits[5].strip().split(";")
            if splits[CheckcolNumValue[0]].strip() in CheckcolNumValue[1]:
                if complex not in Complex2Protdic:
                    Complex2Protdic[complex]=splits[4].strip()
                else:
                    print " more than one complex are present for ", complex
                for prot in proteinsList:
                    if prot not in prot2complex:
                        prot2complex[prot]=[complex]
                    else:
                        prot2complex[prot].append(complex)

    return Complex2Protdic,prot2complex


def tagComplexPair(complexFile,pairDic,outname,header):
    complexDic={}
    output=open(outname,"w")
    header=header+"\tComplexFound\tProtComplex\tProtCompCount"
    output.write(header+"\n")

    complex2protDic,prot2ComplexDic=ParseCorumComplexFile(complexFile,CheckcolNumValue=[2,["Human"]])

    print "total complexes found for human are ",len(complex2protDic)
    print " total proteins present are ",len(prot2ComplexDic)

    found,total=0,0
    for eachPair in pairDic:
        p1=pairDic[eachPair][0]
        p2=pairDic[eachPair][1]
        total+=1
        if p1 in prot2ComplexDic and p2 in prot2ComplexDic:
            commonComplex=list(set(prot2ComplexDic[p1]).intersection(prot2ComplexDic[p2]))
            if len(commonComplex)!=0:
                found+=1
                AllComplex=";".join(commonComplex)
                allComplexProt=""
                allCompProtCount=""
                for comp in commonComplex:
                    allComplexProt+=complex2protDic[comp]+";"
                    allCompProtCount+=str(len(complex2protDic[comp].split(",")))+";"
            else:
                AllComplex="NF"
                allComplexProt="NF"
                allCompProtCount="NF"
        else:
            AllComplex="NA"
            allComplexProt="NA"
            allCompProtCount="NA"
        wline=pairDic[eachPair]+[AllComplex,allComplexProt,allCompProtCount]
        output.write("\t".join(wline)+"\n")
        complexDic[eachPair]=wline
    output.close()

    print " total pairs found as complex are ",found, " from total ",total
    return complexDic,header

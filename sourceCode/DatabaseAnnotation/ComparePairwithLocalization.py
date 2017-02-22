__author__ = 'surya'

############################################################################################################################################
############################################################################################################################################
## check the pair with the GO cellular componenets

def DetermineCellularComponentForPair(prot2loc,pairdic,out,header="",location=False,locIndex=0):
    protDic={}
    outfile = open(out, "w")
    outfile.write(header+"\n")
    found=0
    for eachPair in pairdic:
        p1 = pairdic[eachPair][0]
        p2 = pairdic[eachPair][1]
        if p1 in prot2loc and p2 in prot2loc:
            if location:
                loc1=[loc.strip(" \"") for loc in prot2loc[p1][locIndex].split(";")]
                loc2=[loc.strip(" \"") for loc in prot2loc[p2][locIndex].split(";")]
            else:
                loc1=prot2loc[p1]
                loc2=prot2loc[p2]
            if loc1!="" and loc2!="":
                common=list(set(loc1).intersection(loc2))
                if len(common)!=0:
                    comlocName=";".join(common)
                else:
                    comlocName="NF"
            else:
                if loc1=="":
                    loc1=[]
                if loc2=="":
                    loc2=[]
                comlocName="NA"
                common=[]
            # print common
        else:
            comlocName = "NA"
            common = []
        wline=pairdic[eachPair]+[str(len(common)),comlocName]
        protDic[eachPair] = wline
        line2write="\t".join(wline) + "\n"
        outfile.write(line2write)


    outfile.close()
    return header,protDic

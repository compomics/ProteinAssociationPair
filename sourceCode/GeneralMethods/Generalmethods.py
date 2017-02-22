__author__ = 'surya'


### create dictionary
# with more than one val

## make dictionary
def createDic(file1,key,val_list,header,multiple=True,addDublicate=False,sep="\t",checkCol=False,CheckcolNumValue=[0,0],checkList=False):
    if multiple:
        if addDublicate:
            print("select only one of multiple or adddublicate, both cannot be used")
            return {}
    else:
        if isinstance(val_list,list):
            print ("with sigle entry please provide only one value not list")
            return {}

    dic={}
    with open(file1) as filename:
        if header:
            next(filename)
        for line in filename:
            if line!="\n":
                splits=line.split(sep)
                name=splits[key].strip()
                columnChecked=False
                if checkCol:
                    if checkList:
                        if splits[CheckcolNumValue[0]].strip() in CheckcolNumValue[1]:
                            columnChecked=True
                    else:
                        if splits[CheckcolNumValue[0]].strip() == CheckcolNumValue[1]:
                            columnChecked=True
                else:
                    columnChecked=True
                if columnChecked:
                    if name not in dic:
                        if multiple:
                            dic[name]=[]
                            for index in val_list:
                                dic[name].append(splits[index].strip())
                        else:
                            if addDublicate:
                                dic[name]=[splits[val_list].strip()]
                            else:
                                dic[name]=splits[val_list].strip()
                    else:
                        if addDublicate:
                            if splits[val_list].strip() not in dic[name]:
                                dic[name].append(splits[val_list].strip())
    return dic



## create a dic with two columns as the key rather than one

def createDicWith2Col(file1,key1,key2,val_list,header,multiple=True,addDublicate=False,sep="\t",
                      lineVal=False,checkCol=False,CheckcolNumValue=[0,0],checkList=False):
    """

    :param file1:
    :param key1:
    :param key2:
    :param val_list:
    :param header:
    :param multiple:
    :param addDublicate:
    :param sep:
    :param lineVal:
    :param checkCol:
    :param CheckcolNumValue: [columnNumber,ValueOfColumn]
    :return:
    """

    if multiple:
        if addDublicate:
            print("select only one of multiple or adddublicate, both cannot be used")
            return {}
    else:
        if isinstance(val_list,list):
            print ("with sigle entry please provide only one value not list")
            return {}

    dic={}
    with open(file1) as filename:
        if header:
            next(filename)
        for line in filename:
            splits=line.strip().split(sep)
            name=splits[key1].strip()+"_"+splits[key2].strip()
            columnChecked=False
            if checkCol:
                if checkList:
                    if splits[CheckcolNumValue[0]].strip() in CheckcolNumValue[1]:
                        columnChecked=True
                else:
                    if splits[CheckcolNumValue[0]].strip() == CheckcolNumValue[1]:
                        columnChecked=True
            else:
                columnChecked=True
            if columnChecked:
                if name not in dic:
                    if lineVal:
                        dic[name]=line.strip().split(sep)
                    else:
                        if multiple:
                            dic[name]=[]
                            for index in val_list:
                                dic[name].append(splits[index].strip())
                        else:
                            if addDublicate:
                                dic[name]=[splits[val_list].strip()]
                            else:
                                dic[name]=splits[val_list].strip()
                else:
                    if addDublicate and not lineVal:
                        if splits[val_list].strip() not in dic[name]:
                            dic[name].append(splits[val_list].strip())
                # else:
                #     print (name,"is a duplicate entry")

    return dic

def CreateDicWithAllColAsVal(file1,key,header,sep="\t",checkCol=False,CheckcolNumValue=[0,0],checkList=False):
    dic={}
    dublicate=0
    notfoundColumn=0
    columnChecked=False
    with open(file1) as filename:
        if header:
            next(filename)
        for line in filename:
            splits=line.strip().split(sep)
            name=splits[key].strip()
            if checkCol:
                if checkList:
                    if splits[CheckcolNumValue[0]].strip() in CheckcolNumValue[1]:
                        columnChecked=True
                else:
                    if splits[CheckcolNumValue[0]].strip() == CheckcolNumValue[1]:
                        columnChecked=True
            else:
                columnChecked=True
            if columnChecked:
                if name not in dic:
                    dic[name]=[]
                    for index in range(0,len(splits)):
                        if index!=key:
                            dic[name].append(splits[index].strip())
                else:
                    dublicate+=1
            else:
                notfoundColumn+=1
    print "total dublicate keys found are ",dublicate
    if checkCol:
        print "total lines where ",CheckcolNumValue[1], " was not found are ", notfoundColumn
    print
    return dic


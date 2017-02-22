__author__ = 'surya'

## import
import scipy.stats
import math
from math import factorial as fac


## to determine the correlation

def determinePearsonCorrelation(array1,array2):

    return scipy.stats.pearsonr(array1,array2) ## (pearson coefficient,2 tailed pvalue)


## to calculate the cosine correlation

def cosine_similarity(v1,v2):
    "compute cosine similarity of v1 to v2: (v1 dot v1)/{||v1||*||v2||)"
    sumxx, sumxy, sumyy = 0, 0, 0
    for i in range(len(v1)):
        x = v1[i]; y = v2[i]
        sumxx += x*x
        sumyy += y*y
        sumxy += x*y
    return sumxy/math.sqrt(sumxx*sumyy)

## to calculate a newScore

def ImprovedCosSim(v1,v2,commonExp):
    "compute cosine similarity of v1 to v2: (v1 dot v1)/{||v1||*||v2||)"
    sumxx, sumxy, sumyy = 0, 0, 0
    for i in range(len(v1)):
        x = v1[i]; y = v2[i]
        sumxx += x*x
        sumyy += y*y
        sumxy += x*y
    top=float(sumxy*commonExp*2)
    down1=float(sumxx*len(v1))
    down2=float(sumyy*len(v2))
    # if commonExp!=0:
    #     print float(top)/(down1+down2),top,down1,down2

    return top/(down1+down2)

def BionomialScore(n,p,k):
    C = fac(n) // fac(k) // fac(n - k)
    bionom=C*math.pow(p,k)*math.pow(1-p,n-k)
    return bionom


## calculate newscoring using probability

def ProbScore(common1,total1,common2,total2):
    return (float(common1)/total1)*(float(common2)/total2)


def jaccardSimilarity(list1,list2):
    minList,maxList=0,0
    if len(list1)!=len(list2):
        print "check the list length are not same"
    else:
        # print list1,list2
        for index in range(len(list1)):
            minList+=min(list1[index],list2[index])
            maxList+=max(list1[index],list2[index])
    return round(float(minList)/float(maxList),2)



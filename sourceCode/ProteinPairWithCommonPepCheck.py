__author__ = 'surya'


######################################################################################################################################################
##################### ## import methods and files and file pathways pathways ###################################################################
######################################################################################################################################################

import os
from Analysis import PathwayCoverage, createPairsByComparing
from GeneralMethods import Generalmethods
from DatabasePreprocessing import protein2pathway, SelectingIntActForSignallingProteins
from PreProcessing import filterings
from DatabaseAnnotation import AddParalogueInformation,CheckPairForComplex,RemovePairHavingCrapProteins, \
    GetDistinctivePeptidesNames, ComparePairwithLocalization,checkForBioGridInteraction
from datetime import datetime
from PostProcessing import GetUniprotName, GetKeratinAndIgAnnotation

PPI_MSPath="C:/Users/surya/Desktop/ProjectData/PPI_MSdata/"
DatabasePath="C:/Users/surya/Desktop/PPI_and_functional_interaction_among_species/Database_data/"

MSprotein2Assay2PepCountFile="RespinData/Pandey/prot2assay2DistinctPep.txt"## this is for pandey dataset
# MSprotein2Assay2PepCountFile="RespinData/protein2EachAssay2Distinct_PeptideName_OnlyHuman.txt"
prot2pathFile="Reactome/ver56_7Jun/UniProt2Reactome_lowestLevel_humanOnly_7Jun.txt"
pathwaysFile="Reactome/ver56_7Jun/ReactomeCompleterListOfPathway_AllSpecies_7June.txt" ## for all species
crapFile="cRAP/crapProt2Names1March.txt"
UniprotFile="uniprot/human/uniprot_20201human_reviewed_6Jun.tab"
# Uniprot2GOFile="uniprot/6JuneFile/749ProteinMapped2uniprotwithGOOntology_15Sep.tab"
Uniprot2GOFile="uniprot/6JuneFile/489UniqueProteinsSelectedOnlyGO_13Oct.tab"## this is for pandey dataset
Uniprot2GeneIdFile="uniprot/6JuneFile/HumanProtein6JuneSelectedInfo.txt"
IntActFile="IntAct/MergingIntAct8Feb16.txt"
BioGridFile="biogrid/BIOGRIDhuman34145tab2.txt" #onlyhuman

Uniprot2EnsemblFile="uniprot/human/UniprotAll_18963_Proteins2EnsemblGeneID_6Jun2016.txt"
# UniprotFile="database/"
paralogFile="ensembl/human/geneParalog15Feb.txt"
complexFile="corum_complexDatabase/coreComplexes_07Dec16.txt"
protein2PepetideFile="RespinData/Pandey/protein2distinctPeptide.txt" ## this is for pandey dataset
# protein2PepetideFile="RespinData/EachProtein2DistinctPepetides9March.txt"
interproFile="InterPro/InterproDomainSearch_20204humanReviewed_OnlyCanonicalFasta_15April.tsv"



######################################################################################################################################################
######################################  cutoff selections ############################################################################################
######################################################################################################################################################
start= datetime.now()
print start
minAssayCount=9
jaccCutoff=0.4
todayDate="16Feb"


######################################################################################################################################################
######################################  use cutoff selection to write files ##########################################################################
######################################################################################################################################################
newPath=PPI_MSPath+"CheckingMSPairForAllMSProteins/newDataFileteredDataset/pandey/UniprotWithRemovedTop10Assay/"
# if not os.path.isdir(newPath):
#     os.makedirs(newPath)
log=[str(datetime.now())+" starting time..."]
logFile = open(newPath + "logFile.txt", "w")
outWrite1=open(newPath+"UniProtPairWith0.1Jacc_"+todayDate+".txt","w")
outWrite2=open(newPath+"UniProtPairUsedForPathwayMappingwith"+str(minAssayCount+1)+"CommAssay"+str(jaccCutoff)+"jacc_"+todayDate+".txt","w")
outWrite3=open(newPath+"UniProtPairNotUsedForPathwayMappingwith"+str(minAssayCount+1)+"CommAssay"+str(jaccCutoff)+"jacc_"+todayDate+".txt","w")
outWrite4=open(newPath+"PathWithUniProtFound"+str(minAssayCount+1)+"CommAssay"+str(jaccCutoff)+"jacc_"+todayDate+".txt","w")
outWrite5=open(newPath+"proteinPair2With0.1_CommonAssayNames.txt","w")
selectedProteinFile=open(newPath+"SelectedProteinPresentInUniprot.txt","w")


paraOutname=newPath + "UniProteinAllPairWithParalogue"+str(jaccCutoff)+"jacc_"+todayDate+".txt"
complexOutname=newPath+"UniProteinAllPairWithParalogueComplex"+str(jaccCutoff)+"jacc_"+todayDate+".txt"
PepOutname=newPath+"UniProteinAllPairWithParalogueComplexPepCount"+str(jaccCutoff)+"jacc_"+todayDate+".txt"
crapOut=newPath+"UniProteinAllPairWithParalogueComplexPepCountCrapLabel"+str(jaccCutoff)+"jacc_"+todayDate+".txt"
uniNamesOut=newPath+"UniProtAllPairWithParaCompxCrapName"+str(jaccCutoff)+"jacc_"+todayDate+".txt"
BiogridOut=newPath+"UniAllPairWithParaCompxCrpNameBiogrid"+str(jaccCutoff)+"jacc_"+todayDate+".txt"
keratinOut=newPath+"UniAllPairWithParaCompxCrpNameBGkerIg"+str(jaccCutoff)+"jacc_"+todayDate+".txt"
bp=newPath+"UniAllPairWithParaCompxCrpNameBGKerIgIntBP"+str(jaccCutoff)+"jacc_"+todayDate+".txt"
mf=newPath+"UniAllPairWithParaCompxCrpNameBGKerIgIntBPMF"+str(jaccCutoff)+"jacc_"+todayDate+".txt"
cc=newPath+"UniAllPairWithParaCompxCrpNameBGKerIgIntBPMFCC"+str(jaccCutoff)+"jacc_"+todayDate+".txt"


######################################################################################################################################################
###################################### tissue based  protein filtrations #############################################################################
######################################################################################################################################################

## get all the crapProteins lists

CrapProtIdsDic=Generalmethods.createDic(DatabasePath+crapFile,1,0,header=False,multiple=False)
print "total crap proteins are ", len(CrapProtIdsDic)

## get all the uniprot proteins and filter the mS proteins
UniprotDic=Generalmethods.createDic(DatabasePath+UniprotFile,0,val_list=[1,3],header=True)
#get protein to geneID mapping
Uni2GenIDDic=Generalmethods.createDic(DatabasePath+Uniprot2GeneIdFile,5,val_list=[0,2,3,10],header=True)

print "total uniprot proteins are ", len(UniprotDic)

## get the GO ontologies
Uni2GODic=Generalmethods.createDic(DatabasePath+Uniprot2GOFile,0,[2,3,4],header=True)

## get the dic for protein 2 Assay 2 pepcount
FilProt2Ass2PepNameListDic,AssayList,protein2AssayPepString,logP=filterings.CheckAssayWithPeptideIds\
                                                (PPI_MSPath+MSprotein2Assay2PepCountFile,AssayLen=minAssayCount,
                                                           proteincheck=True,proteinList=UniprotDic,MinPep=1,
                                                 selectedProteinFile=selectedProteinFile)
log+=[logP]

## get all intact pair for all proteins
IntActPair= SelectingIntActForSignallingProteins.selectPairsIntAct(DatabasePath+IntActFile,FilProt2Ass2PepNameListDic) # dic with protein pairs and score as value

######################################################################################################################################################
###################################### pathway selection #############################################################################################
######################################################################################################################################################

## get all pathways with there names

PathId2NameDic=Generalmethods.createDic(DatabasePath+pathwaysFile,0,1,header=False,multiple=False,checkCol=True,CheckcolNumValue=[2,"Homo sapiens"])
ln=str(len(PathId2NameDic))+"total pathways present\n"
log.append(ln)
print ln

## get all the proteins for the leaf pathway and dic from prot to pathway
prot2pathDic,ReactomeDic= protein2pathway.dicFromProtein2Pathway(DatabasePath+prot2pathFile)#,UniprotDic) ## checking with all selected pathways

ln="total protein found in hierarchy are "+str(len(prot2pathDic))+"\n"
log.append(ln)
print ln

######################################################################################################################################################
######################################  Run protein co-occurrence  #############################################################################################
######################################################################################################################################################


# select the pair which have atleast one com assay with more than 2 peptide
## determine the correaltion btw the proteins
pairList,protAfterselection,logP2=createPairsByComparing.createPairsWithPepCheck(FilProt2Ass2PepNameListDic,minAssayCount,
                                                    IntActPair,jaccCutoff,outWrite1,1,protein2AssayPepString,outWrite5)
log += [logP2]

# check for how many proteins have their pathway present
protWithPath=list(set(protAfterselection.keys()).intersection(prot2pathDic.keys()))
ln = "total protein with pathways are " + str(len(protWithPath)) + " from total " + str(len(protAfterselection)) + "\n"
print ln
log.append(ln)

print datetime.now()-start,"total time for finding pairs"
######################################################################################################################################################
######################################  check for biological relevance ###############################################################################
######################################################################################################################################################

## check the pair with the pathways
pair2pathTags,pathHeader=PathwayCoverage.checkTheCoverageForLeafPathway(prot2pathDic,pairList,outWrite4,outWrite2,outWrite3,PathId2NameDic,prot2PathDicAvailable=True,ReactDic=ReactomeDic)



## first look for the proteins in uniprot if they have chaged names
## add paralogue information

paraTaggedPair,Pheader=AddParalogueInformation.findParalogue(DatabasePath+paralogFile,DatabasePath+Uniprot2EnsemblFile,
                                                             pair2pathTags,paraOutname,pathHeader)#,old2newProtName)

print "Done with paralog tagging"

## than determine if they are part of any complex
#
complexPairdDic,CHeader=CheckPairForComplex.tagComplexPair(DatabasePath+complexFile,paraTaggedPair,complexOutname,Pheader)

print "done with complex tagging"
## before using further check if proteins names are changed in uniprot or protein which are removed

## than select common peptide between each pair of proteins

PeptideDic,PpHeader=GetDistinctivePeptidesNames.tagPeptideCount(complexPairdDic,PPI_MSPath+protein2PepetideFile,PepOutname,CHeader)
print "done with peptide tagging"

## labelle crap proteins
crapPairDic,cheader=RemovePairHavingCrapProteins.tagCraps(CrapProtIdsDic,PeptideDic,crapOut,PpHeader)

## add protein name
UniprotNameDic,Uheader=GetUniprotName.AddProteinName(crapPairDic,UniprotDic,uniNamesOut,cheader)

## get the protein interaction from Biogrid
biogridPairDic,bioHeader=checkForBioGridInteraction.checkForProteinList(UniprotNameDic,BiogridOut,Uni2GenIDDic,
                                                                        DatabasePath+BioGridFile,Uheader)
## get the Ig and keratin annotation
keratinPairDic,kerHeader=GetKeratinAndIgAnnotation.labelPairWithFilteredProteins(biogridPairDic,keratinOut,bioHeader)

## get GO
bpheader=kerHeader+"\tBPComCont\tBPCom"
bpheadOut,bpPairDic=ComparePairwithLocalization.DetermineCellularComponentForPair(Uni2GODic,keratinPairDic,bp,bpheader,location=True,locIndex=0)
print "finished writing BP",len(bpPairDic)
mfheader=bpheader+"\tMFComCont\tMFCom"
mfheadOut,mfPairDic=ComparePairwithLocalization.DetermineCellularComponentForPair(Uni2GODic,bpPairDic,mf,mfheader,location=True,locIndex=1)
print "finished writing mf",len(mfPairDic)
ccheader=mfheader+"\tCCComCont\tCcCom"
ccheadOut,ccPairDic=ComparePairwithLocalization.DetermineCellularComponentForPair(Uni2GODic,mfPairDic,cc,ccheader,location=True,locIndex=2)
print "finished writing cc",len(ccPairDic)

ln=str(datetime.now()-start) +" total time"
log.append(ln)

for eachL in log:
    logFile.write(eachL)
logFile.close()

print datetime.now()-start,"total time"
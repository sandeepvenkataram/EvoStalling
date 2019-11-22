library(ggplot2)
library(plyr)
library(reshape2)
library(png)
library(dplyr)
library(gridExtra)
library(egg)
library(RColorBrewer)
##############################
#
# Global Variables
#
##############################
dataInputDir = "inputData/"
outputDir = "PlotsAndTables/"

sourceOrganisms = c("E","S","V","Y","P","A")
sourceOrganismsOrdered = c("E","S","Y","V","A","P")
populationTypes = c("High","High","High","Low","Low","Low")
possibleMutTypes = c("Chr amp","NCOD","NCRNA","NSYN","STOPGAIN","STOPLOSS","SYN")
renamedMutTypes = c("Chr amp","Noncoding","Other","Nonsynonymous","Other","Other","Synonymous")
renamedMutTypesOrder = c("Nonsynonymous","Synonymous","Other","Chr amp","Noncoding")
mutImportanceOrder = c("Chr amp","FRAMESHIFT","STOPGAIN","STOPLOSS","NSYN","SYN","NCOD","NCRNA","")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#CC79A7","#999999")
speciesColors = cbPalette[1:6]
names(speciesColors) <- sourceOrganismsOrdered


replicateColors <- c("#117733","#88CCEE","#DDCC77")

allMutFileName <- paste(dataInputDir,"allMuts.tab",sep="")
selectedMutFileName <- paste(dataInputDir,"selectedMuts.tab",sep="")
lowSelectedMutFileName <- paste(dataInputDir,"selectedMuts_low.tab",sep="")
ultraLowSelectedMutFileName <- paste(dataInputDir,"selectedMuts_ultralow.tab",sep="")


translationGenes = c("accD", "acnB", "acs", "adk", "alaS", "arfA", "arfB", "argS", "asnS", "aspS", "bipA", "cheA", "cheW", "csrA", "cysN", "cysS", "dbpA", "deaD", "def", "der", "dtd", "efp", "epmA", "epmB", "epmC", "era", "ettA", "fmt", "frr", "fusA", "glnS", "gltX", "gluQ", "gshB", "hflX", "hfq", "higB", "hisS", "hpf", "hslR", "iap", "ileS", "infA", "infB", "infC", "lepA", "leuS", "lhr", "lit", "lysS", "lysU", "metG", "mgtL", "miaA", "mqsR", "obgE", "pheS", "pheT", "prfA", "prfB", "prfC", "prfH", "prmA", "prmB", "prmC", "proS", "raiA", "ratA", "rbbA", "rbfA", "relE", "rfaH", "rhlE", "ridA", "rimI", "rimJ", "rimK", "rimL", "rimM", "rimO", "rimP", "rlmA", "rlmB", "rlmC", "rlmD", "rlmE", "rlmF", "rlmG", "rlmH", "rlmI", "rlmJ", "rlmL", "rlmM", "rlmN", "rluA", "rluB", "rluC", "rluD", "rluE", "rluF", "rmf", "rnc", "rne", "rng", "roxA", "rph", "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplI", "rplJ", "rplK", "rplL", "rplM", "rplN", "rplO", "rplP", "rplQ", "rplR", "rplS", "rplT", "rplU", "rplV", "rplW", "rplX", "rplY", "rpmA", "rpmB", "rpmC", "rpmD", "rpmE", "rpmF", "rpmG", "rpmI", "rpmJ", "rpsA", "rpsB", "rpsC", "rpsD", "rpsE", "rpsF", "rpsG", "rpsH", "rpsI", "rpsJ", "rpsK", "rpsL", "rpsM", "rpsN", "rpsO", "rpsP", "rpsQ", "rpsR", "rpsS", "rpsT", "rpsU", "rsfS", "rsgA", "rsmA", "rsmB", "rsmC", "rsmD", "rsmE", "rsmF", "rsmG", "rsmH", "rsmI", "rsmJ", "rsuA", "rutC", "secM", "selB", "serS", "smpB", "sra", "srmB", "sspB", "tar", "tdcF", "thrS", "tig", "tnaC", "trpS", "truC", "tsaC", "tsf", "tufA", "tufB", "typA", "tyrS", "uof", "valS", "yafO", "yafQ", "yajL", "ybaK", "ybcJ", "ybeY", "ycaO", "yceD", "ychF", "yciH", "yeaK", "yeiP", "yhaV", "yhbY", "yhgF", "yihI", "ykgM", "yoeB", "yqgF", "yqjD") # from geneontology.org, searching for all E. coli genes that came up in searches for "translation" and "ribosom"

motilityGenes = c("motA", "fliE", "frdC", "glgS", "malE", "fliS", "motB", "dgcJ", "cdgI", "ydiV", "cheY", "recA", "pyrR", "yggR", "pdeH", "fliH", "ycgR", "fliG", "fliT", "flgG", "flgC", "flgB", "hofN", "mqsR", "fliF", "flgL", "frdB", "ybjN", "yfiR", "dgcC", "frdA", "fliO", "dgcQ", "tsr", "plaP", "yfgJ", "yciG", "uspE", "flgJ", "flgF", "flgD", "flgA", "ycdY", "ycdX", "fdrA", "dgcI", "cheA", "cheZ", "flgK", "fliM", "flgE", "fliC", "fliJ", "uspG", "ecpR", "chaC", "fliN", "pgrR", "dgcM", "fliL", "dgcN", "fliZ", "fliK", "fliI", "fliD", "dosC", "dgcF", "dgcZ", "flgI", "flgH") #"motility"

cytokinesisGenes = c("ftsZ", "pbpG", "minD", "ftsP", "slmA", "zapB", "ftsL", "dapE", "zapA", "ftsQ", "ftsI", "ydaT", "ftsA", "amiB", "sulA", "ftsW", "queE", "cbtA", "yohD", "dedD", "ygeR", "yqjA", "zipA", "yabI", "envC", "sdiA", "zapC", "ftsN", "ytfB", "minC", "damX", "mepM", "dacB", "mreB", "yghB", "amiC", "cpoB", "nlpD", "ftsB", "ftsK", "minE", "zapD", "yihA", "yciB", "rlpA") #"cytokinesis"

#manually called chromosomal amplifications of tufA
chrAmptufAPopulations = c("V_1","V_3","V_5","P_4","P_6","A_1","A_2","A_3","A_4","A_5","A_6")
chrAmptufAFirstObservedTimepoint = c(100,700,600,100,200,100,600,100,100,100,100)
chrAmptufALastObservedTimepoint = c(1000,1000,800,1000,1000,1000,900,1000,1000,1000,1000)
chrAmptufACopyNumber = c(4,2,2,35,6,2,2,2,2,2,2)
chrAmptufAGenes = rep("tufA",NROW(chrAmptufAPopulations))
chrAmptufAmutType = rep("Chr amp",NROW(chrAmptufAPopulations))
chrAmpDF = data.frame(popTitle = gsub("_","",chrAmptufAPopulations), xmin = chrAmptufAFirstObservedTimepoint, xmax = chrAmptufALastObservedTimepoint, copynum = chrAmptufACopyNumber, gene="tufA", mutType = "Chr amp")

varWeightedAverage = function(meanVals,varVals){
  
  newMean = sum(meanVals/varVals)/sum(1/varVals)
  newVar = 1/sum(1/varVals)
  return(c(newMean,newVar))
  
}


###############################
#
# Calculate Fitness from competition assays
#
###############################

inputFile = read.table(paste(dataInputDir,"Gen1000Competition_lnFreqRatios.tab",sep=""),sep="\t",header=TRUE)
ancestorFile = read.table(paste(dataInputDir,"AncCompetition_lnFreqRatios.tab",sep=""),sep="\t",header=TRUE)
reconstructedFile = read.table(paste(dataInputDir,"ReconstructedStrainCompetition_lnFreqRatios.tab",sep=""),sep="\t",header=TRUE)
EvsWTFile = read.table(paste(dataInputDir,"E_vs_WT_Competition_lnFreqRatios.tab",sep=""),sep="\t",header=TRUE)



getFitnessData <- function(myInputFile){
	fitnessDF = data.frame(Name = character(),Fitness = numeric(), stderror = numeric(), variance = numeric())
	myExperiments = as.character(unique(myInputFile$Experiment))
	myXVals = as.character(colnames(myInputFile))[5:NCOL(myInputFile)]
	myNumericXVals = as.numeric(gsub('t','',myXVals)) * log2(10000) #10,000 generations per growth cycle, so convert to generations

	for(experiment in myExperiments){
		reducedInputFile = myInputFile[myInputFile$Experiment == experiment,]
		reducedInputFileMelted = melt(reducedInputFile,id=c("Flask","Experiment","Marker","Replicate"))
		reducedInputFileMelted$variable = as.numeric(gsub('t','',as.character(reducedInputFileMelted$variable)))
		reducedInputFileMelted = reducedInputFileMelted[!is.na(reducedInputFileMelted$value),]
		
		experimentfitnessDF = data.frame(Name = character(),Fitness = numeric(), stderror = numeric(), variance = numeric())
		
		#calculate Fitness for each flask
		for (flask in reducedInputFile$Flask){
			flaskInput = reducedInputFile[reducedInputFile$Flask == flask,]
			myModel = lm(as.numeric(flaskInput[5:NCOL(flaskInput)])~myNumericXVals)
			myModelSummary = summary(myModel)
			sEstimate = myModelSummary$coefficients[2]
			stderrEstimate = myModelSummary$coefficients[4]
			varianceEstimate = (stderrEstimate * sqrt(myModelSummary$df[2]+1))**2 #back calculate variance from stderr and df
			experimentfitnessDF = rbind(experimentfitnessDF,data.frame(Name=as.character(flask),Fitness = sEstimate, stderror = stderrEstimate, variance = varianceEstimate))
		}
		
		#calculate Fitness overall by weighted averaging flask
		
		estimatesToUse = experimentfitnessDF[!is.na(match(experimentfitnessDF$Name,reducedInputFile$Flask)),]
		
		newVariance = estimatesToUse$variance[1]
		newAverage = estimatesToUse$Fitness
		newStdErr = estimatesToUse$stderror
		if(NROW(estimatesToUse)>1){
			#if there is missing variance data when trying to average, set the variance to the average of the other variances (e.g. if there are only 2 valid timepoints of data for that replicate)
			estimatesToUse$variance[is.na(estimatesToUse$variance)] = mean(estimatesToUse$variance,na.rm=TRUE)
			weightedAvg = varWeightedAverage(estimatesToUse$Fitness, estimatesToUse$variance)
			newVariance = weightedAvg[2]
			newAverage = weightedAvg[1]
			newStdErr = sqrt(newVariance) / sqrt(NROW(estimatesToUse)-1)
		}
		experimentfitnessDF = rbind(experimentfitnessDF,data.frame(Name=paste(experiment,"allFlasks_averaged",sep="_"),Fitness = newAverage, stderror = newStdErr, variance = newVariance))

		fitnessDF = rbind(fitnessDF,experimentfitnessDF)
	}
	return(fitnessDF)
}

fitnessDF = getFitnessData(inputFile)
write.table(fitnessDF,file=paste(outputDir,"rawFitnessData.tab",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

ancFitnessDF = getFitnessData(ancestorFile)

#The expected fitness value of competition between two E strains is 0. The observed fitness is used to calculate variance of the fitness estimation
ancFitnessDF$Fitness[ancFitnessDF$Name == "E vs E_allFlasks_averaged"]=0
ancFitnessDF$variance[ancFitnessDF$Name == "E vs E_allFlasks_averaged"]=mean((ancFitnessDF$Fitness[1:4])^2)
ancFitnessDF$stderror[ancFitnessDF$Name == "E vs E_allFlasks_averaged"]=sqrt(mean((ancFitnessDF$Fitness[1:4])^2))

write.table(ancFitnessDF,file=paste(outputDir,"ancFitnessData.tab",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

reconstructedFitnessDF = getFitnessData(reconstructedFile)
write.table(reconstructedFitnessDF,file=paste(outputDir,"reconstructedFitnessData.tab",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

EvsWTFitnessDF = getFitnessData(EvsWTFile)
write.table(EvsWTFitnessDF,file=paste(outputDir,"EvsWTFitnessData.tab",sep=""),sep="\t",row.names=FALSE,quote=FALSE)




###############################
#
# Extract and filter timecourse data
#
###############################

timecourseDataDir <-paste(dataInputDir,"annotated_timecourse_files/",sep="")
thr_cnt = 10
thr_fixfreq = 0.95
thr_detectfreq = 0.1

initialTimecourseDataProcessing <- function(outputFileName){
  numRows = 10000
  numPops = NROW(sourceOrganisms)*6
  validPops = c()
  numTimepoints = 11
  timepointArray = c(0:10)*100
  
  readCounts <- array(NA,dim=c(numRows,numTimepoints,numPops))
  readDepths <- array(NA,dim=c(numRows,numTimepoints,numPops))
  readFreqs <- array(NA,dim=c(numRows,numTimepoints,numPops))
  
  mutPos <- rep(NA,numRows)
  mutChange = rep(NA,numRows)
  mutLocationString = rep(NA,numRows)
  mutAnnotationString = rep(NA,numRows)
  curMutID = 1
  popNumber = 1
  
  ### import data and initial filtering
  for(founderStrain in sourceOrganismsOrdered){
    for(cultureNumber in 1:6)  {
      myPop = paste(founderStrain,cultureNumber,sep="_")
      validPops = c(validPops, myPop)
      infile<-file(paste(timecourseDataDir,myPop,"_merged_timecourse_annovar_withAnnotations.txt",sep=""),"r")
      while(TRUE){
        myLine = readLines(infile,n=1)
        if(length(myLine)==0){
          break
        }
        mySplit = strsplit(myLine,",\t")[[1]]
        myTimepoints = as.numeric(strsplit(mySplit[4]," ")[[1]])
        myReadCounts = as.numeric(strsplit(mySplit[5]," ")[[1]])
        myReadDepths = as.numeric(strsplit(mySplit[6]," ")[[1]])
        myReadFreqs = myReadCounts/myReadDepths
        myNumTimepoints = NROW(myTimepoints)
        
        if(sum(myReadCounts)<thr_cnt){ #skip mutation if it has too few reads
          next
        }
        
        testVals = abs(.5-myReadFreqs[(myNumTimepoints-1):myNumTimepoints]) > (thr_fixfreq - .5)
        if(myReadFreqs[1]>0 & sum(testVals,na.rm=TRUE) < sum(!is.na(testVals)) ){ #mutations present in the beginning should fix or be lost by gen 900 / 1000, remove if this is not the case
          next
        }
        if(sum(myReadCounts > 0 ) <=1){ #skip mutation if it is only present in one timepoint
          next
        }
        
        myMutString = paste(mySplit[2],mySplit[3])
        myMutID = curMutID
        mutIDLoc = which(mutLocationString == myMutString)
        if(NROW(mutIDLoc)==1){
          myMutID = mutIDLoc
        }else{
          mutChange[curMutID] = mySplit[3]
          mutLocationString[curMutID] = myMutString
          if(NROW(mySplit)>6){
            mutAnnotationString[curMutID] = mySplit[7]
          }
          
          curMutPos = as.numeric(mySplit[2]) # lift over coordinates to be compatible with reference MG1655 genome
          if(curMutPos>4175944){ #We are missing tufB in the modified genomes. Add it back in for the lifted over coordinates
            curMutPos = curMutPos + 1185
          }
          if(curMutPos>3470744 & founderStrain == "P"){ # The pseudomonas tufA is 9bp longer than every other tufA due to a 3AA insertion around AA 195. 
            print("here!")
            curMutPos = curMutPos - 9
          }
          mutPos[curMutID] = curMutPos
          curMutID = curMutID + 1
        }
        readCounts[myMutID,(myTimepoints/100 + 1),popNumber] = myReadCounts
        readDepths[myMutID,(myTimepoints/100 + 1),popNumber] = myReadDepths
        readFreqs[myMutID,(myTimepoints/100 + 1),popNumber] = myReadFreqs
      }
      close(infile)
      popNumber = popNumber + 1
    }
  }
  
  ### bad timepoints to remove from analysis
  
  badPops = c("E_1","E_2","Y_3","P_3","P_2","P_2","A_1")
  badTimes = c(600,500,600,600,800,1000,700)
  
  removeBadTimes <- function(pop, tp){
    readCounts[,tp,pop] <<- NA
    readDepths[,tp,pop] <<- NA
    readFreqs[,tp,pop] <<- NA
  }
  
  for(i in c(1:NROW(badPops))){
    removeBadTimes(which(validPops==badPops[i]), badTimes[i]/100+1)
  }
  
  
  ### identify ancestral and founder-specific mutations
  
  ancMuts = apply(readFreqs[,1,],1,function(x) sum(x,na.rm=TRUE)/18) > thr_fixfreq
  founderMuts = rep(FALSE,numRows)
  founderMuts[which(apply(readFreqs[,1,1:6],1,function(x) sum(x,na.rm=TRUE)/6) > thr_fixfreq)]=TRUE
  founderMuts[which(apply(readFreqs[,1,7:12],1,function(x) sum(x,na.rm=TRUE)/6) > thr_fixfreq)]=TRUE
  founderMuts[which(apply(readFreqs[,1,13:18],1,function(x) sum(x,na.rm=TRUE)/6) > thr_fixfreq)]=TRUE
  founderMuts[which(apply(readFreqs[,1,19:24],1,function(x) sum(x,na.rm=TRUE)/6) > thr_fixfreq)]=TRUE
  founderMuts[which(apply(readFreqs[,1,25:30],1,function(x) sum(x,na.rm=TRUE)/6) > thr_fixfreq)]=TRUE
  founderMuts[which(apply(readFreqs[,1,31:36],1,function(x) sum(x,na.rm=TRUE)/6) > thr_fixfreq)]=TRUE
  founderMuts = founderMuts & ! ancMuts
  
  print("Number of ancestral muts")
  print(sum(ancMuts))
  print("Number of founder-specific muts")
  print(sum(founderMuts))
  #mutations that are present at any frequency at t100 in 10+ populations are likely mapping artifacts
  badMuts = apply(readFreqs[,2,],1,function(x) sum(x>0,na.rm=TRUE)) > 11
  print("Number of muts present at any freq at t100 in 11+ populations")
  print(sum(badMuts))
  
  goodMuts = array(FALSE,c(numRows,numPops)) #check to make sure a mutation is present at some minimum frequency in at least one timepoint and detected in at least 2 consecutive timepoints, otherwise it is likely an artifact
  checkFreqTraj <- function(x){ 
    t2 = sum(myReadCounts[1:(myNumTimepoints-1)] > 0 & myReadCounts[2:myNumTimepoints] > 0, na.rm=TRUE) >0 #check if mutation is present in at least 2 consecutive timepoints
    #check if mutation over a threshold frequency in at least one timepoint
    numTOver10 = sum(x > thr_detectfreq,na.rm=TRUE)
    return(t2 & numTOver10 > 0)
  }
  
  for(popNumber in c(1:numPops)){
    goodMuts[,popNumber] = apply(readFreqs[,,popNumber],1, checkFreqTraj)
  }
  
  validMutIDs = which(apply(goodMuts,1,sum)>0 & (! ancMuts) & (! founderMuts) & ! grepl("junction",mutLocationString) & !badMuts)
  
  
  ### write mutations to data table. ignore junctions as they are impossible to quantify frequency accurately and are frequently artifacts
  mutationLines = c()
  for(mutID in validMutIDs){
    myMutPos = gsub("\t","",mutPos[mutID])
    myMutChange = gsub("\t","",mutChange[mutID])
    annoSplit = strsplit(gsub("\t","",gsub(", MG1655","",mutAnnotationString[mutID])),":")[[1]]
    
    myMutType = ""
    myCDSPos = NA
    myAAPos = NA
    myAAChange = ""
    myGene = ""
    myUpstreamGene = ""
    myDownstreamGene = ""
    
    #if nonsynonymous or synonymous
    if(grepl("synonymous", annoSplit[1])){
      myMutType = "SYN"
      myGene = gsub(".*_","",annoSplit[2])
      annoSplit[4] = gsub("c.","",annoSplit[4])
      annoSplit[5] = gsub("p.","",annoSplit[5])
      annoSplit[5] = gsub(", .*","",annoSplit[5])
      myCDSPos = as.numeric(gsub("[A-Za-z]","",annoSplit[4]))
      myAAPos = as.numeric(gsub("[A-Za-z]","",annoSplit[5]))
      myAAChange = gsub("[0-9]+","->",annoSplit[5])
    }
    
    if(grepl("nonsynonymous", annoSplit[1])){
      myMutType = "NSYN"
    }
    #if stop gain / loss
    if(grepl("stop", annoSplit[1])){
      myMutType = toupper(gsub(" .*","",annoSplit[1]))
      myGene = gsub(".*_","",annoSplit[2])
      annoSplit[4] = gsub("c.","",annoSplit[4])
      annoSplit[5] = gsub("p.","",annoSplit[5])
      myCDSPos = as.numeric(gsub("[A-Za-z]","",annoSplit[4]))
      myAAPos = as.numeric(gsub("[A-Za-z]","",annoSplit[5]))
      myAAChange = gsub("[0-9]+","->",annoSplit[5])
    }
    
    
    #if frameshift
    
    if(grepl("frameshift", annoSplit[1])){
      myMutType = "FRAMESHIFT"
      myGene = gsub(".*_","",annoSplit[2])
      annoSplit[4] = gsub("c.","",annoSplit[4])
      annoSplit[5] = gsub("p.","",annoSplit[5])
      myCDSPos = as.numeric(gsub("[A-Za-z_]","",annoSplit[4]))
    }
    
    #if upstream or downstream
    
    if(grepl("stream", annoSplit[1])){
      myMutType = "NCOD"
      annoSplit = strsplit(gsub("\t","",mutAnnotationString[mutID])," ")[[1]]
      locationSplit = strsplit(annoSplit[1],";")[[1]]
      locusSplit = strsplit(annoSplit[2],";")[[1]]
      if(grepl("downstream", annoSplit[1])){
        downGeneStr = locusSplit[locationSplit=="downstream"]
        myDownstreamGene = gsub(".*_","",gsub(",.*","",downGeneStr))
        myGene = myDownstreamGene
      }
      
      if(grepl("upstream", annoSplit[1])){ 
        upGeneStr = locusSplit[locationSplit=="upstream"]
        myUpstreamGene = gsub(".*_","",upGeneStr)
        myGene = myUpstreamGene #annotate gene as upstream gene, even if downstream gene exists, since it is more likely to regulate the gene it is upstream of
      }
    }
    
    #if RNA
    if(grepl("RNA", annoSplit[1])){
      myMutType = "NCRNA"
      myGene = gsub(".*_","",annoSplit[2])
    }
    
    for(validPopID in which(goodMuts[mutID,])){
      freqTraj = readFreqs[mutID,,validPopID]
      depthTraj = readDepths[mutID,,validPopID]
      mutationLines = c(mutationLines, paste(myGene, mutID, myMutPos, myMutChange, myMutType, myCDSPos, myAAPos, myAAChange, myUpstreamGene, myDownstreamGene, validPops[validPopID], paste(freqTraj,collapse="\t"), paste(depthTraj,collapse="\t"), sep="\t"))
    }	
  }
  outputFileHandle <- file(allMutFileName,"w") #write this unfiltered set to file temporarily so we can read it back as a table
  writeLines(mutationLines,outputFileHandle)
  close(outputFileHandle)
  
  allMutFile<-read.table(allMutFileName,sep="\t",header=FALSE)
  
  
  #merge mutations that are too close together as they are likely from a single event
  for(myPop in unique(allMutFile$V11)){
    myMutPosArray = sort(as.numeric(allMutFile$V3[allMutFile$V11==myPop]))
    diffTest = myMutPosArray[1:NROW(myMutPosArray)-1] - myMutPosArray[2:NROW(myMutPosArray)] > -10
    if(sum(diffTest)>0){
      curGroupStart = -1
      linesToRemove = c()
      for(i in c(1:NROW(diffTest))){
        if(diffTest[i]){
          if(curGroupStart ==-1){ #if this is a new run of mutations, start
            curGroupStart = i
          }
          if(i+1 <= NROW(diffTest) & diffTest[i+1]){ #if there are more mutations in the group, just continue
            next
          }else{ #else we are done, fix the mutation list
            myPositionsInGroup = myMutPosArray[curGroupStart:(i+1)]
            validLines = which(allMutFile$V11 == myPop & allMutFile$V3 %in% myPositionsInGroup)
            matchVals = match(allMutFile$V5[validLines],mutImportanceOrder)
            bestMatch = which(matchVals == min(matchVals))[1]
            linesToRemove = c(linesToRemove,validLines[-bestMatch])
            curGroupStart = -1
          }
        }
      }
      if(NROW(linesToRemove)>0){
        allMutFile = allMutFile[-linesToRemove,]
      }
    }
  }
  allMutFile = allMutFile[! c(1:NROW(allMutFile)) %in% grep(";.*;",allMutFile$V4), ] #remove mutations that are more than biallelic, as they are likely mapping artifacts
  
  write.table(allMutFile,file=allMutFileName,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}
initialTimecourseDataProcessing()

allMutFile<-read.table(allMutFileName,sep="\t",header=FALSE)
selectedMuts = allMutFile[apply(allMutFile[,12:22],1,function(x){max(x,na.rm=TRUE)-min(x,na.rm=TRUE)>0.2}),] #selected mutations are those that changed at least 20% frequency across the experiment
selectedMutsLowThresh = allMutFile[apply(allMutFile[,12:22],1,function(x){max(x,na.rm=TRUE)-min(x,na.rm=TRUE)>0.1}),] #to see the robustness of this cutoff, we lower it to 10%

write.table(selectedMuts, file=selectedMutFileName,sep="\t",col.names=FALSE,quote=FALSE,row.names=FALSE)
write.table(selectedMutsLowThresh, file=lowSelectedMutFileName,sep="\t",col.names=FALSE,quote=FALSE,row.names=FALSE)





###############################
#
# Figure 1, S1 and S2
#
###############################

fitnessDF = read.table(paste(outputDir,"rawFitnessData.tab",sep=""),sep="\t",header=TRUE)
ancFitnessDF = read.table(paste(outputDir,"ancFitnessData.tab",sep=""),sep="\t",header=TRUE) #anc fitness is calculating fitness of E strain relative to the different founders, so need to negate the fitness estimates to get founder fitness relative to E

ancFitnessDFTrunc = ancFitnessDF[grep("allFlasks_averaged", ancFitnessDF$Name),]
fitnessDFTrunc = fitnessDF[grep("allFlasks_averaged", fitnessDF$Name),]

ancGrowthRateDF = read.table(paste(dataInputDir,"AncGrowthRateEstimates.txt",sep=""),sep="\t",header=TRUE)
ancGrowthRateDF$SEMGR = ancGrowthRateDF$Stdev.Growth.Rate/sqrt(3)


#convert fitness estimates to percent per generation
ancFitnessDFTrunc$Fitness = ancFitnessDFTrunc$Fitness*100
ancFitnessDFTrunc$stderror = ancFitnessDFTrunc$stderror*100
ancFitnessDFTrunc$variance = ancFitnessDFTrunc$variance*100^2

fitnessDFTrunc$Fitness = fitnessDFTrunc$Fitness*100
fitnessDFTrunc$stderror = fitnessDFTrunc$stderror*100
fitnessDFTrunc$variance = fitnessDFTrunc$variance*100^2


fitnessFounderList = gsub("_.*","",fitnessDFTrunc$Name)
fitnessFounderListFactor = factor(fitnessFounderList, levels = sourceOrganismsOrdered)

ancFitnessFounderList = gsub(" .*","",ancFitnessDFTrunc$Name)
ancFitnessFounderListFactor = factor(ancFitnessFounderList, levels = sourceOrganismsOrdered)
ancFitnessDFTrunc$Founder = ancFitnessFounderList

# calculate average fitness gain across all populations derived from each founder
mydf = data.frame(species = character(), evoFitness = numeric(), evoFitnessSEM=numeric(),ancFitness=numeric(), ancFitnessSEM = numeric())
for(organism in sourceOrganismsOrdered){
  observedFitnessVals = fitnessDFTrunc$Fitness[fitnessFounderList == organism]
  observedFitnessVarVals = fitnessDFTrunc$variance[fitnessFounderList == organism]
  weightedFitnessAndVar = varWeightedAverage(observedFitnessVals,observedFitnessVarVals)
  mydf = rbind(mydf, data.frame(species = organism, evoFitness = weightedFitnessAndVar[1],evoFitnessSEM = sqrt(weightedFitnessAndVar[2])/sqrt(sum(fitnessFounderList == organism)), ancFitness = ancFitnessDFTrunc$Fitness[ancFitnessDFTrunc$Founder==organism], ancFitnessSEM = ancFitnessDFTrunc$stderror[ancFitnessDFTrunc$Founder==organism]))
}
mydf$species = factor(mydf$species, levels = sourceOrganismsOrdered)

p1<-ggplot(mydf)+geom_point(aes(-ancFitness,evoFitness,col=species),size=3,alpha=.7)+scale_y_continuous(breaks =c(0,10,20,30,40,50))+theme_bw()+xlab("Founder fitness, %")+ylab("Fitness gain, %")+theme(legend.title=element_blank(),legend.position="none",axis.text=element_text(colour="black"))+geom_abline(slope = -1,linetype="dashed")+scale_color_manual(values=speciesColors)+geom_text(aes(x=c(0,2,-5,-19,-34,-35),y=c(6,2,3,15,22.5,31.4),label=sourceOrganismsOrdered,col=sourceOrganismsOrdered))+geom_errorbar(aes(x=-ancFitness,ymin=evoFitness -evoFitnessSEM, ymax = evoFitness + evoFitnessSEM),width=0)+geom_errorbarh(aes(y=evoFitness,xmin=-ancFitness -ancFitnessSEM, xmax = -ancFitness + ancFitnessSEM),height=0)

ggsave(p1,file=paste(outputDir,"Figure1.svg",sep=""),width=80,height=80,units="mm")
ggsave(p1,file=paste(outputDir,"Figure1.png",sep=""),width=80,height=80,units="mm")


###
mergedData = merge(ancGrowthRateDF,ancFitnessDFTrunc,by="Founder")
founderFactor = factor(mergedData$Founder, levels = sourceOrganismsOrdered)

p1<-ggplot(mergedData)+geom_point(aes(x=-Fitness,y=Avg.Growth.Rate, col = founderFactor),size=3)+xlab("Fitness, %")+ylab("Growth rate, 1/min")+scale_color_manual(values=speciesColors,breaks = sourceOrganismsOrdered)+geom_errorbar(aes(x=-Fitness, ymin = Avg.Growth.Rate - SEMGR, ymax = Avg.Growth.Rate + SEMGR))+geom_errorbarh(aes(y=Avg.Growth.Rate, xmin = -Fitness - stderror, xmax = -Fitness + stderror))+theme_classic()+theme(axis.text=element_text(size=10),axis.title=element_text(size=14),legend.position="none")+geom_text(aes(x=c(-1,1.4,-4,-20,-33.3,-36),y=c(0.0205,0.024,0.02,0.0185,0.0126,0.0128),label=sourceOrganismsOrdered,col=sourceOrganismsOrdered))

ggsave(p1,file=paste(outputDir,"FigureS1.svg",sep=""),width=120,height=120,units="mm")
ggsave(p1,file=paste(outputDir,"FigureS1.png",sep=""),width=120,height=120,units="mm")


###
p1<-ggplot()+geom_crossbar(aes(y = fitnessDFTrunc$Fitness/2, ymin=rep(0,NROW(fitnessDFTrunc)),ymax=fitnessDFTrunc$Fitness,x=c(1:10,12:21,34:43,23:32,56:65,45:54),fill=fitnessFounderListFactor,col=fitnessFounderListFactor),width=.5)+scale_fill_manual(values=speciesColors,breaks = sourceOrganismsOrdered)+scale_color_manual(values=speciesColors, breaks = sourceOrganismsOrdered)+geom_errorbar(aes(x=c(1:10,12:21,34:43,23:32,56:65,45:54), ymin = fitnessDFTrunc$Fitness - fitnessDFTrunc$stderror, ymax = fitnessDFTrunc$Fitness + fitnessDFTrunc$stderror))+geom_hline(aes(yintercept=0),size=.5)+scale_x_continuous(breaks=c(5.5,16.5,38.5,27.5,60.5,49.5),labels=sourceOrganisms)+theme_bw()+xlab("")+ylab("Fitness, %")+theme(axis.text=element_text(size=10),axis.title=element_text(size=14),axis.ticks.x=element_blank(),legend.position="none", panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())

ggsave(p1,file=paste(outputDir,"FigureS2.svg",sep=""),width=80,height=80, units="mm")
ggsave(p1,file=paste(outputDir,"FigureS2.png",sep=""),width=80,height=80, units="mm")



###############################
#
# Figure 2, S3
#
###############################

figure2PopsToUse = c("E_5","S_5","Y_5","V_6","A_5","P_3")

allMutFile<-read.table(allMutFileName,sep="\t",header=FALSE)

selectedMuts = read.table(selectedMutFileName,sep="\t",header=FALSE)
selectedMutsLowThresh = read.table(lowSelectedMutFileName,sep="\t",header=FALSE)

#putatively adaptive mutations are those selected mutations that target a gene that was hit multiple times

mutHitCounts = plyr::count(selectedMuts$V1)
selectedMutsMultipleHits = selectedMuts[!is.na(match(selectedMuts$V1,mutHitCounts$x[mutHitCounts$freq>1])),]
selectedMutsMultipleHits = selectedMutsMultipleHits[selectedMutsMultipleHits$V1 != "",]

mutHitCounts = plyr::count(selectedMutsLowThresh$V1)
selectedMutsLowMultipleHits = selectedMutsLowThresh[!is.na(match(selectedMutsLowThresh$V1,mutHitCounts$x[mutHitCounts$freq>1])),]
selectedMutsLowMultipleHits = selectedMutsLowMultipleHits[selectedMutsLowMultipleHits$V1 != "",]



trajectoryPlotTableFunction <- function(selMuts){
  myTab = selMuts[,c(1,2,11:22)]
  meltedTable = melt(myTab,id.vars=c("V1","V2","V11"))
  meltedTable <- na.omit(meltedTable)
  meltedTable$variable = (c(0:10)*100)[match(meltedTable$variable, c("V12","V13","V14","V15","V16","V17","V18","V19","V20","V21","V22"))]
  meltedTable$locusClass = "Generic"
  meltedTable$locusClass[meltedTable$V1 %in% translationGenes] = "TM-Specific"
  meltedTable$popTitle = as.character(meltedTable$V11)
  meltedTable$popTitle =gsub("_","",meltedTable$popTitle)
  meltedTable$founder = factor(gsub("_.*","",meltedTable$V11),levels=sourceOrganismsOrdered)
  meltedTable$popNumber = gsub("^[A-Z]_","",meltedTable$V11)
  meltedTable$popTitle = factor(meltedTable$popTitle, levels = meltedTable$popTitle[!duplicated(meltedTable$V11)][ order(meltedTable$popNumber[!duplicated(meltedTable$V11)], meltedTable$founder[!duplicated(meltedTable$V11)])])
  
  return(meltedTable)
}

 #divide the plot into 10 positions between 0 and 100% and order gene labels to be close to their actual trajectories without having labels overlap each other
trajectoryPlotLabelFunction <- function(a){
  a$pasted = paste(a$V1,a$V11)
  labelDB = data.frame(popTitle = factor(levels=levels(a$popTitle)), timepoint = numeric(), freq = numeric(), gene = character())
  for(x in unique(a$pasted)){
    myRow = which(a$pasted == x & a$variable == max(a$variable[a$pasted==x]))
    myFreq = a$value[myRow]
    if(myFreq > 0.5){
      myFreq = myFreq - 0.1
    }else{
      myFreq = myFreq + 0.1
    }
    newData = data.frame(popTitle = a$popTitle[myRow], timepoint = a$variable[myRow], freq = myFreq, gene = a$V1[myRow])
    labelDB<-rbind(labelDB,newData)
  }
  
  validPositions = seq(.1,.9,.1)
  labelDB$effectivePosition = NA
  for(pop in unique(labelDB$popTitle)){
    myRows = which(labelDB$popTitle == pop)
    myRows = myRows[order(labelDB$freq[labelDB$popTitle==pop])]
    myPositionArray = rep(NA,9)
    for(row in myRows){
      idealPosition = floor(labelDB$freq[row]*10)
      if(idealPosition == 0){
        idealPosition = 1
      }
      if(is.na(myPositionArray[idealPosition])){
        myPositionArray[idealPosition] = row
        labelDB$effectivePosition[row] = validPositions[idealPosition]
      }else{
        for(counter in c(1:9)){
          if(c(idealPosition+counter) %in% c(1:NROW(myPositionArray))){
            if(is.na(myPositionArray[idealPosition+counter])){
              myPositionArray[idealPosition+counter] = row
              labelDB$effectivePosition[row] = validPositions[idealPosition+counter]
              break
            }
          }
          if(c(idealPosition-counter) %in% c(1:NROW(myPositionArray))){
            if(is.na(myPositionArray[idealPosition-counter])){
              myPositionArray[idealPosition-counter] = row
              labelDB$effectivePosition[row] = validPositions[idealPosition-counter]
              break
            }}
        }
      }
    }
    
  }
  return(labelDB)
  
}

meltedTable = trajectoryPlotTableFunction(selectedMuts)
meltedTable = meltedTable[meltedTable$V11 %in% figure2PopsToUse,]
meltedTable$popTitle = factor(meltedTable$popTitle, levels = meltedTable$popTitle[!duplicated(meltedTable$V11)][ order(meltedTable$founder[!duplicated(meltedTable$V11)], meltedTable$popNumber[!duplicated(meltedTable$V11)])])


labelDB = trajectoryPlotLabelFunction(meltedTable[meltedTable$V1 %in% translationGenes,])

locusColors <- c("darkgrey",brewer.pal(4,"Dark2")[c(2,1)])
names(locusColors) <- c("Generic","TM-Specific","XXX")

chrAmpDFLocal = chrAmpDF[chrAmpDF$popTitle %in% gsub("_","",figure2PopsToUse),]
p1<- ggplot(meltedTable)+
  geom_rect(data = chrAmpDFLocal, fill=brewer.pal(4,"Dark2")[2], col = brewer.pal(4,"Dark2")[2],alpha=.4,linetype=0,aes(xmin = as.numeric(as.character(xmin)), xmax = as.numeric(as.character(xmax)), ymin = 0, ymax = 100))+
  geom_line(data =meltedTable[meltedTable$locusClass == "Generic",],aes(x=variable,y=value*100,group=V2,col=locusClass,size=.6 - 0.3*as.integer(! (V1 %in% selectedMutsMultipleHits$V1))))+
  geom_line(data =meltedTable[meltedTable$locusClass != "Generic",],aes(x=variable,y=value*100,group=V2,col=locusClass,size=.6 - 0.3*as.integer(! (V1 %in% selectedMutsMultipleHits$V1))))+
  scale_color_manual(values=locusColors)+theme_classic()+
  theme(legend.position="none",strip.text=element_blank(),strip.background=element_blank(),panel.border = element_rect(fill=NA,colour = "black"), axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=14,colour="black"))+
  xlab("Generation")+ylab("Frequency, %")+
  ylim(c(0,100))+xlim(c(0,1000))+
  geom_text(data=labelDB,aes(x=timepoint-100,y=effectivePosition*100,group=gene,label=gene,size=5),col="black")+
  geom_text(aes(x=30,y=90,label=popTitle,linetype=NA),col="black",size=6)+
  facet_wrap(popTitle~.,ncol=3)+
  geom_hline(aes(yintercept=100),col="grey")+geom_hline(aes(yintercept=0),col="grey")

ggsave(p1,file=paste(outputDir,"Figure2.svg",sep=""),width=180,height=90, units="mm")
ggsave(p1,file=paste(outputDir,"Figure2.png",sep=""),width=180,height=90, units="mm")



meltedTable = trajectoryPlotTableFunction(selectedMuts)
labelDB = trajectoryPlotLabelFunction(meltedTable[meltedTable$V1 %in% translationGenes,])

p1<- ggplot(meltedTable)+
  geom_rect(data = chrAmpDF, fill=brewer.pal(4,"Dark2")[2], col = brewer.pal(4,"Dark2")[2],alpha=.4,linetype=0,aes(xmin = as.numeric(as.character(xmin)), xmax = as.numeric(as.character(xmax)), ymin = 0, ymax = 100))+
  geom_line(data =meltedTable[meltedTable$locusClass == "Generic",],aes(x=variable,y=value*100,group=V2,col=locusClass,size=.6 - 0.3*as.integer(! (V1 %in% selectedMutsMultipleHits$V1))))+
  geom_line(data =meltedTable[meltedTable$locusClass != "Generic",],aes(x=variable,y=value*100,group=V2,col=locusClass,size=.6 - 0.3*as.integer(! (V1 %in% selectedMutsMultipleHits$V1))))+
  scale_color_manual(values=locusColors)+theme_classic()+
  theme(legend.position="none",strip.text=element_blank(),strip.background=element_blank(),panel.border = element_rect(fill=NA,colour = "black"), axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=14,colour="black"))+
  xlab("Generation")+ylab("Frequency, %")+
  ylim(c(0,100))+xlim(c(0,1000))+
  geom_text(data=labelDB,aes(x=timepoint-100,y=effectivePosition*100,group=gene,label=gene,size=5),col="black")+
  geom_text(aes(x=30,y=90,label=popTitle,linetype=NA),col="black",size=6)+
  facet_wrap(popTitle~.,ncol=6)+
  geom_hline(aes(yintercept=100),col="grey")+geom_hline(aes(yintercept=0),col="grey")
ggsave(p1,file=paste(outputDir,"FigureS4.svg",sep=""),width=360,height=270, units="mm")
ggsave(p1,file=paste(outputDir,"FigureS4.png",sep=""),width=360,height=270, units="mm")

meltedTable = trajectoryPlotTableFunction(selectedMutsLowThresh)
labelDB = trajectoryPlotLabelFunction(meltedTable[meltedTable$V1 %in% translationGenes,])
p1<- ggplot(meltedTable)+
  geom_rect(data = chrAmpDF, fill=brewer.pal(4,"Dark2")[2], col = brewer.pal(4,"Dark2")[2],alpha=.4,linetype=0,aes(xmin = as.numeric(as.character(xmin)), xmax = as.numeric(as.character(xmax)), ymin = 0, ymax = 100))+
  geom_line(data =meltedTable[meltedTable$locusClass == "Generic",],aes(x=variable,y=value*100,group=V2,col=locusClass,size=.6 - 0.3*as.integer(! (V1 %in% selectedMutsMultipleHits$V1))))+
  geom_line(data =meltedTable[meltedTable$locusClass != "Generic",],aes(x=variable,y=value*100,group=V2,col=locusClass,size=.6 - 0.3*as.integer(! (V1 %in% selectedMutsMultipleHits$V1))))+
  scale_color_manual(values=locusColors)+theme_classic()+
  theme(legend.position="none",strip.text=element_blank(),strip.background=element_blank(),panel.border = element_rect(fill=NA,colour = "black"), axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=14,colour="black"))+
  xlab("Generation")+ylab("Frequency, %")+
  ylim(c(0,100))+xlim(c(0,1000))+
  geom_text(data=labelDB,aes(x=timepoint-100,y=effectivePosition*100,group=gene,label=gene,size=5),col="black")+
  geom_text(aes(x=30,y=90,label=popTitle,linetype=NA),col="black",size=6)+
  facet_wrap(popTitle~.,ncol=6)+
  geom_hline(aes(yintercept=100),col="grey")+geom_hline(aes(yintercept=0),col="grey")
ggsave(p1,file=paste(outputDir,"FigureS5.svg",sep=""),width=360,height=270, units="mm")
ggsave(p1,file=paste(outputDir,"FigureS5.png",sep=""),width=360,height=270, units="mm")



meltedTable = trajectoryPlotTableFunction(selectedMuts)
meltedTable$locusClass[meltedTable$V1 %in% cytokinesisGenes] = "Cytokinesis"
labelDB = trajectoryPlotLabelFunction(meltedTable[meltedTable$V1 %in% c(cytokinesisGenes),])
myCols = locusColors
names(myCols) <- c("Generic","TM-Specific","Cytokinesis")
p1<- ggplot(meltedTable)+
  geom_rect(data = chrAmpDF, fill=brewer.pal(4,"Dark2")[2], col = brewer.pal(4,"Dark2")[2],alpha=.4,linetype=0,aes(xmin = as.numeric(as.character(xmin)), xmax = as.numeric(as.character(xmax)), ymin = 0, ymax = 100))+
  geom_line(data =meltedTable[meltedTable$locusClass == "Generic",],aes(x=variable,y=value*100,group=V2,col=locusClass,size=.6 - 0.3*as.integer(! (V1 %in% selectedMutsMultipleHits$V1))))+
  geom_line(data =meltedTable[meltedTable$locusClass != "Generic",],aes(x=variable,y=value*100,group=V2,col=locusClass,size=.6 - 0.3*as.integer(! (V1 %in% selectedMutsMultipleHits$V1))))+
  scale_color_manual(values=myCols)+theme_classic()+
  theme(legend.position="none",strip.text=element_blank(),strip.background=element_blank(),panel.border = element_rect(fill=NA,colour = "black"), axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=14,colour="black"))+
  xlab("Generation")+ylab("Frequency, %")+
  ylim(c(0,100))+xlim(c(0,1000))+
  geom_text(data=labelDB,aes(x=timepoint-100,y=effectivePosition*100,group=gene,label=gene,size=5),col="black")+
  geom_text(aes(x=30,y=90,label=popTitle,linetype=NA),col="black",size=6)+
  facet_wrap(popTitle~.,ncol=6)+
  geom_hline(aes(yintercept=100),col="grey")+geom_hline(aes(yintercept=0),col="grey")
ggsave(p1,file=paste(outputDir,"FigureS6.svg",sep=""),width=360,height=270, units="mm")
ggsave(p1,file=paste(outputDir,"FigureS6.png",sep=""),width=360,height=270, units="mm")


###############################
#
# Figure 3
#
###############################


allMutFile<-read.table(allMutFileName,sep="\t",header=FALSE)

selectedMuts = read.table(selectedMutFileName,sep="\t",header=FALSE)
selectedMutsLowThresh = read.table(lowSelectedMutFileName,sep="\t",header=FALSE)
selectedMutsUltraLowThresh = read.table(ultraLowSelectedMutFileName,sep="\t",header=FALSE)

#putatively adaptive mutations are those selected mutations that target a gene that was hit multiple times

mutHitCounts = plyr::count(selectedMuts$V1)
selectedMutsMultipleHits = selectedMuts[!is.na(match(selectedMuts$V1,mutHitCounts$x[mutHitCounts$freq>1])),]
selectedMutsMultipleHits = selectedMutsMultipleHits[selectedMutsMultipleHits$V1 != "",]

mutHitCounts = plyr::count(selectedMutsLowThresh$V1)
selectedMutsLowMultipleHits = selectedMutsLowThresh[!is.na(match(selectedMutsLowThresh$V1,mutHitCounts$x[mutHitCounts$freq>1])),]
selectedMutsLowMultipleHits = selectedMutsLowMultipleHits[selectedMutsLowMultipleHits$V1 != "",]

mutHitCounts = plyr::count(selectedMutsUltraLowThresh$V1)
selectedMutsUltraLowMultipleHits = selectedMutsUltraLowThresh[!is.na(match(selectedMutsUltraLowThresh$V1,mutHitCounts$x[mutHitCounts$freq>1])),]
selectedMutsUltraLowMultipleHits = selectedMutsUltraLowMultipleHits[selectedMutsUltraLowMultipleHits$V1 != "",]



Figure3Function <- function(mutSet, geneList, listName, plotName, speciesSet){
  mutIdentityDF = data.frame(founder = character(), gene = character(), mutationClass = character(), locusClass = character())
  
  
  for(gene in unique(as.character(mutSet$V1))){
  	mutClass = ""
  	if(NROW(unique(populationTypes[match(gsub("_.","",mutSet$V11[mutSet$V1 == gene]),sourceOrganismsOrdered)]))==2){
  		mutClass = "Generic"
  	}else if(NROW(unique(gsub("_.","",mutSet$V11[mutSet$V1 == gene])))==1){
  		mutClass = "Founder specific"
  	}else{
  		mutClass = "Fitness specific"
  	}
  	locusClass = "Generic"
  	if(gene %in% geneList){
  		locusClass = listName
  	}
  	
  	mutIdentityDF = rbind(mutIdentityDF, data.frame(founder = gsub("_.","",mutSet$V11[mutSet$V1 == gene]),gene = gene, mutationClass = mutClass, locusClass = locusClass))
  }
  
  tufAType = "Generic"
  if("tufA" %in% geneList){
    tufAType = listName
  }
  
  mutIdentityDF = rbind(mutIdentityDF, data.frame(founder = gsub("_.","",chrAmptufAPopulations),gene = "tufA", mutationClass = "Fitness specific", locusClass = tufAType))
  
  mutCountDF = data.frame(founder = character(), locusClass = character(), count = numeric())
  for(founder in sourceOrganismsOrdered){
    for(locusClass in unique(mutIdentityDF$locusClass)){
      myCount = sum(mutIdentityDF$founder == founder & mutIdentityDF$locusClass == locusClass,na.rm=TRUE)
      mutCountDF = rbind(mutCountDF,data.frame(founder = founder, locusClass = locusClass, count = myCount))
    }
  }
  mycols = locusColors
  names(mycols) <- c("Generic","TM-Specific",listName)
  if(listName == "TM-Specific"){
    names(mycols) <- c("Generic","TM-Specific","XXX")
  }
  mutCountDF$founder <- factor(mutCountDF$founder, levels = sourceOrganismsOrdered)
 
  
  
  
  myfun <- function(x){min(which(x),na.rm=TRUE)}
  getCumDist <- function(mutDF, freqThreshold){
    outputDF = data.frame(freqThresh = numeric(), timepoint = numeric(), numMuts = numeric())
    
    firstGenOverThreshAll = (as.numeric(apply(mutDF[,12:22]>freqThreshold,1, myfun))-1)*100
    myCounts = plyr::count(firstGenOverThreshAll[!is.na(firstGenOverThreshAll)])
    validTimes = c(0:10)*100
    numCounts = rep(0,11)
    for(i in myCounts$x){
      numCounts[validTimes == i] = myCounts$freq[myCounts$x == i]
    }
    outputDF = rbind(outputDF, data.frame(freqThresh = freqThreshold, timepoint = validTimes, numMuts = cumsum(numCounts)))
    
    return(outputDF)
  }
  
  outputDF = getCumDist(mutSet[mutSet$V1 %in% geneList & gsub("_.","",mutSet$V11) %in% speciesSet,],0.95)
  outputDF$Type = listName
  outputDF2 = getCumDist(mutSet[! mutSet$V1 %in% geneList & gsub("_.","",mutSet$V11) %in% speciesSet,],0.95)
  outputDF2$Type = "Generic"
  merged = rbind(outputDF,outputDF2)
  
  p1<-ggplot(mutCountDF)+
    geom_bar(aes(x=founder,y=count/6,group=paste(founder,locusClass),fill=locusClass),stat="identity",position="stack")+
    theme_bw()+
    scale_fill_manual(values = mycols)+
    xlab("Founder")+ylab("Detected mutations per population")+
    theme(legend.title=element_blank(),legend.text=element_text(size=10),legend.position="bottom",axis.text=element_text(size=10),axis.title=element_text(size=12),panel.grid.major.x = element_blank())+
    annotate(geom = 'text', label = 'A', x = -Inf, y = Inf, hjust = -.5, vjust = 1.2,size=6,fontface="bold")
  
  p2<-ggplot(merged)+geom_line(aes(x=timepoint,y=numMuts/18,group=Type,col=Type),size=2)+
    scale_color_manual(values=mycols)+
    theme_bw()+
    xlab("Generation")+ylab("Fixed mutations per population")+
    theme(legend.title=element_blank(),legend.position="bottom",axis.text=element_text(size=10),axis.title=element_text(size=12),legend.text=element_text(size=10))+
    annotate(geom = 'text', label = 'B', x = -Inf, y = Inf, hjust = -.5, vjust = 1.2,size=6,fontface="bold")
  
  p3 = arrangeGrob(p1, p2, nrow = 1)
  ggsave(p3,file=paste(outputDir,plotName,".svg",sep=""),width=160,height=100, units="mm")
  ggsave(p3,file=paste(outputDir,plotName,".png",sep=""),width=160,height=100, units="mm")
}

Figure3Function(selectedMutsMultipleHits, translationGenes, "TM-Specific","Figure3",c("V","P","A"))
Figure3Function(selectedMutsMultipleHits, cytokinesisGenes, "Cytokinesis","FigureS7",sourceOrganisms)

reconstructedFitnessDF = read.table(paste(outputDir,"reconstructedFitnessData.tab",sep=""),sep="\t",header=TRUE)
reconstructedFitnessDF$Fitness = reconstructedFitnessDF$Fitness*100
reconstructedFitnessDF$stderror = reconstructedFitnessDF$stderror*100
reconstructedFitnessDF$variance = reconstructedFitnessDF$variance*100^2
reconstructedFitnessDFTrunc = reconstructedFitnessDF[grep("allFlasks_averaged", reconstructedFitnessDF$Name),]


###############################
#
# Figure 4
#
###############################

## heatmap of selected muts with multiple hits

selectedMuts = read.table(selectedMutFileName,sep="\t",header=FALSE)

mutHitCounts = plyr::count(selectedMuts$V1)
selectedMutsMultipleHits = selectedMuts[!is.na(match(selectedMuts$V1,mutHitCounts$x[mutHitCounts$freq>1])),]
selectedMutsMultipleHits = selectedMutsMultipleHits[selectedMutsMultipleHits$V1 != "",]

mutIdentityDF = data.frame(founder = character(), gene = character(), mutationClass = character(), locusClass = character())


for(gene in unique(as.character(selectedMutsMultipleHits$V1))){
  mutClass = ""
  if(NROW(unique(populationTypes[match(gsub("_.","",selectedMutsMultipleHits$V11[selectedMutsMultipleHits$V1 == gene]),sourceOrganismsOrdered)]))==2){
    mutClass = "Generic"
  }else if(NROW(unique(gsub("_.","",selectedMutsMultipleHits$V11[selectedMutsMultipleHits$V1 == gene])))==1){
    mutClass = "Founder specific"
  }else{
    mutClass = "Fitness specific"
  }
  locusClass = "Generic"
  if(gene %in% translationGenes){
    locusClass = "TM-Specific"
  }
  
  mutIdentityDF = rbind(mutIdentityDF, data.frame(founder = gsub("_.","",selectedMutsMultipleHits$V11[selectedMutsMultipleHits$V1 == gene]),gene = gene, mutationClass = mutClass, locusClass = locusClass))
}
mutIdentityDF = rbind(mutIdentityDF, data.frame(founder = gsub("_.","",chrAmptufAPopulations),gene = "tufA amp", mutationClass = "Fitness specific", locusClass = "TM-Specific"))

mutCountDF = data.frame(founder = character(), locusClass = character(), count = numeric())
for(founder in sourceOrganismsOrdered){
  for(locusClass in unique(mutIdentityDF$locusClass)){
    myCount = sum(mutIdentityDF$founder == founder & mutIdentityDF$locusClass == locusClass)
    mutCountDF = rbind(mutCountDF,data.frame(founder = founder, locusClass = locusClass, count = myCount))
  }
}



#randomization using only TM genes
myDF2 = mutIdentityDF[ mutIdentityDF$gene %in% c(translationGenes,"tufA amp"),]
myDF2$founder = factor(myDF2$founder,levels=sourceOrganismsOrdered)

castedDF2 = dcast(myDF2, gene ~ founder)
secondSort = rowSums(castedDF2[2:NCOL(castedDF2)]>0)
thirdSort = apply(castedDF2,1,function(x){sum((as.numeric(x[2:NCOL(castedDF2)])>0) * c(64,32,16,8,4,2))})

multimutCountsByFounder = as.numeric(colSums(castedDF2[,2:NCOL(castedDF2)]))
numRandomizations = 10000
geneEntropyDFTM = data.frame(gene = character(), pvalue = numeric())
for (gene in castedDF2$gene[order(secondSort,thirdSort)]){
  myRow = as.numeric(castedDF2[castedDF2$gene == gene,2:NCOL(castedDF2)])
  totalNumHits= sum(myRow)
  observedFreqs = myRow / totalNumHits
  observedEntropy = -sum(observedFreqs[observedFreqs>0] * log(observedFreqs[observedFreqs>0]))
  randomEntropyDist = c()
  for(i in c(1:numRandomizations)){
    randomCounts = as.numeric(rmultinom(1,totalNumHits,multimutCountsByFounder))
    randomFreqs = randomCounts / totalNumHits
    randomEntropy = -sum(randomFreqs[randomFreqs>0] * log(randomFreqs[randomFreqs>0]))
    randomEntropyDist = c(randomEntropyDist,randomEntropy)
  }
  geneEntropyDFTM = rbind(geneEntropyDFTM,data.frame(gene = gene, pvalue = as.character(sum(randomEntropyDist < observedEntropy)/numRandomizations)))
}
geneEntropyDFTM$bonP = p.adjust(as.numeric(as.character(geneEntropyDFTM$pvalue)),method="bonferroni")
geneEntropyDFTM$bhP = p.adjust(as.numeric(as.character(geneEntropyDFTM$pvalue)),method="BH")


totalGeneCounts = as.character(rowSums(castedDF2[,2:NCOL(castedDF2)]))
origGeneList = as.character(castedDF2$gene)
newGeneList = origGeneList

pvalueThreshold = 0.05
newGeneList = paste(newGeneList," (",as.character(totalGeneCounts),")",sep="")

#magical reordering to get everything in the right order for the plot
meltedDF2TM = melt(castedDF2,id.var=c("gene"))
geneLevels = as.character(origGeneList[order(secondSort,thirdSort, decreasing=TRUE)])
reordering = match(geneLevels,origGeneList) #origGeneList[reordering] == geneLevels

meltedDF2TM$gene <- as.character(meltedDF2TM$gene)
meltedDF2TM$gene = newGeneList[match(meltedDF2TM$gene,origGeneList)]
meltedDF2TM$gene = factor(meltedDF2TM$gene,levels=rev(newGeneList[reordering]))

geneType = rep("Generic",NROW(castedDF2))
geneType[castedDF2$gene[order(secondSort,thirdSort, decreasing=TRUE)] %in% c(translationGenes,"tufA amp")] = "TM-Specific"

meltedDF2TM$geneType = geneType[match(meltedDF2TM$gene,newGeneList[reordering])]
meltedDF2TM$variable <- factor(as.character(meltedDF2TM$variable),rev(sourceOrganismsOrdered))

geneFaceTM = rep("plain",NROW(geneLevels))
geneFaceTM[geneLevels %in% as.character(geneEntropyDFTM$gene[as.numeric(as.character(geneEntropyDFTM$bhP)) < pvalueThreshold])]= "bold"
meltedDF2TM = rbind(meltedDF2TM, data.frame(gene = rep(unique(meltedDF2TM$gene),3), variable = factor(rep(c("E","S","Y"),each=NROW(unique(meltedDF2TM$gene))),rev(sourceOrganismsOrdered)),value = 0, geneType = "TM-Specific" ))
p1<-ggplot(meltedDF2TM)+geom_tile(aes(gene, variable,fill=value),color="white") +
  scale_fill_gradientn(colors = c("white","grey20"),breaks=c(7:1),guide="legend",limits = c(0,7))+theme_bw()+
  ggtitle("TM-Specific")+ylab("Founder")+
  theme(axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=10),axis.text=element_text(size=10), axis.text.x=element_text(angle=60,hjust=1,face=rev(geneFaceTM)),plot.title=element_text(size=14,hjust = .5), legend.position="none")+
  geom_text(aes(x=2,y=6,label="A"),size=10,fontface="bold")


#randomization using only Generic genes

myDF2 = mutIdentityDF[! mutIdentityDF$gene %in% c(translationGenes,"tufA amp"),]
myDF2$founder = factor(myDF2$founder,levels=sourceOrganismsOrdered)

castedDF2 = dcast(myDF2, gene ~ founder)
secondSort = rowSums(castedDF2[2:NCOL(castedDF2)]>0)
thirdSort = apply(castedDF2,1,function(x){sum((as.numeric(x[2:NCOL(castedDF2)])>0) * c(64,32,16,8,4,2))})


multimutCountsByFounder = as.numeric(colSums(castedDF2[,2:NCOL(castedDF2)]))
numRandomizations = 10000
geneEntropyDFGeneric = data.frame(gene = character(), pvalue = numeric())
for (gene in castedDF2$gene[order(secondSort,thirdSort)]){
  myRow = as.numeric(castedDF2[castedDF2$gene == gene,2:NCOL(castedDF2)])
  totalNumHits= sum(myRow)
  observedFreqs = myRow / totalNumHits
  observedEntropy = -sum(observedFreqs[observedFreqs>0] * log(observedFreqs[observedFreqs>0]))
  randomEntropyDist = c()
  for(i in c(1:numRandomizations)){
    randomCounts = as.numeric(rmultinom(1,totalNumHits,multimutCountsByFounder))
    randomFreqs = randomCounts / totalNumHits
    randomEntropy = -sum(randomFreqs[randomFreqs>0] * log(randomFreqs[randomFreqs>0]))
    randomEntropyDist = c(randomEntropyDist,randomEntropy)
  }
  geneEntropyDFGeneric = rbind(geneEntropyDFGeneric,data.frame(gene = gene, pvalue = as.character(sum(randomEntropyDist < observedEntropy)/numRandomizations)))
}
geneEntropyDFGeneric$bonP = p.adjust(as.numeric(as.character(geneEntropyDFGeneric$pvalue)),method="bonferroni")
geneEntropyDFGeneric$bhP = p.adjust(as.numeric(as.character(geneEntropyDFGeneric$pvalue)),method="BH")

totalGeneCounts = as.character(rowSums(castedDF2[,2:NCOL(castedDF2)]))
origGeneList = as.character(castedDF2$gene)
newGeneList = origGeneList

pvalueThreshold = 0.05
newGeneList = paste(newGeneList," (",as.character(totalGeneCounts),")",sep="")

#magical reordering to get everything in the right order for the plot
meltedDF2Generic = melt(castedDF2,id.var=c("gene"))
geneLevels = as.character(origGeneList[order(secondSort,thirdSort, decreasing=TRUE)])
reordering = match(geneLevels,origGeneList) #origGeneList[reordering] == geneLevels

meltedDF2Generic$gene <- as.character(meltedDF2Generic$gene)
meltedDF2Generic$gene = newGeneList[match(meltedDF2Generic$gene,origGeneList)]
meltedDF2Generic$gene = factor(meltedDF2Generic$gene,levels=rev(newGeneList[reordering]))

geneType = rep("Generic",NROW(castedDF2))
geneType[castedDF2$gene[order(secondSort,thirdSort, decreasing=TRUE)] %in% c(translationGenes,"tufA amp")] = "TM-Specific"

meltedDF2Generic$geneType = geneType[match(meltedDF2Generic$gene,newGeneList[reordering])]
meltedDF2Generic$variable <- factor(as.character(meltedDF2Generic$variable),rev(sourceOrganismsOrdered))

geneFaceGeneric = rep("plain",NROW(geneLevels))
geneFaceGeneric[geneLevels %in% as.character(geneEntropyDFGeneric$gene[as.numeric(as.character(geneEntropyDFGeneric$bhP)) < pvalueThreshold])]= "bold"

p2<-ggplot(meltedDF2Generic)+geom_tile(aes(gene, variable,fill=value),color="white") +
  scale_fill_gradientn(colors = c("white","grey20"),breaks=c(7:1),guide="legend",limits = c(0,7))+
  ggtitle("Generic")+xlab("Gene")+labs(fill="Number of\nmutations")+theme_bw()+
  theme(axis.title=element_text(size=12),axis.title.x=element_text(hjust=.25),legend.title=element_text(size=12),legend.text=element_text(size=10),axis.text=element_text(size=10), axis.text.x=element_text(angle=60,hjust=1,face=rev(geneFaceGeneric)),plot.title=element_text(size=14, hjust = .5), axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+geom_text(aes(x=2,y=6,label="B"),size=10,fontface="bold")



p3<- ggarrange(p1,p2,widths=c(.65,2))

ggsave(p3,file=paste(outputDir,"Figure4.svg",sep=""),width=180,height=100, units="mm")
ggsave(p3,file=paste(outputDir,"Figure4.png",sep=""),width=180,height=100, units="mm")


###############################
#
# Figure 5
#
###############################
reconstructedFitnessDF = read.table(paste(outputDir,"reconstructedFitnessData.tab",sep=""),sep="\t",header=TRUE)

reconstructedFitnessDF$Fitness = reconstructedFitnessDF$Fitness*100
reconstructedFitnessDF$stderror = reconstructedFitnessDF$stderror*100
reconstructedFitnessDF$variance = reconstructedFitnessDF$variance*100^2


reconstructedFitnessDFTrunc = reconstructedFitnessDF[grep("allFlasks_averaged", reconstructedFitnessDF$Name),]

reconstructedFitnessDFTrunc2 = reconstructedFitnessDFTrunc[c(3:8),]
myNames = gsub("_allFlasks_averaged","",reconstructedFitnessDFTrunc2$Name)
#myNames = gsub("BK[2-7] ","",myNames)
myNamesFactor = factor(myNames, levels = myNames)
myBKs = gsub(" .*","",reconstructedFitnessDFTrunc2$Name[grep("allFlasks_averaged", reconstructedFitnessDFTrunc2$Name)])
myBKsFactor = factor(myBKs, levels = sourceOrganismsOrdered)
p1<-ggplot(reconstructedFitnessDFTrunc2)+geom_bar(aes(x=myBKsFactor,y=Fitness, fill = myBKsFactor),stat="identity",position="dodge")+xlab("Founder")+ylab("Fitness, %")+scale_fill_manual(values=speciesColors,breaks = sourceOrganismsOrdered)+geom_errorbar(aes(x=myBKsFactor, ymin = Fitness - stderror, ymax = Fitness + stderror),width=.2)+theme_bw()+theme(axis.text=element_text(size=16),axis.title=element_text(size=20),legend.title=element_blank(),legend.position="none")
ggsave(p1,file=paste(outputDir,"Figure5.svg",sep=""),width=140,height=140,units="mm")
ggsave(p1,file=paste(outputDir,"Figure5.png",sep=""),width=140,height=140,units="mm")


###############################
#
# Figure S3 - coverage plots
#
###############################
coveragePlotDF = data.frame(sample = character(), position = numeric(), coverage = numeric())
for(i in sourceOrganismsOrdered){
  for(j in c(1:6)){
    for(k in c(0:10)*100){
      
      filename = paste(dataInputDir,"coverage_data/",i,"_",j,"_",k,"_averageGenomeCov.tab",sep="")
      #if(file.exists(filename) & !file.exists(paste(outputDir,"coverage_plots/",i,"_",j,"_",k,"_covPlot.png",sep=""))){
      if(file.exists(filename)){
        myDF = read.table(filename, sep="\t")
        myDF=myDF[,c(2:3)]
        names(myDF) <- c("position","coverage")
        myDF$sample = paste(i,j,", t",k,sep="")
        coveragePlotDF = rbind(coveragePlotDF,myDF)
        #p1<-ggplot(myDF)+geom_point(aes(x=x,y=y)) + theme_bw()+xlab("Genomic Coordinate")+ylab("1kb average coverage")+ geom_vline(xintercept=c(3420000,3520000), color = "red", linetype="longdash")+ggtitle(paste(i,'-',j,', t',k,sep=""))+theme(axis.title=element_text(size=16),axis.text=element_text(size=12),plot.title=element_text(size=20))
        #ggsave(p1,file=paste(outputDir,"coverage_plots/",i,"_",j,"_",k,"_covPlot.png",sep=""),width=120,height=80,units="mm")
      }
    }
  }
}

#plots to include
plotsToUse = c("V1, t100", "V3, t700", "V5, t600", "A1, t100", "A2, t600", "A3, t100", "A4, t700", "A5, t100", "A6, t100", "P4, t100", "P6, t200")
coveragePlotDFTrunc = coveragePlotDF[coveragePlotDF$sample %in% plotsToUse,]

p1<-ggplot(coveragePlotDFTrunc)+geom_point(aes(x=position,y=coverage)) + theme_bw()+xlab("Genomic Coordinate")+ylab("1kb average coverage")+ geom_vline(xintercept=c(3420000,3520000), color = "red", linetype="longdash")+theme(axis.title=element_text(size=16),axis.text=element_text(size=12,colour="black"),plot.title=element_text(size=20),strip.text=element_text(size=20),strip.background=element_blank(),panel.border = element_rect(fill=NA,colour = "black"))+facet_wrap(sample~.,ncol=3,scales="free_y")

ggsave(p1,file=paste(outputDir,"FigureS3.png",sep=""),width=120*3,height=80*4,units="mm")
ggsave(p1,file=paste(outputDir,"FigureS3.svg",sep=""),width=120*3,height=80*4,units="mm")

###############################
#
# Figure S8
#
###############################
threeWayCompetitionData = read.table(paste(dataInputDir,"ThreeWayCompetition_rawFreqData.tab",sep=""),sep="\t",header=TRUE)
threeWayCompetitionData<-threeWayCompetitionData[!is.na(threeWayCompetitionData$Frequency),]
threeWayCompetitionData$Title = paste(threeWayCompetitionData$Strain,threeWayCompetitionData$Replicate,sep="-")
threeWayCompetitionData$Title<-factor(threeWayCompetitionData$Title,levels = paste(rep(sourceOrganismsOrdered,3),rep(c(1:3),each=6),sep="-"))

p1<-ggplot(threeWayCompetitionData,aes(x=t,y=Frequency,group=Marker,col=Marker))+geom_line(size=2)+xlab("Time (days)")+ylab("Frequency")+facet_wrap(Title~.,nrow=3)+theme_bw()+scale_color_manual(values=replicateColors)+ theme(legend.position="bottom",strip.text=element_text(size=13),strip.background=element_blank(),axis.text=element_text(size=10),axis.title=element_text(size=14))+scale_x_continuous(breaks=c(1,2,3))

ggsave(p1,file=paste(outputDir,"FigureS8.svg",sep=""),width=180,height=135,units="mm")
ggsave(p1,file=paste(outputDir,"FigureS8.png",sep=""),width=180,height=135,units="mm")



###############################
#
# Table S1
#
###############################

selectedMuts = read.table(selectedMutFileName,sep="\t",header=FALSE)

TableS1FileName <- paste(outputDir,"TableS1.tab",sep="")
selectedMuts$multipleHits = !is.na(match(selectedMuts$V1,mutHitCounts$x[mutHitCounts$freq>1])) & selectedMuts$V1 != ""
selectedMuts<-as.data.frame(selectedMuts)

selectedMuts = bind_rows(selectedMuts,data.frame(V1 = "tufA",V5="ChrAmp",V11=chrAmptufAPopulations,multipleHits=TRUE))  

names(selectedMuts)<-c("Gene","Mut ID","MG1655 Position","DNA Mutation","Mutation Type","CDS position","AA position","AA mutation","Upstream Gene","Downstream Gene","Population ID","Gen 0 frequency","Gen 100 frequency","Gen 200 frequency","Gen 300 frequency","Gen 400 frequency","Gen 500 frequency","Gen 600 frequency","Gen 700 frequency","Gen 800 frequency","Gen 900 frequency","Gen 1000 frequency","Gen 0 read depth","Gen 100 read depth","Gen 200 read depth","Gen 300 read depth","Gen 400 read depth","Gen 500 read depth","Gen 600 read depth","Gen 700 read depth","Gen 800 read depth","Gen 900 read depth","Gen 1000 read depth","Multiply hit locus")
  
selectedMuts2 = apply(selectedMuts,2,function(x){x[is.na(x)]=""
return(x)})

write.table(selectedMuts2, file=TableS1FileName,sep="\t",col.names=TRUE,quote=FALSE,row.names=FALSE)



###############################
#
# Numbers and Stats
#
###############################

### Stats
print("Founder fitness")
print(ancFitnessDFTrunc)

fitnessFounderList = gsub("_.*","",fitnessDFTrunc$Name)
fitnessFounderListFactor = factor(fitnessFounderList, levels = sourceOrganismsOrdered)
ancFitnessFounderList = gsub(" .*","",ancFitnessDFTrunc$Name)
ancFitnessFounderListFactor = factor(ancFitnessFounderList, levels = sourceOrganismsOrdered)
ancFitnessDFTrunc$Founder = ancFitnessFounderList
x = ancFitnessDFTrunc$Fitness[match(fitnessFounderList,ancFitnessFounderList)]
y = fitnessDFTrunc$Fitness
fractionFitnessGain = y/x

print("Average % fitness recovery of evolved V, A and P populations relative to E assuming fitness transitivity")
print(mean(fractionFitnessGain[fitnessFounderList %in% c("V","A","P")]))


print("number of evolved populations that did not significantly increase in fitness, B-H correction")
print(sum(p.adjust(sort(pt(-abs(fitnessDFTrunc$Fitness/fitnessDFTrunc$stderror),df=3)),method="BH")>=0.05))

print("Total num adaptive mutations excepting chr amplifications")
print(NROW(selectedMutsMultipleHits))

print("Total num adaptive mutations")
print(NROW(selectedMutsMultipleHits)+NROW(chrAmptufAPopulations))

print("Num genetic targets for adaptation")
print(NROW(unique(c(as.character(selectedMutsMultipleHits$V1),"tufA"))))

#expected number of adaptive mutations and genetic targets for those mutations
geneLengths = read.table(paste(dataInputDir,"MG1655_geneLengths.tab",sep=""),sep="\t",header=FALSE)
multinomResultsGenes = c()
multinomResultsMuts = c()
for(i in c(1:10000)){
  multRes =   rmultinom(1,NROW(selectedMutsMultipleHits),geneLengths$V4)
  multinomResultsGenes = c(multinomResultsGenes, sum(multRes>1))
  multinomResultsMuts = c(multinomResultsMuts, sum(multRes[multRes>1]))
}
print("Expected FDR (%) of adapted mutations")
print(100*mean(multinomResultsMuts)/(NROW(selectedMutsMultipleHits)))

multinomResultsGenes = c()
multinomResultsMuts = c()
for(i in c(1:10000)){
  multRes =   rmultinom(1,NROW(selectedMutsLowMultipleHits),geneLengths$V4)
  multinomResultsGenes = c(multinomResultsGenes, sum(multRes>1))
  multinomResultsMuts = c(multinomResultsMuts, sum(multRes[multRes>1]))
}
print("Expected FDR (%) of adapted mutations w/ low threshold for calling variants")
print(100*mean(multinomResultsMuts)/(NROW(selectedMutsLowMultipleHits)))


print("Number of TM-specific adaptive mutations")
print(sum(selectedMutsMultipleHits$V1 %in% translationGenes) + NROW(chrAmptufAPopulations))

print("Number of TM-specific adapted genes")
print(NROW(unique(c(as.character(selectedMutsMultipleHits$V1[selectedMutsMultipleHits$V1 %in% translationGenes]),"tufA"))))

print("Number of TM annotated genes in the genome")
print(NROW(translationGenes))

multinomResultsGenes = c()
multinomResultsMuts = c()
for(i in c(1:10000)){
  multRes =   rmultinom(1,NROW(selectedMuts),geneLengths$V4)
  multinomResultsGenes = c(multinomResultsGenes, sum(multRes>1 & geneLengths$V1 %in% translationGenes))
  multinomResultsMuts = c(multinomResultsMuts, sum(multRes[multRes>1 & geneLengths$V1 %in% translationGenes]))
}
print("Num randomizations out of 10k with excess observed TM-specific adapted genes than simulated")
print(sum(multinomResultsMuts <= sum(selectedMutsMultipleHits$V1 %in% translationGenes) + NROW(chrAmptufAPopulations)))

print("Num translation genes")
print(NROW(translationGenes))


print("Percentage of E. coli genome that is the CDS of TM annotated genes")
positionValues = rep(0,4641652)
for(gene in translationGenes){
  myRow = geneLengths[geneLengths$V1 == gene,]
  if(NROW(myRow)==1){
    positionValues[c(myRow$V2:myRow$V3)]=1
  }
}
print(sum(positionValues)/4641652)

print("Num populations with TM-specific mutations")
print(NROW(unique(c(as.character(selectedMutsMultipleHits$V11[as.character(selectedMutsMultipleHits$V1) %in% translationGenes]),chrAmptufAPopulations ))))

print("Num fixed TM-specific / generic mutations")

myfun <- function(x){min(which(x))}
selectedMutsMultipleHitsLowFitnessFoundersOnly = selectedMutsMultipleHits[gsub("_.","",selectedMutsMultipleHits$V11) %in% c("P","A","V"),]
firstGenOverThreshAll = (as.numeric(apply(selectedMutsMultipleHitsLowFitnessFoundersOnly[,12:22]>0.95,1, myfun))-1)*100
firstGenOverThreshAll[!is.finite(firstGenOverThreshAll)]=NA

TMGeneFixCounts = data.frame(founder= character(), population = character(), totalTMMutCount=numeric(), TMFixMutCount = numeric(), TMFirstFixMutCount = numeric(), NonTMMutCount = numeric(), NonTMFixMutCount = numeric(), NonTMFirstFixMutCount=numeric())
for(pop in sort(unique(selectedMutsMultipleHitsLowFitnessFoundersOnly$V11))){
  subtable = selectedMutsMultipleHitsLowFitnessFoundersOnly[selectedMutsMultipleHitsLowFitnessFoundersOnly$V11==pop,]
  firstGenThreshVals = firstGenOverThreshAll[selectedMutsMultipleHitsLowFitnessFoundersOnly$V11==pop]
  TMGeneFixCounts = rbind(TMGeneFixCounts, data.frame(founder = gsub("_.","",pop), population = pop, totalTMMutCount = sum(subtable$V1 %in% translationGenes), TMFixMutCount=  sum(subtable$V1 %in% translationGenes & !is.na(firstGenThreshVals)),TMFirstFixMutCount= sum(subtable$V1 %in% translationGenes & firstGenThreshVals == min(firstGenThreshVals,na.rm=TRUE),na.rm=TRUE), NonTMMutCount = sum(!(subtable$V1 %in% translationGenes)), NonTMFixMutCount = sum(!(subtable$V1 %in% translationGenes) & !is.na(firstGenThreshVals)), NonTMFirstFixMutCount=sum(!(subtable$V1 %in% translationGenes) & firstGenThreshVals == min(firstGenThreshVals,na.rm=TRUE),na.rm=TRUE)))
}
print(colSums(TMGeneFixCounts[,3:8])) #total number of TM muts, total num fixed, and total num fixed first

print("Timing of fixation of TM-specific and generic mutations, mean and sem")
print("Generic genes")
print(mean(firstGenOverThreshAll[which(! selectedMutsMultipleHitsLowFitnessFoundersOnly$V1 %in% translationGenes)],na.rm=TRUE)) #Generic genes
print(sd(firstGenOverThreshAll[which(! selectedMutsMultipleHitsLowFitnessFoundersOnly$V1 %in% translationGenes)],na.rm=TRUE) / sqrt(sum(! selectedMutsMultipleHitsLowFitnessFoundersOnly$V1 %in% translationGenes & !is.na(firstGenOverThreshAll)))) #SEM Generic genes
print("TM genes")
print(mean(firstGenOverThreshAll[which(selectedMutsMultipleHitsLowFitnessFoundersOnly$V1 %in% translationGenes)],na.rm=TRUE)) #TM genes
print(sd(firstGenOverThreshAll[which(selectedMutsMultipleHitsLowFitnessFoundersOnly$V1 %in% translationGenes)],na.rm=TRUE) / sqrt(sum( selectedMutsMultipleHitsLowFitnessFoundersOnly$V1 %in% translationGenes & !is.na(firstGenOverThreshAll)))) #SEM TM genes

print("Num TM-specific mutations that fix after gen 600")
print(sum(firstGenOverThreshAll[which(selectedMutsMultipleHitsLowFitnessFoundersOnly$V1 %in% translationGenes)] > 600 , na.rm=TRUE))

print("Num Generic mutations that fix after gen 600")
print(sum(firstGenOverThreshAll[which(! selectedMutsMultipleHitsLowFitnessFoundersOnly$V1 %in% translationGenes)] > 600 , na.rm=TRUE))

print("Mean fitness gain of Y populations")
print(mean(fitnessDFTrunc$Fitness[fitnessFounderList=="Y"]))

print("Avg num TM-specific fixations in V, A and P populations")
print(as.numeric(colSums(TMGeneFixCounts[,3:8])[2])/18)

print("Average fitness deficit (% per generation) of evolved V, A and P populations relative to E assuming fitness transitivity")
print(mean((x-y)[fitnessFounderList %in% c("V","A","P")]))

print("Fitness gains of reconstructed mutations")
print(reconstructedFitnessDFTrunc[c(1:2),])

print("TM-specific founder entropy by gene randomizations")
print(geneEntropyDFTM)

print("Generic founder entropy by gene randomizations")
print(geneEntropyDFGeneric)

print("Fitness defect of E relative to wt E. coli")
EvsWTFitnessDF$Fitness = EvsWTFitnessDF$Fitness*100
EvsWTFitnessDF$stderror = EvsWTFitnessDF$stderror*100
EvsWTFitnessDF$variance = EvsWTFitnessDF$variance*100^2
print(EvsWTFitnessDF[NROW(EvsWTFitnessDF),])

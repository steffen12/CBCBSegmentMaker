#startTime = proc.time();

#Read in files
transcriptomeLibrary = "TxDb.Hsapiens.UCSC.hg19.knownGene";
genomeLibrary = "BSgenome.Hsapiens.UCSC.hg19";
k = 100

library(transcriptomeLibrary, character.only=TRUE);
library(parallel);

# Calculate the number of cores
no_cores <- detectCores()
cl <- makeCluster(no_cores)

#Create Output Directory if it doesn't exist
if(!dir.exists(file.path("SegmentMakerFinalOutput")))
{
  dir.create("SegmentMakerFinalOutput");
}

#Function to return a FASTA string for outputting given a region
addFASTASegmentEntry <- function(region, hg19, currentStrand, currentChromosome, geneID, commonTranscripts, FASTAString)
{
  regionID = paste(">", currentChromosome, "|", geneID,"|", paste(csepommonTranscripts, collapse=","),"|", sep="");
  sequence = "";
  regionGRange = GRanges(seqnames=currentChromosome, ranges=region, strand=currentStrand)
  sequence = paste(getSeq(hg19, regionGRange), collapse = "")
  regionRangesString = paste("[", start(region), ",", end(region), "]-", collapse="", sep="");
  FASTAString = paste(FASTAString, regionID, regionRangesString, "\n", sequence, "\n", sep="");
  return(FASTAString);
}

#This function finds the forks and merges in the transcriptome graph
findStartAndEndStops <- function(transcriptsIndices, geneLength)
{
  numTranscripts = length(transcriptsIndices);
  startStops <- vector(mode="logical", length=geneLength); #Vector element will be true if index comes after a merge
  endStops <- vector(mode="logical", length=geneLength); #Vector element will be false if index comes before a fork
  indices <- rep(2, numTranscripts);#start at 2 because 1 is -1 (root node)
  for(h in 1:geneLength)
  {
    #For each transcript, get a list of the next exon for each exon in a list and another list of the previous exon bin for each exon in a list
    #If for exon #2, there are two different exons that follow on different transcripts, then exon #2 is the start of a fork. If there are two different exons that precede on different transcripts, then exon #2 is the end of a merge.
    currentElements = sapply(1:numTranscripts, function(x) transcriptsIndices[[x]][indices[x]]); #Try vapply to make it faster
    currentIndices = which(currentElements == h);
    nextElements = sapply(1:length(currentIndices), function(x) transcriptsIndices[[currentIndices[x]]][indices[currentIndices[x]]+1]);
    previousElements = sapply(1:length(currentIndices), function(x) transcriptsIndices[[currentIndices[x]]][indices[currentIndices[x]]-1]);
    if(max(nextElements) != min(nextElements) || max(nextElements) == 99999) #All aren't the same
    {
      endStops[h] = TRUE;
    }
    if(max(previousElements) != min(previousElements) || min(previousElements) == -1)
    {
      startStops[h] = TRUE;
    }
    indices[currentIndices] = indices[currentIndices] + 1;
  }
  return(list(startStops, endStops))
}

processGene(geneID, sgf_ucsc, geneGraph, outputFile)
{
  print("GENE");
  print(geneID);
  FASTAString = ""; #FASTAString will contain the segment entries for the whole gene and will be written to the file output at the end
  #Construct all the graphs of the gene
  geneJunctionGraph <- sgf_ucsc[(geneID(sgf_ucsc) == geneID) & (SGSeq::type(sgf_ucsc) == "J")];
  geneGraph <- sgf_ucsc[(geneID(sgf_ucsc) == geneID) & (SGSeq::type(sgf_ucsc) == "E")];
  geneTranscripts <- sort(unique(unlist(txName(geneGraph))));
  emptyExon = any(unlist(lapply(txName(geneGraph), function(element) identical(character(0), element))));
  print(geneTranscripts);
  previousSegments = GRanges();
  
  if(emptyExon)
  {
    print(paste("Error on Gene", geneID));
  }
  else if(length(geneTranscripts) == 1) #If the gene only has one transcript, don't process it normally, just output that transcript as it is
  {
    regionID = paste(">", currentChromosome, ":", geneID,"|",geneTranscripts,"|");
    sequence = "";
    sequence = paste(getSeq(hg19, geneGraph), collapse = "")
    regionRangesString = paste("[", start(ranges(geneGraph)), ",", end(ranges(geneGraph)), "]-", collapse="");
    FASTAString = paste(regionID, regionRangesString, "\n", sequence, "\n");
  }
  else
  {
    #Get Indices of exons used in the transcripts. -1 is a placeholder for the root, and 99999 is a placeholder for the leaf
    transcriptsIndices = lapply(geneTranscripts, function(geneTranscript) c(-1, which(any(geneTranscript == txName(geneGraph))), 99999));#Include placeholders for root and leaf
    transcripts = lapply(transcriptsIndices, function(transcriptIndices) geneGraph[transcriptIndices[c(-1, -length(transcriptIndices))]]);
    startEndStopsList = findStartAndEndStops(transcriptsIndices, length(geneGraph)); #Finds gene indexes of forks and merges
    startStops = startEndStopsList[[1]]; #Unpackage startStops and endStops
    endStops = startEndStopsList[[2]];
    for(i in 1:length(transcripts))
    {
      transcript = ranges(transcripts[[i]]);
      transcriptIndices = transcriptsIndices[[i]][c(-1, -1*length(transcriptsIndices[[i]]))];
      currentTranscript = geneTranscripts[i]; #Current Transcript Name
      currentStrand = as.character(strand(transcripts[[1]][1]));
      rangeLength = 0;
      startIndex = 1; #Start at the first exon
      endIndex = 1;
      #To make the first region, have the first exon be the starting exon, and keep moving down the transcript until the sum of the lengths of all exons is greater than or equal to k
      while(rangeLength + width(transcript[endIndex]) < k && (endIndex < length(transcript))) #If region is smaller than k, should it even be made?
      {
        rangeLength = rangeLength + width(transcript[endIndex]);
        endIndex = endIndex + 1
      }
      #startIndex is the index of the first exon in the region, startBase is the base in that exon which the region begins (first base in this case)
      #endIndex is the index of the last exon in the region, endBase is the base in that exon which the region ends
      startExon = transcript[startIndex];
      startBase= start(startExon);
      endExon = transcript[endIndex];
      endBase = start(endExon) + (k - 1);
      #Subtract from endBase so that endBase is k bp downstream from the start of the transcript and the region represented is k bp long
      if(startIndex != endIndex) 
      {
        endBase = endBase - sum(width(transcript[startIndex:(endIndex-1)]));
      }
      prevStartIndex = startIndex;
      prevStartBase = start(startExon);
      while(endIndex <= length(transcript)) #Keep running until the program attempts to move onto an exon after the last exon
      {
        startExon = transcript[startIndex];
        endExon = transcript[endIndex];
        #These values will be used later, need to be stored now
        startIndexBefore = startIndex; #Start and end index before + 1
        endIndexBefore = endIndex;
        if(end(startExon)-startBase < end(endExon)-endBase) #If for the region, the distance from startBase to the end of the first exon in the region is greater than the distance from endBase to the end of the last exon in the region
        {
          endSegment = FALSE;
          startSegment = TRUE;
          startIndex = startIndex + 1;
          endBase = endBase+end(startExon)-startBase + 1; #+1?
          startBase = start(transcript[startIndex]);
        }
        else if (end(startExon)-startBase > end(endExon)-endBase)  #If for the region, the distance from endBase to the end of the last exon in the region is greater than the distance from startBase to the end of the first exon in the region
        {
          endSegment = TRUE;
          startSegment = FALSE;
          endIndex = endIndex + 1;
          startBase = startBase+(end(endExon) - endBase + 1);
          if(endIndex <= length(transcript))
          {
            endBase = start(transcript[endIndex]);
          }
        }
        else #If for the region, the distance from endBase to the end of the last exon in the region equals the distance from startBase to the end of the first exon in the region
        {
          endSegment = TRUE;
          startSegment = TRUE;
          startIndex = startIndex + 1;
          endIndex = endIndex + 1;
          if(endIndex <= length(transcript))
          {
            endBase = start(transcript[endIndex]);
          }
          startBase = start(transcript[startIndex]);
        }
        #Check previous values of ended and started
        ended = (endSegment && endStops[transcriptIndices[endIndexBefore]]); #True if previous last exon of the region was the start of a fork
        started = (startSegment && startStops[transcriptIndices[startIndexBefore+1]]); #True if previous first exon of the region followed a merge
        if(ended) #If the previous region was segmented once the last exon hit the start of a fork
        {
          region = transcript[prevStartIndex:(endIndexBefore)]; #Region is created from the startIndex of the previous first exon of the region
          currentSegmentTranscripts <- sort(Reduce(intersect, (txName(transcripts[[i]][prevStartIndex:endIndexBefore])))); #Find the intersection of all transcripts of all the elements
          #txName(transcripts[[i]][prevStartIndex])[[1]] #sort(unique(unlist(txName(transcripts[[i]][prevStartIndex:(endIndexBefore)]))));
          start(region[1]) = prevStartBase;
          prevStartBase = startBase;#+1 #Add one to move to next base, as previous one was already covered
          prevStartIndex = startIndex;
          if(currentSegmentTranscripts[1] == currentTranscript) #The currentSegmentTranscripts array should be in the same order for duplicate segments created, so this line ensures that a segment is made only if the currentTranscript is the first transcript in the list
          {
            FASTAString = addFASTASegmentEntry(region, hg19, currentStrand, currentChromosome, geneID, currentSegmentTranscripts, FASTAString)
          }
        }
        else if(started) #If the previous region was segmented once the first exon hit the end of a merge
        {
          region = transcript[prevStartIndex:endIndexBefore]; #Region is created from the startIndex of the previous first exon of the region
          currentSegmentTranscripts <- sort(Reduce(intersect, (txName(transcripts[[i]][prevStartIndex:endIndex])))); #Find the intersection of all transcripts of all the elements
          start(region[1]) = prevStartBase;
          end(region[length(region)]) = endBase-1; #-1?
          prevStartBase = startBase;
          prevStartIndex = startIndex;
          if(currentSegmentTranscripts[1] == currentTranscript) #The currentSegmentTranscripts array should be in the same order for duplicate segments created, so this line ensures that a segment is made only if the currentTranscript is the first transcript in the list
          {
            FASTAString = addFASTASegmentEntry(region, hg19, currentStrand, currentChromosome, geneID, currentSegmentTranscripts, FASTAString)
          }
        }
      }
    }
  }
  write(substr(FASTAString, 1, nchar(FASTAString)-2), file=outputFile, append=TRUE, sep = "");
}

#Call this function once per chromosome to process it and make segments
processChromosome <- function(currentChromosome, k, transcriptomeLibrary, genomeLibrary, addFASTASegmentEntry, findStartAndEndStops)
{
  library(SGSeq);
  library(genomeLibrary, character.only=TRUE);
  library(transcriptomeLibrary, character.only=TRUE);
  txdb <- eval(parse(text = transcriptomeLibrary));
  hg19 <- eval(parse(text = genomeLibrary));
  txdb = restoreSeqlevels(txdb);
  txdbChr <- keepSeqlevels(txdb, currentChromosome);
  sgf_ucsc = convertToSGFeatures(convertToTxFeatures(txdbChr));
  print(currentChromosome)
  outputFile = paste(getwd(), "/SegmentMakerFinalOutput/", currentChromosome, "_Output.fa", sep="");
  close(file(outputFile, open="w")); #Clear output file
  numGenes = max(geneID(sgf_ucsc));
  print(numGenes)
  startGene = 1
  
  for (geneID in startGene:numGenes)
  {
    try(processGene(geneID, sgf_ucsc, geneGraph, outputFile)); #Gene may fail due to issues with SGSeq package
  }
}

txdb <- eval(parse(text = transcriptomeLibrary));
txdb <- restoreSeqlevels(txdb);
chromosomes = seqlevels(txdb)
#nullOutput = bplapply(chromosomes, processChromosome, BPPARAM = MulticoreParam(), txdb=txdb);
nullOutput = parLapply(cl, chromosomes, processChromosome, k=k, transcriptomeLibrary=transcriptomeLibrary, genomeLibrary=genomeLibrary, addFASTASegmentEntry=addFASTASegmentEntry, findStartAndEndStops=findStartAndEndStops); #Code is parallelized

#print(proc.time() - startTime);

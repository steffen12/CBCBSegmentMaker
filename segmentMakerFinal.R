#startTime = proc.time();

transcriptomeLibrary = "TxDb.Hsapiens.UCSC.hg19.knownGene";
genomeLibrary = "BSgenome.Hsapiens.UCSC.hg19";
k = 100

library(transcriptomeLibrary, character.only=TRUE);
library(parallel);

# Calculate the number of cores
no_cores <- detectCores()
cl <- makeCluster(no_cores)

if(!dir.exists(file.path("SegmentMakerFinalOutput")))
{
  dir.create("SegmentMakerFinalOutput");
}

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

findStartAndEndStops <- function(transcriptsIndices, geneLength)
{
  numTranscripts = length(transcriptsIndices);
  startStops <- vector(mode="logical", length=geneLength);
  endStops <- vector(mode="logical", length=geneLength);
  indices <- rep(2, numTranscripts);#start at 2 because 1 is -1 (root node)
  for(h in 1:geneLength)
  {
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
  FASTAString = "";
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
  else if(length(geneTranscripts) == 1)
  {
    regionID = paste(">", currentChromosome, ":", geneID,"|",geneTranscripts,"|");
    sequence = "";
    sequence = paste(getSeq(hg19, geneGraph), collapse = "")
    regionRangesString = paste("[", start(ranges(geneGraph)), ",", end(ranges(geneGraph)), "]-", collapse="");
    FASTAString = paste(regionID, regionRangesString, "\n", sequence, "\n");
  }
  else
  {
    transcriptsIndices = lapply(geneTranscripts, function(geneTranscript) c(-1, which(any(geneTranscript == txName(geneGraph))), 99999));#Include placeholders for root and leaf
    transcripts = lapply(transcriptsIndices, function(transcriptIndices) geneGraph[transcriptIndices[c(-1, -length(transcriptIndices))]]);
    startEndStopsList = findStartAndEndStops(transcriptsIndices, length(geneGraph)); #Finds gene indexes of forks and merges
    startStops = startEndStopsList[[1]];
    endStops = startEndStopsList[[2]];
    for(i in 1:length(transcripts))
    {
      transcript = ranges(transcripts[[i]]);
      transcriptIndices = transcriptsIndices[[i]][c(-1, -1*length(transcriptsIndices[[i]]))];
      currentTranscript = geneTranscripts[i]; #Current Transcript Name
      currentStrand = as.character(strand(transcripts[[1]][1])); #Make cause issues
      rangeLength = 0;
      startIndex = 1;
      endIndex = 1;
      while(rangeLength + width(transcript[endIndex]) < k && (endIndex < length(transcript)))
      {
        rangeLength = rangeLength + width(transcript[endIndex]);
        endIndex = endIndex + 1
      }
      startExon = transcript[startIndex];
      startBase= start(startExon);
      endExon = transcript[endIndex];
      endBase = start(endExon) + (k - 1);
      if(startIndex != endIndex) 
      {
        endBase = endBase - sum(width(transcript[startIndex:(endIndex-1)]));
      }
      prevStartIndex = startIndex;
      prevStartBase = start(startExon);
      while(endIndex <= length(transcript))
      {
        startExon = transcript[startIndex];
        endExon = transcript[endIndex];
        startIndexBefore = startIndex; #Start and end index before + 1
        endIndexBefore = endIndex;
        if(end(startExon)-startBase < end(endExon)-endBase) #Start is closer
        {
          endSegment = FALSE;
          startSegment = TRUE;
          startIndex = startIndex + 1;
          endBase = endBase+end(startExon)-startBase + 1; #+1?
          startBase = start(transcript[startIndex]);
        }
        else if (end(startExon)-startBase > end(endExon)-endBase) #end is closer
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
        else #start and end are equally close
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
        ended = (endSegment && endStops[transcriptIndices[endIndexBefore]]);
        started = (startSegment && startStops[transcriptIndices[startIndexBefore+1]]); #Check if next one is start
        if(ended)
        {
          region = transcript[prevStartIndex:(endIndexBefore)];
          currentSegmentTranscripts <- sort(Reduce(intersect, (txName(transcripts[[i]][prevStartIndex:endIndexBefore]))));#txName(transcripts[[i]][prevStartIndex])[[1]] #sort(unique(unlist(txName(transcripts[[i]][prevStartIndex:(endIndexBefore)]))));
          start(region[1]) = prevStartBase;
          prevStartBase = startBase;#+1 #Add one to move to next base, as previous one was already covered
          prevStartIndex = startIndex;
          if(currentSegmentTranscripts[1] == currentTranscript)
          {
            FASTAString = addFASTASegmentEntry(region, hg19, currentStrand, currentChromosome, geneID, currentSegmentTranscripts, FASTAString)
          }
        }
        else if(started)
        {
          region = transcript[prevStartIndex:endIndexBefore];
          currentSegmentTranscripts <- sort(Reduce(intersect, (txName(transcripts[[i]][prevStartIndex:endIndex]))));#txName(transcripts[[i]][prevStartIndex])[[1]] #sort(unique(unlist(txName(transcripts[[i]][prevStartIndex:(endIndexBefore)]))));
          start(region[1]) = prevStartBase;
          end(region[length(region)]) = endBase-1; #-1?
          prevStartBase = startBase;
          prevStartIndex = startIndex;
          if(currentSegmentTranscripts[1] == currentTranscript)
          {
            FASTAString = addFASTASegmentEntry(region, hg19, currentStrand, currentChromosome, geneID, currentSegmentTranscripts, FASTAString)
          }
        }
      }
    }
  }
  write(substr(FASTAString, 1, nchar(FASTAString)-2), file=outputFile, append=TRUE, sep = "");
}

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
  
  for (geneID in startGene:numGenes)#make numGenes
  {
    try(processGene(geneID, sgf_ucsc, geneGraph, outputFile));
  }
}

txdb <- eval(parse(text = transcriptomeLibrary));
txdb <- restoreSeqlevels(txdb);
chromosomes = seqlevels(txdb)
#nullOutput = bplapply(chromosomes, processChromosome, BPPARAM = MulticoreParam(), txdb=txdb);
nullOutput = parLapply(cl, chromosomes, processChromosome, k=k, transcriptomeLibrary=transcriptomeLibrary, genomeLibrary=genomeLibrary, addFASTASegmentEntry=addFASTASegmentEntry, findStartAndEndStops=findStartAndEndStops);

#print(proc.time() - startTime);

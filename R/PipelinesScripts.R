#' DADA2 pipeline function
#'
#' This function implements the DADA2 sequence QC and annotation pipeline for amplicon sequence libraries. Requires a directory of fastq files that are named with _R1_001.fastq.gz or _R2_001.fastq.gz, formatted with searchable sample IDs incorporated into the file names. Depends on stringr and DADA2. Returns a phyloseq object with metadata.
#' @param in.dir directory path to fastq files
#' @param ref directory path to reference for taxomony annotation
#' @param ID string to search file names for sample ID values
#' @param truncLen where to truncate sequence forward and reverse reads. Default value is c(200, 150). See DADA2 documentation
#' @param trimLeft removes bases from left, used to remove primers left over. Default value is c(18, 19). See DADA2 documentation
#' @param maxN sets maximum number of N values in sequences that are retained. Default value is 0. See DADA2 documentation
#' @param maxEE sets the maximum number of expected errors for each read. Default value is c(2,2). See DADA2 documentation
#' @param truncQ paramter to tune sequence truncation. Default is 2. See DADA2 documentation
#' @param rm.phix logical to remove PhiX. Default is TRUE. See DADA2 documentation
#' @param compress logical. Default is TRUE. See DADA2 documentation
#' @param multithread logical. Default is TRUE. See DADA2 documentation
#' @param metadf metadata for sequence annotation
#' @keywords DADA2 pipeline
#' @export
#' @examples
#' DAD2Apipe()
DADA2pipe<-function(in.dir, ref, ID, truncLen=c(200,150), trimLeft = c(18,19), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE){
  require(dada2)
  require(stringr)
    fnFs<-sort(list.files(in.dir, pattern="_R1_001.fastq.gz", full.names=T))
    fnRs<-sort(list.files(in.dir, pattern="_R2_001.fastq.gz", full.names=T))
    sample.names<-str_match(fnFs, ID)
    filtFs<-file.path(in.dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs<-file.path(in.dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    out<-filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen, trimLeft, maxN, maxEE, truncQ, rm.phix, compress, multithread)
    errF <- learnErrors(filtFs, multithread)
    errR <- learnErrors(filtRs, multithread)
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names

    dadaFs <- dada(derepFs, err=errF, multithread, pool=TRUE)
    dadaRs <- dada(derepRs, err=errR, multithread, pool=TRUE)

    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
    seqtab <- makeSequenceTable(mergers)
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread, verbose=TRUE)

    taxa <- assignTaxonomy(seqtab.nochim, ref, multithread)

    ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
                   sample_data(samdf),
                   tax_table(taxa))
    ps
  }


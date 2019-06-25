
#Function to get the annotated ORF information
#' @title Output gene structure information.
#' @description Take a gtf or gff file to output several GRangesLists for the range information for transcripts, CDS and untranslation regions (UTRs) of the gene of interest to the global environment. These GRangesLists will be used for other functions  of RiboPlotter to define the features of genes.
#' @param annotation Path to the annotation file in either gtf or gff format.
#' @param format Is the annotation file in gtf or gff format?
#' @param dataSource Which organization annotate the genome? For example, TAIR or Araport for Arabidopsis.
#' @param organism Which organism is it? For example, Arabidopsis.
#'
#' @return multiple GRangesLists with transcript, CDS and UTR information.
#' @export
#'
#' @examples
#' \dontrun{
#' gene.structure(annotation="PathToArabidopsis.gff",format="gff",dataSource="Araport",organism="Arabidopsis thaliana")
#' }
#'
gene.structure <- function(annotation,format="gtf",dataSource="",organism=""){
  txdb <- makeTxDbFromGFF(annotation,format=format, dataSource=dataSource,organism=organism)
  exonsByTx <- exonsBy(txdb,by='tx',use.names=T)
  exonsByGene <- exonsBy(txdb,by='gene')
  exonsByTxGene <- exonsBy(txdb,by=c('tx','gene'),use.names=T)
  txByGene <- transcriptsBy(txdb,by='gene')
  cdsByTx <- cdsBy(txdb, by="tx",use.names=T)
  fiveUTR <- fiveUTRsByTranscript(txdb,use.names=T)
  threeUTR <- threeUTRsByTranscript(txdb,use.names=T)
  cds <- cdsBy(txdb, by=c("tx","gene"),use.names=TRUE)

  assign("txdb", txdb, envir = .GlobalEnv)
  assign("exonsByTx", exonsByTx, envir = .GlobalEnv)
  assign("exonsByGene", exonsByGene, envir = .GlobalEnv)
  assign("exonsByTxGene", exonsByTxGene, envir = .GlobalEnv)
  assign("txByGene", txByGene, envir = .GlobalEnv)
  assign("cdsByTx", cdsByTx, envir = .GlobalEnv)
  assign("fiveUTR", fiveUTR, envir = .GlobalEnv)
  assign("threeUTR", threeUTR, envir = .GlobalEnv)
  assign("cds",cds, envir = .GlobalEnv)
}

#Function to get uORF information
#' @title Output gene structure information for uORFs.
#' @description uORFs are not usually in the regular annotation files. Therefore you would need to create a gtf or gff of the uORFs of interest. This function takes a gtf or gff file to output several GRangesLists for the range information for transcripts, CDS and untranslation regions (UTRs) of the uORFs to the global environment. These GRangesLists will be used for other functions of RiboPlotter to define the features of the uORF
#' @param uorf_annotation Path to the uORF annotation file in either gtf or gff format.
#' @param format Is the uORF annotation file in gtf or gff format?
#' @param dataSource Which organization annotate the genome? For example, TAIR or Araport for Arabidopsis.
#' @param organism Which organism is it? For example, Arabidopsis.
#'
#' @return Multiple GRangesLists for uORFs with transcript, CDS and UTR information.
#' @export
#'
#' @examples
#' \dontrun{
#' uorf.structure(annotation="PathToArabidopsis.gff",format="gff",dataSource="Araport",organism="Arabidopsis thaliana")
#' }
#'
uorf.structure <- function(uorf_annotation,format="gtf",dataSource="",organism=""){
  txdb_u <- makeTxDbFromGFF(uorf_annotation,format=format, dataSource=dataSource,organism=organism)
  exonsByTx_u <- exonsBy(txdb_u,by='tx',use.names=T)
  exonsByGene_u <- exonsBy(txdb_u,by='gene')
  exonsByTxGene_u <- exonsBy(txdb_u,by=c('tx','gene'),use.names=T)
  txByGene_u <- transcriptsBy(txdb_u,by='gene')
  cdsByTx_u <- cdsBy(txdb_u, by="tx",use.names=T)
  cds_u <- cdsBy(txdb_u, by=c("tx","gene"),use.names=TRUE)
  assign("txdb_u", txdb_u, envir = .GlobalEnv)
  assign("exonsByTx_u", exonsByTx_u, envir = .GlobalEnv)
  assign("exonsByGene_u", exonsByGene_u, envir = .GlobalEnv)
  assign("exonsByTxGene_u", exonsByTxGene_u, envir = .GlobalEnv)
  assign("txByGene_u", txByGene_u, envir = .GlobalEnv)
  assign("cdsByTx_u", cdsByTx_u, envir = .GlobalEnv)
  assign("cds_u",cds_u, envir = .GlobalEnv)
}

#
#' @title Function to get bam file path, ribo-seq p-site information and axis labels information.
#' @description This function obtains the RNA-seq bam file path, ribo-seq p-site file information, and the labels for y axes of the RNA-seq and ribo-seq plots
#' @param rna1 The first RNA-seq bam file path.
#' @param ribo1 The first ribo-seq tsv file path.
#' @param rna2 The second RNA-seq bam file path.
#' @param ribo2 The second ribo-seq tsv file path.
#' @param RNAlab1 The y-axis label for the first RNA-seq datasets.
#' @param RNAlab2 The y-axis label for the second RNA-seq datasets.
#' @param Ribolab1 The y-axis label for the first ribo-seq datasets.
#' @param Ribolab2 The y-axis label for the second ribo-seq datasets.
#'
#' @return Assign pathes or tsv files to global environment required for downstream analysis
#' @export
#'
#' @examples
#' \dontrun{
#' uorf.structure(uorf_annotation="/Volumes/BACKUP/project2/TAIR10.29.gtf",dataSource="Araport",organism="Arabidopsis thaliana")
#' }
rna_bam.ribo <- function(rna1,ribo1,rna2,ribo2,RNAlab1="RNA_sample1",RNAlab2="RNA_sample2",Ribolab1="Ribo_sample1",Ribolab2="Ribo_sample2"){
  #get path to RNASeq Bam file
  RNAseqBam1 <- rna1
  RNAseqBam2 <- rna2
  #get ribo-seq all p-site information
  riboR1 <- read.delim(file=ribo1,header=F,stringsAsFactors=F,sep="\t")
  colnames(riboR1) <- c("count", "chr", "position", "strand")
  riboR2 <- read.delim(file=ribo2,header=F,stringsAsFactors=F,sep="\t")
  colnames(riboR2) <- c("count", "chr", "position", "strand")
  assign("RNAseqBam1", RNAseqBam1, envir = .GlobalEnv)
  assign("riboR1", riboR1, envir = .GlobalEnv)
  assign("RNAseqBam2", RNAseqBam2, envir = .GlobalEnv)
  assign("riboR2", riboR2, envir = .GlobalEnv)
  assign("RNAlab1", RNAlab1, envir = .GlobalEnv)
  assign("RNAlab2", RNAlab2, envir = .GlobalEnv)
  assign("Ribolab1", Ribolab1, envir = .GlobalEnv)
  assign("Ribolab2", Ribolab2, envir = .GlobalEnv)
}


#' @title plotRanges is a function to plot GRanges
#' @description Plot each isoform based on gtf/gff information
#' @param isoform Integer, the number of isoform
#' @param uORF The name of the uORF
#' @param shortest3UTR The plotGeneModel will provide this information automatically. Providing the shortest 3UTR range allows this function to decide how to plot the shape of the 3'UTR.
#' @param ybottom The ylim for this function.
#' @param main Title for the plot. The gene name is the default.
#' @param colCDS Color for CDS
#' @param col3 Color for 3'UTR
#' @param col5 Color for 5'UTR
#' @return plot each isoform

plotRanges <- function(isoform,uORF=NULL,shortest3UTR,ybottom = 0,main = deparse(substitute(x)),colCDS = "black",col3="white",col5="lightgrey") {
  if(isoform %in% names(cdsByTx)) {
    height <- 0.1
    xlim=ranges(unlist(exonsByTx[isoform]))
    xlimCds=ranges(unlist(cdsByTx[isoform]))

    # plot 5'UTR
    if (isoform %in% names(fiveUTR)) {
      xlim5=ranges(unlist(fiveUTR[isoform]))
      rect(start(xlim5), ybottom, end(xlim5), ybottom + height, col = col5, border = "black")
    }
    # plot lines between exons to represent introns
    if (length(unlist(exonsByTx[isoform]))>1) {
      GAPS <- gaps(unlist(exonsByTx[isoform]),start=NA)
      segments(x0 = start(GAPS),
               y0 = ybottom+height/2,
               x1 = start(GAPS)+width(ranges(GAPS))/2,
               y1 = ybottom+height,
               col = "black",lwd=1)
      segments(x0 = start(GAPS)+width(ranges(GAPS))/2,
               y0 = ybottom+height,
               x1 = end(GAPS),
               y1 = ybottom+height/2,
               col = "black",lwd=1)
    }
    # check strand info
    Frame <- ifelse(as.character(runValue(strand(exonsByTx[isoform])))=="+", firstInFramePSitePerExonPositive(isoform), firstInFramePSitePerExonNegative(isoform))
    # plot
    rect(start(xlimCds), ybottom, end(xlimCds), ybottom + height, col =c("black","black","black")[Frame] , border = "black")
    if (!is.null(uORF)) {
      #The next if check if the isoform contains the range of the uORF, if not, do not plot the uORF in the gene model.
      if(sum(width(setdiff(cdsByTx_u[paste0(uORF,".1")],exonsByTx[isoform])))==0) {
        uORF=paste0(uORF,".1")
        xlim_uORF=ranges(unlist(cdsByTx_u[uORF]))
        rect(start(xlim_uORF), ybottom, end(xlim_uORF), ybottom + height, col ="yellow" , border = "black")
      }
    }
    # print(isoform %in% names(threeUTR))
    ########################
    # Plot 3'UTR with an arrow shap
    if (isoform %in% names(threeUTR)) {
      xlim3=ranges(sort(unlist(threeUTR[isoform])))
      Length=length(unlist(threeUTR[isoform]))
      if (shortest3UTR <=50) {
        z=shortest3UTR
      }
      else{
        z=shortest3UTR/3
      }
      if (as.character(runValue(strand(exonsByTx[isoform])))=="+") {
        if (length(unlist(threeUTR[isoform]))==1) {
          polygon(x=c(start(xlim3), end(xlim3)-z,end(xlim3),end(xlim3)-z,start(xlim3)),y=c(ybottom+height,ybottom+height,ybottom+height/2,ybottom,ybottom), col = col3, border = "black")
        }
        else {
          rect(start(xlim3[1:Length-1]), ybottom, end(xlim3[1:Length-1]), ybottom + height, col = col3, border = "black")
          polygon(x=c(start(xlim3[Length]), end(xlim3[Length])-z,end(xlim3[Length]),end(xlim3[Length])-z,start(xlim3[Length])),y=c(ybottom+height,ybottom+height,ybottom+height/2,ybottom,ybottom), col = col3, border = "black")
        }
      }
      if (as.character(runValue(strand(exonsByTx[isoform])))=="-") {
        if (length(unlist(threeUTR[isoform]))==1) {
          polygon(x=c(start(xlim3),start(xlim3)+z,end(xlim3),end(xlim3),start(xlim3)+z),y=c(ybottom+height/2,ybottom+height,ybottom+height,ybottom,ybottom), col = col3, border = "black")
        }
        else {
          rect(start(xlim3[2:Length]), ybottom, end(xlim3[2:Length]), ybottom+height, col = col3, border = "black")
          polygon(x=c(start(xlim3[1]),start(xlim3[1])+z,end(xlim3[1]),end(xlim3[1]),start(xlim3[1])+z),y=c(ybottom+height/2,ybottom+height,ybottom+height,ybottom,ybottom), col = col3, border = "black")
        }
      }
    }
    axis(1)
  }
  else {
    stop("Input transcript is not a coding gene in gtf/gff file.")
  }
} 

#
#' @title Identify highest value of riboseq reads in a give gene.
#' @description Identify highest value of riboseq reads in a give gene. This is for defining the max Y-axis values.
#' @param GeneName Name of gene used
#' @param isoform Which isoforms used
#' @param ribo A data.frame of riboseq reads
#' @param CDSonly True or False, do we only want to plot the CDS range
#' @param Extend Do we want to plot a wider range. Default is 0.
#'
#' @return  The highest value of riboseq reads
# Do not @export
# No @example

p_site_Y_max <- function(GeneName,isoform,ribo,CDSonly=FALSE,Extend=Extend) {
  #CDSonly=T, then only plot the reads in the CDS
  if(paste0(GeneName,".",isoform,sep = "") %in% names(cdsByTx)) {
    CDS <- cds[paste(GeneName,".",isoform,sep = ""),]
    #find ranges of exons
    Exon <- exonsByGene[GeneName,]
    #Extract chromosome number from CDS object
    chr=as.numeric(as.character(seqnames(unlist(CDS))))[1]
    #Extract strand information from CDS object
    txStrand=as.character(strand(unlist(CDS)))[1]
    #Extract the CDS ranges
    cdsRanges = cdsByTx[names(cdsByTx)==paste0(GeneName,".",isoform,sep = ""),]
    #Extract most left position from the Exon object
    txLeft <-min(start(ranges(unlist(Exon))))
    #Extract most right position from the Exon object
    txRight <-max(end(ranges(unlist(Exon))))
    #Extract most left position from the CDS object
    cdsLeft=min(start(ranges(unlist(CDS))))
    #Extract most right position from the CDS object
    cdsRight=max(end(ranges(unlist(CDS))))
    ##Extract start site from CDS object
    cdsStart=ifelse(txStrand=="+",as.numeric(min(start(ranges(CDS)))),as.numeric(max(end(ranges(CDS)))))
    cdsEnd=ifelse(txStrand=="+",as.numeric(max(end(ranges(CDS)))),as.numeric(min(start(ranges(CDS)))))
    #Extract riboseq reads in the region of the transcript
    if (CDSonly==TRUE) {
      RiboRslt <- ribo[ribo[,2]==chr & ribo[,3] > cdsLeft & ribo[,3] < cdsRight & ribo$strand==txStrand,]
    }
    else if(CDSonly==FALSE) {
      RiboRslt <- ribo[ribo[,2]==chr & ribo[,3] > txLeft-Extend & ribo[,3] < txRight+Extend & ribo$strand==txStrand,]
    }
    if(length(RiboRslt$count)==0) {
      5 #have to give it a number otherwise it will stop here
    } else {
      max(RiboRslt$count)
    }
  }
  else {
    stop("Input transcript is not a coding gene in gtf/gff file.")
  }
}

#p_site_plot_p function to plot periodicity according to transcript info
#' @title p_site_plot_p function to plot periodicity according to transcript info
#' @description Plot Ribo-seq peridicity in three colors according to the frame information.
#' @param GeneName Name of gene used
#' @param isoform Which isoforms used
#' @param ribo riboseq file
#' @param CDSonly TRUE or FALSE
#' @param Extend Integer, plot more at both ends 
#' @param YLIM Integer, max value of Y-axis
#' @param axesQ PLOTc or PLOT, determine if the axes will be plotted
#' @return a plot for the Ribo-seq reads with periodicity 

p_site_plot_p <- function(GeneName,isoform,ribo,CDSonly=FALSE,Extend=Extend,YLIM,axesQ) {
  #CDSonly=T, then only plot the reads in the CDS
  if(paste0(GeneName,".",isoform,sep = "") %in% names(cdsByTx)) {
    CDS <- cds[paste(GeneName,".",isoform,sep = ""),]
    #find ranges of exons
    Exon <- exonsByGene[GeneName,]
    #Extract chromosome number from CDS object
    chr=as.numeric(as.character(seqnames(unlist(CDS))))[1]
    #Extract strand information from CDS object
    txStrand=as.character(strand(unlist(CDS)))[1]
    #Extract the CDS ranges
    cdsRanges = cdsByTx[names(cdsByTx)==paste0(GeneName,".",isoform,sep = ""),]
    #Extract most left position from the Exon object
    txLeft <-min(start(ranges(unlist(Exon))))
    #Extract most right position from the Exon object
    txRight <-max(end(ranges(unlist(Exon))))
    #Extract most left position from the CDS object
    cdsLeft=min(start(ranges(unlist(CDS))))
    #Extract most right position from the CDS object
    cdsRight=max(end(ranges(unlist(CDS))))
    ##Extract start site from CDS object
    cdsStart=ifelse(txStrand=="+",as.numeric(min(start(ranges(CDS)))),as.numeric(max(end(ranges(CDS)))))
    cdsEnd=ifelse(txStrand=="+",as.numeric(max(end(ranges(CDS)))),as.numeric(min(start(ranges(CDS)))))
    #Generate the sequences of positions in the transcript
    if(txStrand=="+") {
      sposition = sort(unlist(mapply(seq,start(unlist(cdsRanges)),end(unlist(cdsRanges)))),decreasing=F)
      sseq=seq(1,length(sposition),1)
      s1 = sposition[which(sseq%%3==1)]
      s2 = sposition[which(sseq%%3==2)]
      s3 = sposition[which(sseq%%3==0)]
    }
    else {
      sposition = sort(unlist(mapply(seq,start(unlist(cdsRanges)),end(unlist(cdsRanges)))),decreasing=T)
      sseq=seq(1,length(sposition),1)
      s1 = sposition[which(sseq%%3==1)]
      s2 = sposition[which(sseq%%3==2)]
      s3 = sposition[which(sseq%%3==0)]
    }
    #Extract riboseq reads in the region of the transcript
    if (CDSonly==TRUE) {
      RiboRslt <- ribo[ribo[,2]==chr & ribo[,3] > cdsLeft & ribo[,3] < cdsRight & ribo$strand==txStrand,]
    }
    else if(CDSonly==FALSE) {
      RiboRslt <- ribo[ribo[,2]==chr & ribo[,3] > txLeft-Extend & ribo[,3] < txRight+Extend & ribo$strand==txStrand,]
    }
    RiboRslt$frame <- factor(ifelse(RiboRslt$position%in%s1,0,ifelse(RiboRslt$position%in%s2,1,ifelse(RiboRslt$position%in%s3,2,3))),levels=c(0,1,2,3))
    if(axesQ=="PLOTc"){
      plot(x=RiboRslt$position,y=RiboRslt$count,type="h",
           xlim=c(txLeft-Extend,txRight+Extend),ylim=c(0,YLIM), #c(0,max(c(0,30))),
           col=c("red","#3366FF","#009900","darkgrey")[RiboRslt$frame],
           lwd=1,xaxt = "n",axes=F)
    }
    else if (axesQ=="PLOT") {
      plot(x=RiboRslt$position,y=RiboRslt$count,type="h",
           xlim=c(txLeft-Extend,txRight+Extend),ylim=c(0,YLIM), #c(0,max(c(0,30))),
           col=c("red","#3366FF","#009900","darkgrey")[RiboRslt$frame],
           lwd=1,xaxt = "n")
    }

    axis(side=1, labels=FALSE, tck = -0.01)
    abline(v=cdsStart,lty=2,lwd=1)
    abline(v=cdsEnd,lty=2,lwd=1, col="darkgrey")
  }
  else {
    stop("Input transcript is not a coding gene in gtf/gff file.")
  }
}

#p_site_plot_p2 function to plot uORF periodicity according to uORF info

#' @title p_site_plot_p2 is for plotting uORFs
#' @description p_site_plot_p2 plots uORFs one at a time and in three colors.
#' @param gene gene ID
#' @param uORF uORF ID
#' @param uORF.isoform numberic, uORF isoform, default use 1 
#' @param ribo riboseq file 
#' @param CDSonly TRUE or FALSE
#' @param Extend Integer, plot more at both ends 
#' @param YLIM Integer, max value of Y-axis
#' @return a plot for the Ribo-seq reads with periodicity

p_site_plot_p2 <- function(gene,uORF,uORF.isoform=1,ribo,CDSonly=TRUE,Extend=Extend,YLIM) {
  #CDSonly=T, then only plot the reads in the CDS
  if(paste0(uORF,".",uORF.isoform,sep = "") %in% names(cdsByTx_u)) {
    CDS <- cds_u[paste(uORF,".",uORF.isoform,sep = ""),]
    #find ranges of exons
    Exon <- exonsByGene_u[uORF,]
    #Extract chromosome number from CDS object
    chr=as.numeric(as.character(seqnames(unlist(CDS))))[1]
    #Extract strand information from CDS object
    txStrand=as.character(strand(unlist(CDS)))[1]
    #Extract the CDS ranges
    cdsRanges = cdsByTx_u[names(cdsByTx_u)==paste0(uORF,".",uORF.isoform,sep = ""),]
    #Extract most left position from the Exon object
    txLeft <-min(start(ranges(unlist(Exon))))
    #Extract most right position from the Exon object
    txRight <-max(end(ranges(unlist(Exon))))
    #Extract most left position from the CDS object
    cdsLeft=min(start(ranges(unlist(CDS))))
    #Extract most right position from the CDS object
    cdsRight=max(end(ranges(unlist(CDS))))
    ##Extract start site from CDS object
    cdsStart=ifelse(txStrand=="+",as.numeric(min(start(ranges(CDS)))),as.numeric(max(end(ranges(CDS)))))
    cdsEnd=ifelse(txStrand=="+",as.numeric(max(end(ranges(CDS)))),as.numeric(min(start(ranges(CDS)))))
    #Generate the sequences of positions in the transcript
    if(txStrand=="+") {
      sposition = sort(unlist(mapply(seq,start(unlist(cdsRanges)),end(unlist(cdsRanges)))),decreasing=F)
      sseq=seq(1,length(sposition),1)
      s1 = sposition[which(sseq%%3==1)]
      s2 = sposition[which(sseq%%3==2)]
      s3 = sposition[which(sseq%%3==0)]
    }
    else {
      sposition = sort(unlist(mapply(seq,start(unlist(cdsRanges)),end(unlist(cdsRanges)))),decreasing=T)
      sseq=seq(1,length(sposition),1)
      s1 = sposition[which(sseq%%3==1)]
      s2 = sposition[which(sseq%%3==2)]
      s3 = sposition[which(sseq%%3==0)]
    }
    #Extract riboseq reads in the region of the transcript
    RiboRslt <- ribo[ribo[,2]==chr & ribo[,3] > cdsLeft & ribo[,3] < cdsRight & ribo$strand==txStrand,]
    RiboRslt$frame <- factor(ifelse(RiboRslt$position%in%s1,0,ifelse(RiboRslt$position%in%s2,1,ifelse(RiboRslt$position%in%s3,2,3))),levels=c(0,1,2,3))
    ###@@@@
    generanges <- ranges(unlist(exonsByGene[gene]))
    genelim <- c(min(start(generanges))-Extend, max(end(generanges))+Extend)
    ###@@@@
    plot(x=RiboRslt$position,y=RiboRslt$count,type="n",xlim=genelim,ylim=c(0,YLIM),lwd=1,xlab=NA, ylab=NA,axes=F)
    lines(x=RiboRslt$position,y=RiboRslt$count,type="h",ylab="Count",xlim=c(txLeft,txRight),YLIM,col=c("white","white","white","white")[RiboRslt$frame],lwd=1.1,xaxt = "n")
    lines(x=RiboRslt$position,y=RiboRslt$count,type="h",ylab="Count",xlim=c(txLeft,txRight),YLIM,col=c("red","#3366FF","#009900","dark grey")[RiboRslt$frame],lwd=1,xaxt = "n")
    abline(v=cdsStart,lty=2,lwd=1, col="black")
    abline(v=cdsEnd,lty=2,lwd=1, col="grey")
  }
  else {
    stop("Input transcript is not a coding gene in gtf/gff file.")
  }
}

#plotGeneModel combines both plotRanges and p_site_plot_all functions

#' @title plot plotGeneModel
#' @description plotGeneModel combines both plotRanges and p_site_plot_all functions
#' @param gene gene ID
#' @param uORF uORF ID
#' @param Extend number of nucleotides to extend on both side of the gene model
#' @return plot the gene model

plotGeneModel <- function(gene,uORF,Extend=Extend){
  isoforms <- length(unlist(txByGene[gene]))
  generanges <- ranges(unlist(exonsByGene[gene]))
  # genelim <- c(min(start(generanges)), max(end(generanges)))
  genelim <- c(min(start(generanges))-Extend, max(end(generanges))+Extend)
  isoforms.w.3UTR <- unlist(txByGene[gene])$tx_name[which(unlist(txByGene[gene])$tx_name %in% names(threeUTR))]
  plot.new()
  yAxis <- (isoforms*0.3+0.1)
  plot.window(genelim,c(0,yAxis))
  for (i in unlist(txByGene[gene])$tx_name) {
    w <- which(unlist(txByGene[gene])$tx_name==i)
    if (i %in% names(threeUTR)) {
      shortest3UTR <- min(sapply(isoforms.w.3UTR, function(j) width(tail(unlist(threeUTR[j]),1))))
      plotRanges(isoform=i,uORF,shortest3UTR,ybottom=(yAxis-0.28*w)) #removed
    }
    else {
      plotRanges(isoform=i,uORF,ybottom=(yAxis-0.3*w))
    }
  }
}

#firstInFramePSitePerExonPositive finds out the frame information with respect to isoform CDS start site when gene is on + strand
#' @title firstInFramePSitePerExonPositive
#' @description if the reads are on the positive strand first in frame P-site per exon
#' @param x Transcript ID
#' @return if the reads are on the positive strand first in frame P-site per exon 

firstInFramePSitePerExonPositive <- function(x){
  nCDS =length(unlist(cdsByTx[x]))
  YFGrange = as.numeric(IRanges:::unlist_as_integer(ranges(unlist(cdsByTx[x]))))
  Range = seq(1,length(YFGrange))
  Seq3 = factor((Range-1)%%3,levels=c(0,1,2))
  dfs = data.frame(YFGrange,Range,Seq3)
  dfs3 = data.frame()
  Listdf = split(dfs,rep(1:nCDS,width(unlist(cdsByTx[x]))))
  for (i in 1:nCDS) {
    dfs3 = rbind(dfs3,Listdf[[i]][Listdf[[i]]$Seq3==0,][1,])
  }
  dfs3$frame <- factor((dfs3$YFGrange- dfs3$YFGrange[1])%%3,levels=c(0,1,2))
  return(dfs3$frame)
}

#firstInFramePSitePerExonNegative finds out the frame information with respect to isoform CDS start site when gene is on - strand

#' @title firstInFramePSitePerExonNegative
#' @description if the reads are on the negative strand first in frame P-site per exon
#' @param x Transcript ID
#' @return if the reads are on the negative strand first in frame P-site per exon 

firstInFramePSitePerExonNegative <- function(x){
  nCDS =length(unlist(cdsByTx[x]))
  YFGrange = sort(IRanges:::unlist_as_integer(ranges(unlist(cdsByTx[x]))),decreasing=T)
  Range = seq(1,length(YFGrange))
  Seq3 = factor((Range-1)%%3,levels=c(0,1,2))
  dfs = data.frame(YFGrange,Range,Seq3)
  Listdf = split(dfs,rep(1:nCDS,width(unlist(cdsByTx[x]))))
  dfs3 = data.frame()
  for (i in 1:nCDS) {
    dfs3 = rbind(dfs3,Listdf[[i]][Listdf[[i]]$Seq3==0,][1,])
  }
  dfs3$frame <- factor((dfs3$YFGrange[1]-dfs3$YFGrange)%%3,levels=c(0,1,2))
  return(dfs3$frame)
}

#
#' @title PLOTt, plot RNA-seq and ribo-seq separately.
#' @description Plot both RNA-seq and ribo-seq in two plots with a gene model.
#' @param YFG Gene ID
#' @param data1 Dataset 1 to plot. Default is RNAseqBam1 that was loaded by rna_bam.ribo.
#' @param CDSonly TRUE or FALSE. Only plot CDS region or all riboseq reads in defined area. Default plot all riboseq reads in the defined area.
#' @param Extend Integer. The number of extra nt ploted at the ends of the plots. 
#' @param isoform Integer. Which isoform to plot periodicity.
#' @param uORF Gene ID for uORF
#' @param NAME Name of the gene
#' @return One plot with RNA-seq and ribo-seq separately.
#' @export
#'
# No @examples
PLOTt <-function(YFG,data1=RNAseqBam1,CDSonly=FALSE,Extend=0,isoform,uORF=NULL,NAME="") {
  transcript_id <- unlist(txByGene[YFG])$tx_name
  #Do not set first transcript because some genes do not have isoform 1
  suppressWarnings(first_transcript <- as.numeric(substring(transcript_id,11))[1])
  if(missing(isoform)) {isoform <- "1"}
  stopifnot(paste0(YFG,".",isoform,sep = "") %in% names(cdsByTx))
  par(mfrow=c(3,1),mar=c(0.2,0.3,0.2,0.2),oma=c(3,4,3,2))
  chr <- as.character(seqnames(exonsByGene[YFG])[[1]])[1]
  generanges <- ranges(unlist(exonsByGene[YFG]))
  GR <- GRanges(seqnames=as.character(chr),IRanges(min(start(generanges))-Extend, max(end(generanges))+Extend),strand=strand(unlist(exonsByGene[YFG]))[1])
  #Add ranges for extracting RNAseq reads
  SZ <- GenomicRanges::reduce(unlist(txByGene[YFG]))
  which1 <- resize(SZ,width=width(SZ)+Extend,fix = "end")
  which1 <- resize(which1,width=width(which1)+Extend,fix = "start")

  # # Here determine anti-sense strand
  # if(as.character(strand(GR))=="+"){AntiSenseStrand <- "-"}
  # else if (as.character(strand(GR))=="-") {AntiSenseStrand <- "+"}
  # else{print("strand info missing")}

  what1 <- c("rname", "strand", "pos", "qwidth","seq")
  param <- ScanBamParam(which = which1, what = what1)
  #Layout
  isoforms <- length(unlist(txByGene[YFG]))
  layout(matrix(c(1,1,2,2,3,3),3,2,byrow=TRUE), widths=c(6,6,6), heights=c(2,2,0.5*isoforms))
  readPairs2 <- readGAlignmentPairs(data1, param=param,strandMode = 2)
  readPairs2 <- readPairs2[strand(readPairs2)==as.character(strand(GR))]
  cvg2 <- coverage(readPairs2)
  Gtx2 <- as.numeric(cvg2[[chr]][ranges(GR)])
  max_Y <- max(max(Gtx2))
  max_P <- max(p_site_Y_max(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR1,Extend=Extend))
  plot(Gtx2,type="h",col="aliceblue",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  par(new = T)
  plot(Gtx2,type="l",col="darkgrey",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  lines(x=c(1,length(Gtx2)),y=c(0,0),col="white",lwd=2)
  legend("topleft",legend="Condition_1_RNAseq",bty="n",cex=1.2)
  p_site_plot_p(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR1,Extend=Extend,YLIM=max_P, axesQ="PLOT")
  par(new = T)
  if (!is.null(uORF)) {p_site_plot_p2(gene=YFG,uORF=uORF,CDSonly=TRUE,uORF.isoform=1,ribo=riboR1,Extend=Extend,YLIM=max_P)}
  legend("topleft",legend="Condition_1_Riboseq",bty="n",cex=1.2)
  mtext(paste('Read counts'), side = 2, outer = T, line = 3.1,cex=0.8)
  #Plot gene model
  plotGeneModel(YFG,Extend=Extend,uORF=uORF)
  mtext(paste(YFG,NAME),side=3,line=0.4, cex=1.2, col="black", outer=TRUE)
}

#PLOTt2, plot 2 sets of RNA-seq and ribo-seq for comparison. It also contains a plot with transcript models.
#' @title PLOTt2, plot 2 sets of RNA-seq and ribo-seq for comparison
#' @description seperately plot 2 sets of RNA-seq and ribo-seq for comparison. It also contains a plot with transcript models.
#' @param YFG Gene ID
#' @param data1 Dataset 1 to plot. Default is RNAseqBam1 that was loaded by rna_bam.ribo.
#' @param data2 Dataset 2 to plot. Default is RNAseqBam2 that was loaded by rna_bam.ribo.
#' @param CDSonly TRUE or FALSE. Only plot CDS region or all riboseq reads in defined area. Default plot all riboseq reads in the defined area.
#' @param Extend Integer. The number of extra nt ploted at the ends of the plots. 
#' @param isoform Integer. Which isoform to plot periodicity.
#' @param uORF Gene ID for uORF
#' @param NAME Name of the gene
#' @return 2 plots for RNAseq and Riboseq in 2 different conditions. 
#' @export

PLOTt2 <-function(YFG,data1=RNAseqBam1,data2=RNAseqBam2,CDSonly=FALSE,Extend=0,isoform,uORF=NULL,NAME="") {
  transcript_id <- unlist(txByGene[YFG])$tx_name
  #Do not set first transcript because some genes do not have isoform 1
  suppressWarnings(first_transcript <- as.numeric(substring(transcript_id,11))[1])
  if(missing(isoform)) {isoform <- "1"}
  stopifnot(paste0(YFG,".",isoform,sep = "") %in% names(cdsByTx))
  
  par(mfrow=c(5,1),mar=c(0.2,0.3,0.2,0.2),oma=c(3,4,3,2))
  chr <- as.character(seqnames(exonsByGene[YFG])[[1]])[1]
  generanges <- ranges(unlist(exonsByGene[YFG]))
  GR <- GRanges(seqnames=as.character(chr),IRanges(min(start(generanges))-Extend, max(end(generanges))+Extend),strand=strand(unlist(exonsByGene[YFG]))[1])
  
  #Add ranges for extracting RNAseq reads
  SZ <- GenomicRanges::reduce(unlist(txByGene[YFG]))
  which1 <- resize(SZ,width=width(SZ)+Extend,fix = "end")
  which1 <- resize(which1,width=width(which1)+Extend,fix = "start")
  
  what1 <- c("rname", "strand", "pos", "qwidth","seq")
  param <- ScanBamParam(which = which1, what = what1)
  #Layout
  isoforms <- length(unlist(txByGene[YFG]))
  layout(matrix(c(1,1,2,2,3,3,4,4,5,5),5,2,byrow=TRUE), widths=c(6,6,6,6,6), heights=c(1.5,1.5,1.5,1.5,0.35*isoforms))
  readPairs1 <- readGAlignmentPairs(data1, param=param,strandMode = 2)
  readPairs1 <- readPairs1[strand(readPairs1)==as.character(strand(GR))]
  cvg1 <- coverage(readPairs1)
  Gtx1 <- as.numeric(cvg1[[chr]][ranges(GR)])
  
  readPairs2 <- readGAlignmentPairs(data2, param=param,strandMode = 2)
  readPairs2 <- readPairs2[strand(readPairs2)==as.character(strand(GR))]
  cvg2 <- coverage(readPairs2)
  Gtx2 <- as.numeric(cvg2[[chr]][ranges(GR)])
  
  max_Y <- max(max(Gtx1),max(Gtx2))
  max_P <- max(max(p_site_Y_max(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR1,Extend=Extend)),
               max(p_site_Y_max(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR2,Extend=Extend)))
  
  plot(Gtx1,type="h",col="aliceblue",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  par(new = T)
  plot(Gtx1,type="l",col="darkgrey",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  lines(x=c(1,length(Gtx1)),y=c(0,0),col="white",lwd=2)
  legend("topleft",legend=RNAlab1,bty="n",cex=1.2)
  p_site_plot_p(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR1,Extend=Extend,YLIM=max_P, axesQ="PLOT")
  par(new = T)
  if (!is.null(uORF)) {p_site_plot_p2(gene=YFG,uORF=uORF,CDSonly=TRUE,uORF.isoform=1,ribo=riboR1,Extend=Extend,YLIM=max_P)}
  legend("topleft",legend=Ribolab1,bty="n",cex=1.2)
  mtext(paste('Read counts'), side = 2, outer = T, line = 3.1,cex=0.8)  
  
  plot(Gtx2,type="h",col="aliceblue",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  par(new = T)
  plot(Gtx2,type="l",col="darkgrey",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  lines(x=c(1,length(Gtx2)),y=c(0,0),col="white",lwd=2)
  legend("topleft",legend=RNAlab2,bty="n",cex=1.2)
  p_site_plot_p(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR2,Extend=Extend,YLIM=max_P, axesQ="PLOT")
  par(new = T)
  if (!is.null(uORF)) {p_site_plot_p2(gene=YFG,uORF=uORF,CDSonly=TRUE,uORF.isoform=1,ribo=riboR2,Extend=Extend,YLIM=max_P)}
  legend("topleft",legend=Ribolab2,bty="n",cex=1.2)
  mtext(paste('Read counts'), side = 2, outer = T, line = 3.1,cex=0.8)
  
  #Plot gene model
  plotGeneModel(YFG,Extend=Extend,uORF=uORF)
  mtext(paste(YFG,NAME),side=3,line=0.4, cex=1.2, col="black", outer=TRUE)
}

#' @title PLOTc plot RNA-seq and ribo-seq together for one datasets.
#' @description PLOTc plot RNA-seq and ribo-seq together for one datasets. It also contains a plot with transcript models.
#' 
#' @param YFG Gene ID
#' @param data1 Dataset 1 to plot. Default is RNAseqBam1 that was loaded by rna_bam.ribo.
#' @param CDSonly TRUE or FALSE. Only plot CDS region or all riboseq reads in defined area. Default plot all riboseq reads in the defined area.
#' @param Extend Integer. The number of extra nt ploted at the ends of the plots. 
#' @param isoform Integer. Which isoform to plot periodicity.
#' @param uORF Gene ID for uORF
#' @param NAME Name of the gene
#' @return Both RNAseq and Riboseq plot together for one set of data
#' @export

PLOTc <-function(YFG,data1=RNAseqBam1,CDSonly=FALSE,Extend=0,isoform,uORF=NULL,NAME="") {
  transcript_id <- unlist(txByGene[YFG])$tx_name
  #Do not set first transcript because some genes do not have isoform 1
  suppressWarnings(first_transcript <- as.numeric(substring(transcript_id,11))[1])
  if(missing(isoform)) {isoform <- "1"}
  stopifnot(paste0(YFG,".",isoform,sep = "") %in% names(cdsByTx))
  par(mfrow=c(2,1),mar=c(0.2,0.3,0.2,0.2),oma=c(3,4,3,4))
  chr <- as.character(seqnames(exonsByGene[YFG])[[1]])[1]
  generanges <- ranges(unlist(exonsByGene[YFG]))
  GR <- GRanges(seqnames=as.character(chr),IRanges(min(start(generanges))-Extend, max(end(generanges))+Extend),strand=strand(unlist(exonsByGene[YFG]))[1])
  #Add ranges for extracting RNAseq reads
  SZ <- GenomicRanges::reduce(unlist(txByGene[YFG]))
  which1 <- resize(SZ,width=width(SZ)+Extend,fix = "end")
  which1 <- resize(which1,width=width(which1)+Extend,fix = "start")

  # Here determine anti-sense strand
  if(as.character(strand(GR))=="+"){AntiSenseStrand <- "-"}
  else if (as.character(strand(GR))=="-") {AntiSenseStrand <- "+"}
  else{print("strand info missing")}

  what1 <- c("rname", "strand", "pos", "qwidth","seq")
  param <- ScanBamParam(which = which1, what = what1)
  #Layout
  isoforms <- length(unlist(txByGene[YFG]))
  layout(matrix(c(1,1,2,2),2,2,byrow=TRUE), widths=c(6,6), heights=c(2.5,0.45*isoforms))

  readPairs <- readGAlignmentPairs(data1, param=param,strandMode = 2)
  readPairs <- readPairs[strand(readPairs)==as.character(strand(GR))]
  cvg <- coverage(readPairs)
  Gtx <- as.numeric(cvg[[chr]][ranges(GR)])
  max_Y <- max(max(Gtx))
  max_P <- max(p_site_Y_max(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR1,Extend=Extend))
  plot(Gtx,type="h",col="aliceblue",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  par(new = T)
  plot(Gtx,type="l",col="darkgrey",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  lines(x=c(1,length(Gtx)),y=c(0,0),col="white",lwd=2)
  par(new = T)
  p_site_plot_p(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR1,Extend=Extend,YLIM=max_P, axesQ="PLOTc")
  par(new=TRUE)
  if (!is.null(uORF)) {p_site_plot_p2(gene=YFG,uORF=uORF,CDSonly=TRUE,uORF.isoform=1,ribo=riboR1,Extend=Extend,YLIM=max_P)}
  axis(side = 4)
  mtext(side = 4, line = 3, 'Condition_1_Riboseq Counts')
  mtext(paste('Condition_1_RNAseq ounts'), side = 2, line = 3)
  #Plot gene model
  plotGeneModel(YFG,Extend=Extend,uORF=uORF)
  mtext(paste(YFG,NAME),side=3,line=0.4, cex=1.2, col="black", outer=TRUE)
}

#PLOTc2, plot 2 sets of RNA-seq and ribo-seq for comparison. It also contains a plot with transcript models.
#' @title PLOTc2
#' @description plot 2 sets of RNA-seq and ribo-seq for comparison. It also contains a plot with transcript models.

#' @param YFG Gene ID
#' @param data1 Dataset 1 to plot. Default is RNAseqBam1 that was loaded by rna_bam.ribo.
#' @param data2 Dataset 2 to plot. Default is RNAseqBam2 that was loaded by rna_bam.ribo.
#' @param CDSonly TRUE or FALSE. Only plot CDS region or all riboseq reads in defined area. Default plot all riboseq reads in the defined area.
#' @param Extend Integer. The number of extra nt ploted at the ends of the plots. 
#' @param isoform Integer. Which isoform to plot periodicity.
#' @param uORF Gene ID for uORF
#' @param NAME Name of the gene
#' @return 2 plots for RNAseq and Riboseq in 2 different conditions. 
#' @export

PLOTc2 <-function(YFG,data1=RNAseqBam1,data2=RNAseqBam2,CDSonly=FALSE,Extend=0,isoform,uORF=NULL,NAME="") {
  transcript_id <- unlist(txByGene[YFG])$tx_name
  #Do not set first transcript because some genes do not have isoform 1
  suppressWarnings(first_transcript <- as.numeric(substring(transcript_id,11))[1])
  if(missing(isoform)) {isoform <- "1"}
  stopifnot(paste0(YFG,".",isoform,sep = "") %in% names(cdsByTx))
  par(mfrow=c(3,1),mar=c(0.2,0.3,0.2,0.2),oma=c(3,4,3,4))
  chr <- as.character(seqnames(exonsByGene[YFG])[[1]])[1]
  generanges <- ranges(unlist(exonsByGene[YFG]))
  GR <- GRanges(seqnames=as.character(chr),IRanges(min(start(generanges))-Extend, max(end(generanges))+Extend),strand=strand(unlist(exonsByGene[YFG]))[1])
  #Add ranges for extracting RNAseq reads
  SZ <- GenomicRanges::reduce(unlist(txByGene[YFG]))
  which1 <- resize(SZ,width=width(SZ)+Extend,fix = "end")
  which1 <- resize(which1,width=width(which1)+Extend,fix = "start")

  # Here determine anti-sense strand
  if(as.character(strand(GR))=="+"){AntiSenseStrand <- "-"}
  else if (as.character(strand(GR))=="-") {AntiSenseStrand <- "+"}
  else{print("strand info missing")}

  what1 <- c("rname", "strand", "pos", "qwidth","seq")
  param <- ScanBamParam(which = which1, what = what1)
  #Layout
  isoforms <- length(unlist(txByGene[YFG]))
  layout(matrix(c(1,1,2,2,3,3),3,2,byrow=TRUE), widths=c(6,6,6), heights=c(2.5,2.5,0.45*isoforms))

  readPairs1 <- readGAlignmentPairs(data1, param=param,strandMode = 2)
  readPairs1 <- readPairs1[strand(readPairs1)==as.character(strand(GR))]
  cvg1 <- coverage(readPairs1)
  Gtx1 <- as.numeric(cvg1[[chr]][ranges(GR)])
  
  
  readPairs2 <- readGAlignmentPairs(data2, param=param,strandMode = 2)
  readPairs2 <- readPairs2[strand(readPairs2)==as.character(strand(GR))]
  cvg2 <- coverage(readPairs2)
  Gtx2 <- as.numeric(cvg2[[chr]][ranges(GR)])
  
  max_Y <- max(max(Gtx1),max(Gtx2))
  max_P <- max(p_site_Y_max(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR1,Extend=Extend),
               p_site_Y_max(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR2,Extend=Extend))
  
  
  plot(Gtx1,type="h",col="aliceblue",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  par(new = T)
  plot(Gtx1,type="l",col="darkgrey",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  lines(x=c(1,length(Gtx1)),y=c(0,0),col="white",lwd=2)
  par(new = T)
  p_site_plot_p(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR1,Extend=Extend,YLIM=max_P, axesQ="PLOTc")
  par(new=TRUE)
  if (!is.null(uORF)) {p_site_plot_p2(gene=YFG,uORF=uORF,CDSonly=TRUE,uORF.isoform=1,ribo=riboR1,Extend=Extend,YLIM=max_P)}
  axis(side = 4)
  mtext(side = 4, line = 3, Ribolab1)
  mtext(RNAlab1, side = 2, line = 3)

  plot(Gtx2,type="h",col="aliceblue",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  par(new = T)
  plot(Gtx2,type="l",col="darkgrey",lwd=1,xaxt='n',ylim=c(0,max_Y+2))
  lines(x=c(1,length(Gtx2)),y=c(0,0),col="white",lwd=2)
  par(new = T)
  p_site_plot_p(YFG,CDSonly=CDSonly,isoform=isoform,ribo=riboR2,Extend=Extend,YLIM=max_P, axesQ="PLOTc")
  par(new=TRUE)
  if (!is.null(uORF)) {p_site_plot_p2(gene=YFG,uORF=uORF,CDSonly=TRUE,uORF.isoform=1,ribo=riboR2,Extend=Extend,YLIM=max_P)}
  axis(side = 4)
  mtext(side = 4, line = 3, Ribolab2)
  mtext(RNAlab2, side = 2, line = 3)

  #Plot gene model
  plotGeneModel(YFG,Extend=Extend,uORF=uORF)
  mtext(paste(YFG,NAME),side=3,line=0.4, cex=1.2, col="black", outer=TRUE)
}













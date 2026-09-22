#5 The setup 

library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("Gviz")

#5.1 Reference download and preparation(programmtic renaming)

# 1. Uncompress and Read the raw NCBI file as simple text
con <- gzfile("data/GCF_000146045.2_R64_genomic.fna.gz", "rt")
ygen <- readLines(con)
close(con)  

# 2. Define the new UCSC-style headers (chrI through chrXVI + chrM)
new_headers <- c(paste0(">chr", as.roman(1:16)), ">chrM")

# 3. Find and Replace
header_lines <- grep(">", ygen)

# Look at the header for chromosome 1 before replacing
print(ygen[header_lines[1]])

# Safety Check: Yeast has 16 nuclear chromosomes + 1 mitochondrial
if(length(header_lines) == 17) {
  ygen[header_lines] <- new_headers
  print("After replacement:")
  print(ygen[header_lines[1]])
  writeLines(ygen, "data/S288C_UCSC.fna")
  message("SUCCESS: Reference renamed to UCSC standards.")
} else {
  stop("Error: Chromosome count mismatch! Check your input file.")
}

#11.3.2 Visualizaing the alignment to chromosome 1 with GRanges

#1
# BAM file path
bam_file <- "outputs/alignment.sorted.bam"

# Open connection to BAM file
bam <- BamFile(bam_file)

# Read all alignments
alns <- readGAlignments(bam)

# Convert to GRanges object
gr_alns <- granges(alns)

# Find the start and end boundaries of all mapped contigs combined
range(gr_alns)

# Peek at the individual contig boundaries
ranges(gr_alns)

#2.
aln_track <- AnnotationTrack(gr_alns, name = "Contigs", genome = "sacCer3", chromosome = "chrI")

# Add genome axis for scale
axis_track <- GenomeAxisTrack()

#3. 
# Plot
plotTracks(list(axis_track, aln_track),
           from = min(start(gr_alns)),
           to = max(end(gr_alns)),
           main = "Contig Alignments to Chromosome 1")
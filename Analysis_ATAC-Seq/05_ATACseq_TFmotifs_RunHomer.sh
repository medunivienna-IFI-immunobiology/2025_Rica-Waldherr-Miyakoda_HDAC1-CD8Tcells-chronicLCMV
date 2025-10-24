# title: run Homer to find enriched transcription factor motifs in genomic regions, specifically DARs
# author: Monika Waldherr
# input: bed files containing consensus regions of interest

# commands to run for the 4 different groups
findMotifsGenome.pl TprogWTopenCons.bed mm10 TprogWTmotifsOutput -size given -p 16 -S 10 -mask
findMotifsGenome.pl TprogKOopenCons.bed mm10 TprogKOmotifsOutput -size given -p 16 -S 10 -mask
findMotifsGenome.pl TtermWTopenCons.bed mm10 TtermWTmotifsOutput -size given -p 16 -S 10 -mask
findMotifsGenome.pl TtermKOopenCons.bed mm10 TtermKOmotifsOutput -size given -p 16 -S 10 -mask

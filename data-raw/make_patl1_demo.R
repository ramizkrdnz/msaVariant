## Generate a small synthetic PATL1 ortholog MSA for demo/test data.
## Real ortholog alignments are obtained from OrthoDB / Ensembl
## Compara; this script creates a plausible 8-species, 60-residue
## window centred on PATL1 residues 490..550 for documentation.
##
## The sequences here are illustrative only -- do NOT use for any
## real biological inference.

set.seed(42)

species <- c(
  "PATL1_HUMAN",         # Homo sapiens
  "PATL1_MOUSE",
  "PATL1_RAT",
  "PATL1_PIG",
  "PATL1_BOVIN",
  "PATL1_CHICK",
  "PATL1_XENTR",         # Xenopus tropicalis
  "PATL1_DANRE"          # Danio rerio
)

# A 61-residue stretch (positions 490..550 in the reference) with
# K518 highly conserved across mammals (positions 29 of the stretch),
# variable in lower vertebrates, and gap regions inserted in
# zebrafish.
human  <- "EELRKSGAEKVDPNALQALKQDLEFKKPKEYAESLRMDFAHRTPSAEYPKDQVNQSSPDYK"
mouse  <- "EELRKSGAEKVDPNALQTLKQDLEFKKPKEYTESLRMEFAHRTPSAEYPKDQVNQSSPDYK"
rat    <- "EELRKSGAEKVDPNALQTLKQDLEFKKPKEYTESLRMEFAHRTPSAEYPKDQVNQSSPDYK"
pig    <- "EELRKSGAEKLDPNALQALKQDLEFKRPKEYAESLRMDFAHRTPSAEYPKDQVNQSSPDYK"
cow    <- "EELRKSGAEKLDPNALQALKQDLEFKRPKEYAESLRMDFAHRTPSAEYPKDQVNQSSPDYK"
chick  <- "EELKKSGADKVDSNAVQALKQDLEFKKPKEYAESVRMDFAHKTPSAEYPRDQVNQSSPDYR"
xen    <- "DELKKQGAEKLDANSLQTLKQELEFKKPHEYAESIRMDFTHRTPSAEYPHDQVNQ-SPEYK"
zfish  <- "DDLKKQGADKVDQNSLQTLRQELDFKKPNDYTESIKMDFTHRSPNADYPRDQVTQSSP--K"

# Sanity: all 61 chars
seqs <- list(PATL1_HUMAN = human, PATL1_MOUSE = mouse, PATL1_RAT = rat,
             PATL1_PIG = pig, PATL1_BOVIN = cow, PATL1_CHICK = chick,
             PATL1_XENTR = xen, PATL1_DANRE = zfish)
stopifnot(all(nchar(unlist(seqs)) == 61))

out_path <- "inst/extdata/patl1_orthologs.fasta"
con <- file(out_path, "w")
on.exit(close(con))
for (sp in names(seqs)) {
  cat(sprintf(">%s\n%s\n", sp, seqs[[sp]]), file = con)
}
message("Wrote ", out_path)

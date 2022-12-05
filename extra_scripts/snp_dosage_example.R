# Example how to extract APOE SNP dosages and info
# from Gen-R HRC 1.1 imputed data

library(indepth)

# Read in HRC map file to lookup coordinates of SNPs
# from http://www.haplotype-reference-consortium.org/site
map <- fread("HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz")
setnames(map, "#CHROM", "CHR") #change names

# APOE SNP RS IDs
snp <- c("rs429358","rs7412")

# Path to vcf file
vcf_path <- "~/GENR3/Imputed/HRC_release1.1/chr"

# Extract SNP dosage for all SNPs specified in SNP
dosage.list <- lapply(snp, function(x) snp_dosage(x, map, vcf_path, n = 5756))

# Combine list of SNP dosage
dosage.data <- do.call("cbind", dosage.list)

# Remove the GENR3_ part from participant identifier
# to be consistent with other Gen-R data
dosage.data$IDC <- gsub("GENR3_", "", row.names(dosage.data))

# Extract SNP info for all SNPs specified in SNP
info.data <- lapply(snp, function(x) snp_info(x, map, vcf_path))

# Combine list of SNP dosage
info.data <- do.call("rbind", info.data)


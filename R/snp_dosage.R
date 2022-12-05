#' \code{snp_dosage} extracts SNP dosages from a VCF file given RSID.
#'
#' @param rsid A string containing the SNP rsid
#' @param map A data.frame containing ID, CHR, POS (rsid, chromosome, position)
#' @param vcf_path Path to vcf file, including prefix, excluding ".dose.vcf.gz"
#' @param n maximum sample size (optional, but useful in case of missing SNPs)
#' 
#' @return A data.frame containing allele dosages for the SNP
#' 
#' @examples
#' # Example how to extract APOE SNP dosages
#' # from Gen-R HRC 1.1 imputed data
#' library(indepth)
#'
#' # Read in HRC map file to lookup coordinates of SNPs
#' # from http://www.haplotype-reference-consortium.org/site
#' map <- fread("HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz")
#' setnames(map, "#CHROM", "CHR") #change names
#' 
#' # APOE SNP RS IDs
#' snp <- c("rs429358","rs7412")
#' 
#' # Path to vcf file
#' vcf_path <- "~/GENR3/Imputed/HRC_release1.1/chr"
#' 
#' # Extract SNP dosage for all SNPs specified in SNP
#' dosage.list <- lapply(snp, function(x) snp_dosage(x, map, vcf_path, n = 5756))
#' 
#' # Combine list of SNP dosage
#' dosage.data <- do.call("cbind", dosage.list)
#' 
#' # Remove the GENR3_ part from participant identifier
#' # to be consistent with other Gen-R data
#' dosage.data$IDC <- gsub("GENR3_", "", row.names(dosage.data))
#'
#' @export

snp_dosage <- function(rsid, map, vcf_path, n=1) {
  tryCatch({
    #Extract chromosome from map file
    chromosome <- map[ID==rsid, CHR]
    #Extract position from map file
    position <- map[ID==rsid, POS]
    #Define path to dose file and corresponding tabix
    tabix.file <- 
      Rsamtools::TabixFile(paste0(vcf_path, chromosome, ".dose.vcf.gz"))
    #Define position to extract
    rng <- GenomicRanges::GRanges(as.character(chromosome), position, "*")
    #Read VCF file specified with TabixFile(), but only at the SNP position defined by GRanges
    vcf <- VariantAnnotation::readVcf(tabix.file, "hg19", rng)
    #Extract dosage
    dosage <- as.data.frame(t(VariantAnnotation::geno(vcf)[[2]]))
    names(dosage) <- rsid
    return(dosage)
  }, error=function(err) {
    message("SNP not found: ", rsid)
    dosage <- data.frame(rep(NA, n))
    names(dosage) <- rsid
    return(dosage)
  })
}
#' \code{snp_info} extracts SNP information
#'
#' @param rsid A string containing the SNP rsid
#' @param map A data.frame containing ID, CHR, POS (rsid, chromosome, position)
#' @param vcf_path Path to vcf file, including prefix, excluding ".dose.vcf.gz"
#'
#' @return A data.frame containing effect/reference allele,
#'  effect allele frequency and imputation quality (R^2)
#'  
#' @examples 
#' # Example how to extract APOE SNP dosages and info
#' # from Gen-R HRC 1.1 imputed data
#' 
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
#' # Extract SNP info for all SNPs specified in SNP
#' info.data <- lapply(snp, function(x) snp_info(x, map, vcf_path))
#'
#' # Combine list of SNP dosage
#' info.data <- do.call("rbind", info.data)
#'
#' @export

snp_info <- function(rsid, map, vcf_path) {
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
    #Read VCF file specified with TabixFile(),
    #but only at the SNP position defined by GRanges
    vcf <- VariantAnnotation::readVcf(tabix.file, "hg19", rng)
    snpid <- rsid #RSID
    #Effect Allele
    effect_allele <- as.character(VariantAnnotation::alt(vcf)[[1]][[1]])
    #Reference Allele
    reference_allele <- as.character(VariantAnnotation::ref(vcf))
    #EAF
    eaf <- as.numeric(VariantAnnotation::info(vcf)[[1]])
    #RSQ
    rsq <- as.numeric(VariantAnnotation::info(vcf)[[3]])
    #Combine all information to a data.frame
    data.frame(snpid, effect_allele, reference_allele, eaf, rsq)
  }, error=function(err) {
    dosage <- data.frame(rsid, NA, NA, NA, NA)
    names(dosage) <- 
      c("snpid", "effect_allele", "reference_allele", "eaf", "rsq")
    return(dosage)
  })
}
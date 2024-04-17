


load.basilica <- function(x, type) {
  fixed <- basilica:::get_fixed_signatures(x, types = type, matrix = TRUE)[[1]]
  denovo <- basilica:::get_denovo_signatures(x, types = type, matrix = TRUE)[[1]]
  exposure <- basilica:::get_exposure(x, types = type, matrix = TRUE)[[1]]
  return(list(fixed=fixed, denovo=denovo, exposure=exposure))
}

load.serena <- function(reference_path, exposure_path, organ_type) {
  
  # ------------------------------ REFERENCE -----------------------------------
  # load serena reference data.frame
  reference.T <- read.table(file = reference_path, sep = '\t', header = TRUE)
  reference <- as.data.frame(t(reference.T)) # data.frame (wide)
  
  # ------------------------------ SIGNATURES ----------------------------------
  # common signatures (9)
  common_names <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS13", "SBS17", "SBS18", "SBS127")
  # rare signatures  (19)
  rare_names <- c("SBS4", "SBS6", "SBS7a", "SBS23", "SBS26", "SBS30", "SBS31", "SBS33", "SBS44", 
                  "SBS57", "SBS87", "SBS94", "SBS96", "SBS105", "SBS114", "SBS116", "SBS123", "SBS124", 
                  "SBS136")
  
  common <- reference[common_names, ]
  rare <- reference[rare_names, ]
  
  # ------------------------------- EXPOSURE -----------------------------------
  # load exposure data.frame (all - unnormalized)
  exposure <- read.table(file = exposure_path, sep = '\t', header = TRUE)
  # filter exposure based on organ type (unnormalized)
  exposure <- subset(exposure, organ==organ_type)
  # removing the cohort and organ column  (unnormalized)
  exposure <- exposure[!(names(exposure) %in% c("cohort", "organ"))]
  # removing columns (signatures) where all are zeros (unnormalized)
  exposure <- exposure %>% dplyr::select(where(~ any(. != 0))) # dplyr required
  # normalizing
  exposure <- sweep(exposure, 1, rowSums(exposure), "/")
  # selecting exposures involved in specified organ
  exposure <- exposure[,c(common_names, rare_names)]
  
  return(list(common=common, rare=rare, exposure=exposure))
}

#-------------------------------------------------------------------------------

input.qc <- function(fixed, denovo, BasilicaExposure, common, rare, SerenaExposure, reference) {
#run.QC <- function(basilica, serena, reference) {
  
  a <- all((colnames(fixed) %in% colnames(reference))==TRUE)
  b <- all((colnames(denovo) %in% colnames(reference))==TRUE)
  c <- all((colnames(common) %in% colnames(reference))==TRUE)
  d <- all((colnames(rare) %in% colnames(reference))==TRUE)
  
  if (!(a & b & c & d)) {
    print("ERROR!")
    return(NULL)
  }
  
  fixed <- fixed[, colnames(reference)]
  denovo <- denovo[, colnames(reference)]
  common <- common[, colnames(reference)]
  rare <- rare[, colnames(reference)]
  
  # extracting common samples between basilica and serena exposure matrices
  intersected_samples <- intersect(rownames(BasilicaExposure), rownames(SerenaExposure))
  BasilicaExposure <- BasilicaExposure[intersected_samples, ]
  SerenaExposure <- SerenaExposure[intersected_samples, ]
  
  #return(list(basilica=basilica, serena=serena))
  return(
    list(
    fixed=fixed, denovo=denovo, BasilicaExposure=BasilicaExposure, 
    common=common, rare=rare, SerenaExposure=SerenaExposure
    )
  )
}

#-------------------------------------------------------------------------------

# input
#   x -----> basilica object
#   type --> ("SBS", "DBS")
# output
#   list of:
#   - fixed
#   - denovo
#   - BasilicaExposure
#   - common
#   - rare
#   - SerenaExposure
load.all <- function(basilica_obj, serena_ref_path, serena_exp_path, organ, reference, type) {
  
  basilica <- load.basilica(basilica_obj, type = type)
  
  serena <- load.serena(
    reference_path=serena_ref_path, 
    exposure_path=serena_exp_path, 
    organ_type=organ
    )
  
  #qc <- run.QC(basilica=basilica, serena=serena, reference=basilica::COSMIC_sbs)
  qc <- input.qc(
    fixed = basilica$fixed, 
    denovo = basilica$denovo, 
    BasilicaExposure = basilica$exposure, 
    common = serena$common, 
    rare = serena$rare, 
    SerenaExposure = serena$exposure, 
    reference = reference
    )
  
  return(list(
    fixed=qc$fixed, denovo=qc$denovo, BasilicaExposure=qc$BasilicaExposure, 
    common=qc$common, rare=qc$rare, SerenaExposure=qc$SerenaExposure
    )
  )
}
  
  
  


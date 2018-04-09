#' Convert fcs files to a data frame
#'
#' A function that convert fcs files to a data frame.
#'
#' @param fcsFiles A vector specifying the location of fcs files (relative to
#'   working directory).
#' @param y A vector containing the clinical outcome of each sample. Must have
#'   the same length as fcsFiles. Null for testing data.
#' @param assay Either "FCM" or "CyTOF" to indicate the type of cytometry data.
#' @param b A positive number used to specify the arcsinh transformation. f(x) =
#'   asinh (b*x) where x is the original value and f(x) is the value after
#'   transformation. The suggested value is 1/150 for flow cytometry (FCM) data
#'   and 1/8 for CyTOF data.
#' @param fileSampleSize An integer specifying the number of events sampled from
#'   each fcs file. If NULL, all the events will be pre-processed and wrote out
#'   to the new fcs files.
#' @param excludeTransformParameters A vector specifying the name of parameters
#'   not to be transformed (left at linear scale).
#' @param compFiles A vector specifying the paths of user supplied compensation
#'   matrix for each fcs file. The matrix must be stored in csv files.
#' @param nameDict A vector used to change marker names.Each element in the
#'   vector is the prefered name of a marker. The name of each element is the
#'   marker name used in the fcs file. For example, a vector
#'   c("CD8b"="CD8","cd8"="CD8") will change "CD8b" and "cd8" into "CD8", making
#'   annotations more consistent.
#' @return Returns a data frame containing the preprocessed cytometry data.
#'   Cells from different fcs files are combined into one flow frame. A new
#'   column, xSample, is introduced to indicate the origin of each cell. The
#'   data frame also includes the clinical outcome y.
#' @importFrom flowCore read.flowSet fsApply keyword logTransform transformList
#'   transform compensate arcsinhTransform colnames
#' @importFrom utils read.csv
#' @examples
#' # Find the table containing fcs file names in CytoDx package
#' path <- system.file("extdata",package="CytoDx")
#' # read the table
#' fcs_info <- read.csv(file.path(path,"fcs_info.csv"))
#' # Specify the path to the cytometry files
#' fn <- file.path(path,fcs_info$fcsName)
#' # Read cytometry files using fcs2DF function
#' train_data <- fcs2DF(fcsFiles=fn,
#'                     y=fcs_info$Label,
#'                     assay="FCM",
#'                     b=1/150,
#'                     excludeTransformParameters=
#'                       c("FSC-A","FSC-W","FSC-H","Time"))
#' @export
fcs2DF <- function(fcsFiles,
                y=NULL,
                assay=c("FCM", "CyTOF"),
                b=1/200,
                fileSampleSize=5000,
                compFiles=NULL,
                nameDict=NULL,
                excludeTransformParameters=
                  c("FSC-A","FSC-W","FSC-H","Time","Cell_length")){

  fcs_param <- NULL
  fcs_names <- gsub(".*/","",fcsFiles)
  excludeTransformParameters <- paste(excludeTransformParameters,collapse="|")

  # 1) identify sample :
  fcsFiles <- as.character(fcsFiles)

  # 2) read the fcs files
  fcs <- flowCore::read.flowSet(fcsFiles,transformation="linearize",alter.names=FALSE,truncate_max_range=FALSE)
  if(!is.null(fileSampleSize)){
    fcs <- flowCore::fsApply(fcs, function(f){
      L <- nrow(f@exprs)
      if(L>fileSampleSize){
        f@exprs <- f@exprs[sample(seq_len(L),fileSampleSize),]
      }
      return(f)
    })
  }

  # 3) compensation and transformation
  if (assay == "FCM") {
    # retreiving fluorescent biomarkers
    w <- which(!grepl(excludeTransformParameters,flowCore::colnames(fcs),ignore.case = TRUE))
    biomarker_vector <- flowCore::colnames(fcs)[w]
    # identify the fcs file format
    version <- flowCore::fsApply(fcs, function(frame) {
      return(flowCore::keyword(frame, "FCSversion")[[1]])
    })
    unique_version <- as.numeric(unique(as.vector(version)))

    # compensation, transformation
    ## if the file version is 2.0, then code will log transform the data
    if (unique_version == 2) {
      trans <- flowCore::logTransform()
      translist <- flowCore::transformList(biomarker_vector,trans)
      fcs <- flowCore::transform(fcs, translist)
    } else if (unique_version == 3) {
      # check if user have provided the compmatrix
      if(!is.null(compFiles)){
        compList=lapply(compFiles,function(x){
          if(is.na(x)){return(as.matrix(NA))}else{return(as.matrix(read.csv(x,row.names = 1,check.names = FALSE)))}
        })
        for(id in seq_along(fcs)){
          if(!is.na(compList[[id]][1,1])){
            fcs[[id]]@description$SPILL <- compList[[id]]
          }
        }
      }

      ### function first checks if the spill matrix in the flowframe is an actual spill matrix or an identity matrix
      if (is.null(flowCore::keyword(fcs[[1]], "SPILL")[[1]]) == FALSE) {
        check <- flowCore::fsApply(fcs, function(x) {
          result <- isSymmetric(flowCore::keyword(x, "SPILL")[[1]])
        })
        ### if all flowframes have spill matrices then my function
        ### will apply them to each flowframe in the flowset for compensation
        if (unique(check)[[1]] == FALSE) {
          fcs <- flowCore::fsApply(fcs, function(x) {
            new_frame <- flowCore::compensate(x, flowCore::keyword(x, "SPILL")[[1]])
            return(new_frame)
          })
        }
      }

      trans <- flowCore::arcsinhTransform(transformationId="defaultArcsinhTransform",a=0,b=b,c=0)
      translist <- flowCore::transformList(biomarker_vector, trans)
      fcs <- flowCore::transform(fcs, translist)
    }
  } else if (assay == "CyTOF") {
    w <- which(!grepl(excludeTransformParameters,flowCore::colnames(fcs),ignore.case = TRUE))
    biomarker_vector  <-  flowCore::colnames(fcs)[w]

    trans <- flowCore::arcsinhTransform(transformationId = "defaultArcsinhTransform", a = 0, b = b, c = 0)
    translist <- flowCore::transformList(biomarker_vector, trans)
    fcs <- flowCore::transform(fcs, translist)
  }
  fcs <- set2DF(fcs,fcsFiles,y)
  if(!is.null(nameDict)){
    t1 <- colnames(fcs)
    w <- which(t1%in%names(nameDict))
    t1[w] <- nameDict[t1[w]]
    colnames(fcs) <- t1
  }
  return(fcs)
}


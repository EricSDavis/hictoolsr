#' Merged loop calls with contact frequency from replicate Hi-C experiments
#'
#' Using data from "Phase separation drives aberrant chromatin looping and
#' cancer development" by Ahn et al. 2021 (GEO links GSE143465 and GSE143465),
#' looping interactions were called with SIP (https://github.com/PouletAxel/SIP).
#' Resulting loop calls were merged using functions from `hictoolsr` and
#' interactions were extracted from replicate `.hic` files. The resulting file
#' can be used to demonstrate differential loop calling between two biological
#' conditions.
#'
#' @format a GInteractions object with metadata columns describing
#'   loop properties and interaction counts from replicate
#'   Hi-C data
#'
#' @docType data
#'
#' @usage data("loopCounts")
#'
#' @references Ahn, Jeong Hyun, Eric S. Davis, Timothy A. Daugird,
#'  Shuai Zhao, Ivana Yoseli Quiroga, Hidetaka Uryu, Jie Li,
#'  et al. 2021. “Phase Separation Drives Aberrant Chromatin
#'  Looping and Cancer Development.” Nature 595 (7868): 591–95.
#'
#'  SIP (https://github.com/PouletAxel/SIP)
"loopCounts"

#' Aggregate peak analysis of differential loop results from WT Hi-C
#'
#' Differential loops between WT and FS conditions were called using data from
#' "Phase separation drives aberrant chromatin looping and cancer development"
#' by Ahn et al. 2021 (GEO links GSE143465 and GSE143465). Short interactions
#' were filtered out from all loops, WT-specific loops and FS-specific loops
#' and a 21x21 matrix of KR-normalized pixels at 10-kb resolution were
#' aggregated for each loop group from the replicate-merged Hi-C map of the
#' WT condition. Data is stored in "data/loopApaWtHic.rda".
#'
#' @format The data object contains a list of 3 APA matrices extracted from
#'   a replicate-merged Hi-C map of the WT condition corresponding to all loops,
#'   WT-specific loops, and FS-specific loops.
#'
#' @docType data
#'
#' @usage data("loopApaWtHic")
#'
#' @references Ahn, Jeong Hyun, Eric S. Davis, Timothy A. Daugird,
#'  Shuai Zhao, Ivana Yoseli Quiroga, Hidetaka Uryu, Jie Li,
#'  et al. 2021. “Phase Separation Drives Aberrant Chromatin
#'  Looping and Cancer Development.” Nature 595 (7868): 591–95.
"loopApaWtHic"

#' Aggregate peak analysis of differential loop results from FS Hi-C
#'
#' Differential loops between WT and FS conditions were called using data from
#' "Phase separation drives aberrant chromatin looping and cancer development"
#' by Ahn et al. 2021 (GEO links GSE143465 and GSE143465). Short interactions
#' were filtered out from all loops, WT-specific loops and FS-specific loops
#' and a 21x21 matrix of KR-normalized pixels at 10-kb resolution were
#' aggregated for each loop group from the replicate-merged Hi-C map of the
#' FS condition. Data is stored in "data/loopApaFsHic.rda".
#'
#' @format The data object contains a list of 3 APA matrices extracted from
#'   a replicate-merged Hi-C map of the FS condition corresponding to all loops,
#'   WT-specific loops, and FS-specific loops.
#'
#' @docType data
#'
#' @usage data("loopApaFsHic")
#'
#' @references Ahn, Jeong Hyun, Eric S. Davis, Timothy A. Daugird,
#'  Shuai Zhao, Ivana Yoseli Quiroga, Hidetaka Uryu, Jie Li,
#'  et al. 2021. “Phase Separation Drives Aberrant Chromatin
#'  Looping and Cancer Development.” Nature 595 (7868): 591–95.
"loopApaFsHic"

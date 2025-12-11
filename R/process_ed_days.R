#' Preprocess ED Data
#'
#' This function processes daily patient data to include or exclude Emergency Department (ED) rows based on different inclusion rules.
#'
#' @param daily_data A data frame containing daily patient data with columns `unique_pt_id`, `seqnum`, `day`, `death`, `ALL_DAYS`, and other clinical variables.
#' @param ed_inclusion An integer (1–4) specifying how to preprocess the ED rows in daily_data:
#'   \describe{
#'     \item{1}{ASE Toolkit default: keep 2 days in ED (i.e., rows with `day > -2`.)}
#'     \item{2}{Include all ED rows (no filtering).}
#'     \item{3}{Drop all ED rows (i.e., `day < 1`).}
#'     \item{4}{Drop ED rows, but if any ED blood culture is detected,
#'              then force `bcx_daily = 1` on hospital day 1.}
#'   }
#' @return A modified `daily_data` data frame with ED rows included or excluded
#' @examples
#' # Example daily data frame
#' daily_data <- data.frame(
#'   unique_pt_id = c(1,1,1, 2,2, 3,3,3, 4,4,4,4,4),
#'   seqnum       = c(1111,1111,1111, 2222,2222, 3333,3333,3333, 4444,4444,4444,4444,4444),
#'   day          = c(0,1,2, 0,1, 0,1,2, 0,1,2,3,4),
#'   death        = c(0,0,0, 0,1, 0,0,0, 0,0,0,0,0),
#'   ALL_DAYS     = c(3,3,3, 2,2, 3,3,3, 5,5,5,5,5),
#'   bcx_daily    = c(1,0,0, 0,0, 1,0,0, 0,1,0,0,0),
#'   vasop_daily  = c(0,1,0, 0,0, 0,1,0, 0,0,1,0,0),
#'   imv_daily    = c(0,0,1, 0,1, 0,0,1, 0,0,1,0,0),
#'   lact_daily_hi = c(1.5,2.5,1.8, 1.9,2.1, 1.5,2.5,1.8, 1.5,1.9,2.1,2.1,3.0),
#'   tbili_daily_hi = c(30,35,40, 32,36, 30,35,40, 30,32,35,36,50),
#'   tbili_daily_lo = c(30,35,40, 32,36, 30,35,40, 30,32,35,36,50),
#'   tbili_baseline = c(20,20,20, 25,25, 20,20,20, 25,25,25,25,25),
#'   creat_daily_hi = c(40,45,50, 35,60, 40,45,50, 35,35,45,60,55),
#'   creat_daily_lo = c(40,45,50, 35,60, 40,45,50, 35,35,45,60,55),
#'   creat_baseline = c(20,20,20, 25,25, 20,20,20, 25,25,25,25,25),
#'   plt_daily_hi   = c(150,80,90, 110,70, 150,80,90, 150,140,70,60,50),
#'   plt_daily_lo   = c(150,80,90, 110,70, 150,80,90, 150,140,70,60,50),
#'   plt_baseline   = c(200,200,200, 180,180, 200,200,200, 180,180,180,180,180),
#'   esrd_icd       = c(0,0,0, 0,0, 0,0,0, 0,0,0,0,0),
#'   new_abx_start  = c(0,1,0, 0,1, 0,1,0, 0,1,0,0,0),
#'   abx_daily      = c(0,1,1, 0,1, 0,1,1, 0,1,1,1,1)
#' )
#' process_ed_days(daily_data, ed_inclusion = 1)
#' @export
process_ed_days <- function(daily_data, ed_inclusion = 1) {
  
  if (ed_inclusion == 1) {
    # ASE Toolkit default: keep only 2 days in ED
    daily_data <- daily_data[daily_data$day > -2, ]
  }
  
  if (ed_inclusion == 2) {
    # include all ED days
    return(daily_data)
  }
  
  if (ed_inclusion == 3) {
    # drop all ED days
    daily_data <- daily_data[daily_data$day >= 1, ]
  }
  
  if (ed_inclusion == 4) {
    # drop all ED days, but:
    # if any bcx in ED, then force bcx_daily == 1 on day 1
    
    needed_cols <- c("seqnum", "day", "bcx_daily")
    if (!all(needed_cols %in% names(daily_data))) {
      warning(
        "handle_ed_inclusion(ed_inclusion = 4): missing required columns: ",
        paste(setdiff(needed_cols, names(daily_data)), collapse = ", "),
        ". Skipping ER BCx logic and just truncating to day >= 1."
      )
      daily_data <- daily_data[daily_data$day >= 1, ]
      return(daily_data)
    }
    
    if ("bcx_er" %in% names(daily_data)) {
      ids_with_ed_bcx <- unique(
        daily_data$seqnum[daily_data$bcx_er == 1 & daily_data$day < 1]
      )
    } else {
      # fallback: infer ED bcx from day < 1
      ids_with_ed_bcx <- unique(
        daily_data$seqnum[daily_data$bcx_daily == 1 & daily_data$day < 1]
      )
    }
    
    idx_day1 <- daily_data$seqnum %in% ids_with_ed_bcx & daily_data$day == 1
    daily_data$bcx_daily[idx_day1] <- 1L
    
    # truncate ED days (day < 1)
    
    daily_data <- daily_data[daily_data$day >= 1, ]
  }
  
  return(daily_data)
}

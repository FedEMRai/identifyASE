#' Define New Antimicrobial Definition
#'
#' This function creates the variable `new_abx_start` to indicate the start of new antimicrobial administration.
#'
#' @param daily_data A data frame containing daily patient data with columns `unique_pt_id`, `seqnum`, `day`, `death`, `ALL_DAYS`, `new_drug_cat` and other clinical variables.
#' @return A modified `daily_data` data frame with variable `new_abx_start` added
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
#'   new_drug_cat  = c(0,1,0, 0,1, 0,1,0, 0,1,0,0,0),
#'   abx_daily      = c(0,1,1, 0,1, 0,1,1, 0,1,1,1,1)
#' )
#' define_new_abx(daily_data)
#' @export
define_new_abx <- function(daily_data) {
  # Check for required columns
  req <- c("unique_pt_id", "seqnum", "day", "abx_daily", "new_drug_cat")
  miss <- setdiff(req, names(daily_data))
  if (length(miss) > 0) {
    stop("define_new_abx(): missing required columns: ", paste(miss, collapse = ", "))
  }
  
  daily_data <- daily_data %>%
    dplyr::arrange(unique_pt_id, seqnum, day) %>%
    dplyr::group_by(unique_pt_id, seqnum) %>%
    dplyr::mutate(
      first_day_num = min(day)
    ) %>%
    dplyr::ungroup()
  
  daily_data$new_abx_start <- daily_data$new_drug_cat  # S4: a new abx drug started
  daily_data$new_abx_start[is.na(daily_data$new_drug_cat)] <- 0
  
  daily_data <- daily_data %>%
    dplyr::arrange(unique_pt_id, seqnum, day) %>%
    dplyr::group_by(unique_pt_id, seqnum) %>%
    dplyr::mutate(
      new_abx_start2 = ifelse(
        day != first_day_num &
          (dplyr::lag(abx_daily, 1) == 0 & dplyr::lag(abx_daily, 2) == 0 & abx_daily == 1),
        1, 0
      ), # S3: no abx on the prior two days
      new_abx_start3 = ifelse(
        day == (first_day_num + 1) &
          (dplyr::lag(abx_daily, 1) == 0 & abx_daily == 1),
        1, 0
      )  # S2: no abx on the first day but started on the second day
    ) %>%
    dplyr::ungroup()
  
  daily_data$new_abx_start_1 <- daily_data$new_abx_start
  daily_data$new_abx_start_2 <- daily_data$new_abx_start_1
  daily_data$new_abx_start_2[
    (daily_data$abx_daily == 1 & daily_data$day == daily_data$first_day_num) |
      daily_data$new_abx_start2 == 1] <- 1  # S1: abx started on the first day
  
  daily_data$new_abx_start_3 <- daily_data$new_abx_start_1
  daily_data$new_abx_start_3[
    (daily_data$abx_daily == 1 & daily_data$day == daily_data$first_day_num) |
      daily_data$new_abx_start2 == 1 |
      daily_data$new_abx_start3 == 1
  ] <- 1
  
  daily_data$new_abx_start <- daily_data$new_abx_start_3
  
  # Remove the intermediate variables
  daily_data$first_day_num   <- NULL
  daily_data$new_abx_start2  <- NULL
  daily_data$new_abx_start3  <- NULL
  daily_data$new_abx_start_1 <- NULL
  daily_data$new_abx_start_2 <- NULL
  daily_data$new_abx_start_3 <- NULL
  
  return(daily_data)
}

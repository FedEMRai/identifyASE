#' Define Adult Sepsis Events
#'
#' This function processes daily patient data, slices data around blood culture days, and sequentially evaluates the ASE toolkit criteria to define adult sepsis events. It allows subgroups selection and acute transfer out indication.
#'
#' @param daily_data A data frame containing daily patient data with columns `unique_pt_id`, `seqnum`, `day`, `death`, `ALL_DAYS`, and other clinical variables.
#' @param transferout_id A vector of sequence numbers (`seqnum`) indicating patients who were transferred out to an acute hospital. If no value, indicate transferout_id = NULL. 
#' @param cohort_id A list of two vectors: the first vector contains patient IDs (`unique_pt_id`) and the second vector contains sequence numbers (`seqnum`) for selecting a sub-group (default is NULL).
#' @param abx_days An integer specifying the required number of consecutive antimicrobial days (default is 4 in the absence of death, transfer or discharge, per the ASE toolkit).
#' @param window An integer specifying the number of calendar days on either side of the date of blood culture collection ie. blood culture window period. Must be an integer between 1 and 4 (default is 2, per the ASE toolkit).
#' @param creat_hi_lo_ratio The ratio of high to low creatinine levels to define renal dysfunction (default is 2, per the ASE toolkit).
#' @param creat_hi_cutoff The cutoff value for high creatinine levels (default is 44 µmol/L, per the Risk, Injury, Failure, Loss, and End-stage (RIFLE) renal disease classification system guidelines).
#' @param tbili_hi_cutoff The cutoff value for high bilirubin levels (default is 34.2 µmol/L, per the ASE toolkit).
#' @param tbili_hi_lo_ratio The ratio of high to low bilirubin levels to define liver dysfunction (default is 2, per the ASE toolkit).
#' @param lact_hi_cutoff The cutoff value for high lactate levels (default is 2 mmol/L, per the ASE toolkit).
#' @param plt_lo_cutoff The cutoff value for low platelet counts (default is 100 10^9/L, per the ASE toolkit).
#' @param plt_lo_hi_ratio The ratio of low to high platelet counts to define hematologic dysfunction (default is 0.5, per the ASE toolkit).
#' @param ed_inclusion An integer (1–4) specifying how to preprocess the ED rows in daily_data:
#'   \itemize{
#'     \item `1`: ASE Toolkit default — keep 2 days in ED
#'       (rows with `day > -2`).
#'     \item `2`: Include all ED rows (no filtering).
#'     \item `3`: Drop all ED rows (`day < 1`).
#'     \item `4`: Drop ED rows, but if any ED blood culture is detected,
#'       assign `bcx_daily = 1` on hospital day 1.
#'   }
#' @return A list of sequence numbers of ASE cases categorized by onset type and a data frame containing the data surrounding the blood culture events in the specified window, with additional variables such as indicators for qualifying antimicrobial treatments, the presence of various types of acute organ dysfunctions, and indicators for sepsis onset types.
#' @examples
#' # Get path to example CSV shipped in the package
#' data_path <- system.file(
#'   "extdata", "example_input_data", "daily_data_mock.csv",
#'   package = "identifyASE"
#' )
#' daily_data <- read.csv(data_path, stringsAsFactors = FALSE)
#' define_ase(daily_data, transferout_id = NULL, ed_inclusion = 1)
#' @import dplyr
#' @import purrr
#' @import future
#' @import furrr
#' @export
define_ase <- function(daily_data,
                       transferout_id,
                       cohort_id = NULL,
                       window = 2,
                       abx_days = 4,
                       creat_hi_lo_ratio = 2,
                       creat_hi_cutoff = 44,
                       tbili_hi_cutoff = 34.2,
                       tbili_hi_lo_ratio = 2,
                       lact_hi_cutoff = 2,
                       plt_lo_cutoff = 100,
                       plt_lo_hi_ratio = 0.5,
                       ed_inclusion = 1) {
  
  
  # apply ED data preprocessing before downstream operations
  daily_data <- process_ed_days(daily_data, ed_inclusion)
  
  # define the start of new antimicrobial administration
  daily_data <- define_new_abx(daily_data)
  
  # validate the window parameter
  if (!window %in% 1:4) {
    stop("Error: The 'window' parameter must be an integer between 1 and 4.")
  }
  
  # select a sub-cohort, e.g., ICU
  if(!is.null(cohort_id)) {
    daily_data <- daily_data[daily_data$unique_pt_id %in% unique(cohort_id[[1]]), ]
    daily_data <- daily_data[daily_data$seqnum %in% unique(cohort_id[[2]]), ] %>%
      arrange(unique_pt_id, seqnum, day)
  }
  
  # create acute transfer out indicator
  daily_data$transfer_acute <- ifelse(daily_data$seqnum %in% unique(transferout_id), 1, 0)
  daily_data$transfer_acute[daily_data$day<daily_data$ALL_DAYS] <- 0
  
  daily_data <- daily_data %>%
    group_by(unique_pt_id, seqnum) %>%
    mutate(death_day=ifelse(death==1, ALL_DAYS, NA),
           disch_day=ifelse(death==0, ALL_DAYS, NA)) %>%
    ungroup()
  
  # add bcx window +/-1 day, +/-2 days, +/-3 days, +/-4 days
  
  daily_data <- add_window_day(daily_data, window_day_col = paste0("window_day", window), window = window)
  
  # slice patient daily data around blood culture days
  
  sliced_data_list <- slice_bcx_data(daily_data,
                                     slide_day_before = window,
                                     slide_day_after = window + 4)
  
  # Apply add_window_day to each slice in the sliced_data_list; this function is slow
  
  updated_sliced_data_list <- apply_all_transformations(sliced_data_list,
                                                        window_day_col = paste0("window_day", window),
                                                        aim=1,
                                                        abx_days,
                                                        creat_hi_lo_ratio,
                                                        creat_hi_cutoff,
                                                        tbili_hi_cutoff,
                                                        tbili_hi_lo_ratio,
                                                        lact_hi_cutoff,
                                                        plt_lo_cutoff,
                                                        plt_lo_hi_ratio)
  
  
  final_combined_data <- bind_rows(lapply(updated_sliced_data_list, function(slice_info) slice_info$data))
  final_combined_data <- final_combined_data[!duplicated(final_combined_data), ]
  
  
  #### sepsis seqnum list: bcs +/-(x)days window
  sepsis_v2_com_rev_seqnum <- unique(final_combined_data$seqnum[final_combined_data$sepsis_com_v2==1])
  sepsis_v2_hosp_rev_seqnum <- unique(final_combined_data$seqnum[final_combined_data$sepsis_hosp_v2==1])
  sepsis_v2_hosp_rev_seqnum <- sepsis_v2_hosp_rev_seqnum[!(sepsis_v2_hosp_rev_seqnum %in% sepsis_v2_com_rev_seqnum)]
  sepsis_v2_rev_seqnum <- unique(c(sepsis_v2_com_rev_seqnum, sepsis_v2_hosp_rev_seqnum))
  
  
  
  ase_pid_list <- list(sepsis_v2_com_rev_seqnum, sepsis_v2_hosp_rev_seqnum, sepsis_v2_rev_seqnum)
  names(ase_pid_list) <- c("community-onset", "hospital-onset", "all_ase")
  
  result_list <- list(ase_pid_list, final_combined_data)
  dynamic_name <- paste0("Data_", window, "days")
  names(result_list) <- c("IDs", dynamic_name)
  
  return(result_list)
}


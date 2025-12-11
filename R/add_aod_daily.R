#' Add Acute Organ Dysfunction (AOD) Daily Variables
#'
#' This function calculates and adds daily indicators for various types of acute organ dysfunction (AOD) to the data,
#' based on given criteria for cardiovascular, respiratory, lactate, liver, renal, and hematologic dysfunctions.
#'
#' @param data A data frame containing patient data with columns `unique_pt_id`, `seqnum`, `day`, `vasop_daily`, `imv_daily`, `lact_daily_hi`, `tbili_daily_hi`, `tbili_baseline`, `creat_daily_hi`, `creat_baseline`, `plt_daily_lo`, `plt_baseline`, and `esrd_icd`.
#' @param window_day_col A string specifying the name of the column indicating the presence of blood culture days within a specified window.
#' @param creat_hi_lo_ratio The ratio of high to low creatinine levels to define renal dysfunction (default is 2, per the ASE toolkit).
#' @param creat_hi_cutoff The cutoff value for high creatinine levels (default is 44 µmol/L, per the Risk, Injury, Failure, Loss, and End-stage (RIFLE) renal disease classification system guidelines).
#' @param tbili_hi_cutoff The cutoff value for high bilirubin levels (default is 34.2 µmol/L, per the ASE toolkit).
#' @param tbili_hi_lo_ratio The ratio of high to low bilirubin levels to define liver dysfunction (default is 2, per the ASE toolkit).
#' @param lact_hi_cutoff The cutoff value for high lactate levels (default is 2 mmol/L, per the ASE toolkit).
#' @param plt_lo_cutoff The cutoff value for low platelet counts (default is 100 10^9/L, per the ASE toolkit).
#' @param plt_lo_hi_ratio The ratio of low to high platelet counts to define hematologic dysfunction (default is 0.5, per the ASE toolkit).
#' @return A data frame with new columns for each type of AOD and additional calculations based on the specified window.
#' @examples
#' # Example data frame
#' data <- data.frame(
#'   unique_pt_id = c(1, 1, 1, 2, 2, 2),
#'   seqnum = c(1111, 1111, 1111, 2222, 2222, 2222),
#'   day = c(1, 2, 3, 1, 2, 3),
#'   vasop_daily = c(0, 1, 0, 0, 1, 0),
#'   imv_daily = c(0, 0, 1, 1, 0, 1),
#'   lact_daily_hi = c(1.5, 2.5, 1.8, 1.9, 2.1, 3.0),
#'   tbili_daily_hi = c(30, 35, 40, 32, 36, 38),
#'   tbili_baseline = c(20, 20, 20, 25, 25, 25),
#'   creat_daily_hi = c(40, 45, 50, 35, 60, 55),
#'   creat_baseline = c(20, 20, 20, 25, 25, 25),
#'   plt_daily_lo = c(150, 80, 90, 110, 70, 50),
#'   plt_baseline = c(200, 200, 200, 180, 180, 180),
#'   esrd_icd = c(0, 0, 0, 0, 0, 0)
#' )
#' data <- add_window_day(data, "window_day", 1)
#' add_aod_daily(data, "window_day")
#' @import dplyr
#' @import purrr
#' @import future
#' @import furrr
#' @export
add_aod_daily <- function(data,
                          window_day_col,
                          creat_hi_lo_ratio = 2,
                          creat_hi_cutoff = 44,
                          tbili_hi_cutoff = 34.2,
                          tbili_hi_lo_ratio = 2,
                          lact_hi_cutoff = 2,
                          plt_lo_cutoff = 100,
                          plt_lo_hi_ratio = 0.5
                          ) {
  data <- data %>%
    arrange(unique_pt_id, seqnum, day) %>%
    group_by(unique_pt_id, seqnum) %>%
    mutate(
      aod_cv_daily = ifelse(vasop_daily == 1 & (lag(vasop_daily) == 0 | is.na(lag(vasop_daily))), 1, 0),
      aod_imv_daily = ifelse(imv_daily == 1 & (lag(imv_daily) == 0 | is.na(lag(imv_daily))), 1, 0),
      aod_lactate_daily = ifelse(lact_daily_hi >= lact_hi_cutoff, 1, 0),
      aod_liver_daily = ifelse(tbili_daily_hi >= tbili_hi_lo_ratio * tbili_baseline & tbili_daily_hi >= tbili_hi_cutoff & tbili_daily_hi >= 0 & tbili_baseline >= 0, 1, 0),
      aod_renal_daily = ifelse(creat_daily_hi >= creat_hi_lo_ratio * creat_baseline & creat_daily_hi >= creat_hi_cutoff & creat_baseline >= 0 & creat_daily_hi >= 0 & creat_daily_hi >= creat_baseline & esrd_icd != 1, 1, 0),
      aod_heme_daily = ifelse(plt_baseline >= plt_lo_cutoff & plt_baseline >= 0 & plt_daily_lo < plt_lo_cutoff & plt_daily_lo >= 0 & plt_daily_lo <= plt_lo_hi_ratio * plt_baseline & plt_daily_lo >= 0, 1, 0)
    ) %>%
    mutate(
      aod_cv_daily = replace(aod_cv_daily, is.na(aod_cv_daily), 0),
      aod_imv_daily = replace(aod_imv_daily, is.na(aod_imv_daily), 0),
      aod_lactate_daily = replace(aod_lactate_daily, is.na(aod_lactate_daily), 0),
      aod_liver_daily = replace(aod_liver_daily, is.na(aod_liver_daily), 0),
      aod_renal_daily = replace(aod_renal_daily, is.na(aod_renal_daily), 0),
      aod_heme_daily = replace(aod_heme_daily, is.na(aod_heme_daily), 0)
    ) %>%
    ungroup() %>%
    group_by(unique_pt_id, seqnum) %>%
    mutate(
      tbili_baseline_window = ifelse(any(.data[[window_day_col]] == 1) & any(!is.na(tbili_daily_lo)),
                                     ifelse(all(is.na(tbili_daily_lo[.data[[window_day_col]] == 1])), NA, min(tbili_daily_lo[.data[[window_day_col]] == 1], na.rm = TRUE)), NA),
      creat_baseline_window = ifelse(any(.data[[window_day_col]] == 1) & any(!is.na(creat_daily_lo)),
                                     ifelse(all(is.na(creat_daily_lo[.data[[window_day_col]] == 1])), NA, min(creat_daily_lo[.data[[window_day_col]] == 1], na.rm = TRUE)), NA),
      plt_baseline_window = ifelse(any(.data[[window_day_col]] == 1) & any(!is.na(plt_daily_lo)),
                                   ifelse(all(is.na(plt_daily_lo[.data[[window_day_col]] == 1])), NA, min(plt_daily_hi[.data[[window_day_col]] == 1], na.rm = TRUE)), NA)
    ) %>%
    ungroup() %>%
    mutate(
      aod_liver_daily_window = ifelse(tbili_daily_hi >= tbili_hi_lo_ratio * tbili_baseline_window & tbili_daily_hi >= tbili_hi_cutoff & tbili_daily_hi >= 0 & tbili_baseline_window >= 0, 1, 0),
      aod_renal_daily_window = ifelse(creat_daily_hi >= creat_hi_lo_ratio * creat_baseline_window & creat_daily_hi >= 0 & creat_baseline_window >= 0 & creat_daily_hi >= creat_hi_cutoff & creat_daily_hi >= creat_baseline_window & esrd_icd != 1, 1, 0),
      aod_heme_daily_window = ifelse(plt_baseline_window >= plt_lo_cutoff & plt_baseline_window >= 0 & plt_daily_lo < plt_lo_cutoff & plt_daily_lo >= 0 & plt_daily_lo <= plt_lo_hi_ratio * plt_baseline_window & plt_daily_lo >= 0, 1, 0),
      aod_liver_daily_window = replace(aod_liver_daily_window, is.na(aod_liver_daily_window), 0),
      aod_renal_daily_window = replace(aod_renal_daily_window, is.na(aod_renal_daily_window), 0),
      aod_heme_daily_window = replace(aod_heme_daily_window, is.na(aod_heme_daily_window), 0)
    ) %>%
    mutate(
      aod_any_daily_comm = ifelse(((aod_renal_daily == 1 & .data[[window_day_col]] == 1) |
                                     (aod_liver_daily == 1 & .data[[window_day_col]] == 1) |
                                     (aod_heme_daily == 1 & .data[[window_day_col]] == 1) |
                                     (aod_cv_daily == 1 & .data[[window_day_col]] == 1) |
                                     (aod_imv_daily == 1 & .data[[window_day_col]] == 1) |
                                     (aod_lactate_daily == 1 & .data[[window_day_col]] == 1)), 1, 0),

      aod_any_daily_hosp = ifelse(((aod_renal_daily_window == 1 & .data[[window_day_col]] == 1) |
                                     (aod_liver_daily_window == 1 & .data[[window_day_col]] == 1) |
                                     (aod_heme_daily_window == 1 & .data[[window_day_col]] == 1) |
                                     (aod_cv_daily == 1 & .data[[window_day_col]] == 1) |
                                     (aod_imv_daily == 1 & .data[[window_day_col]] == 1) |
                                     (aod_lactate_daily == 1 & .data[[window_day_col]] == 1)), 1, 0)
    )

  return(data)
}

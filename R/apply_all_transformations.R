#' Apply All Transformations to Each Slice in the List
#'
#' This function applies a series of transformations to each slice in the list, including adding AOD daily variables, identifying qualifying antimicrobial treatments, and defining sepsis.
#'
#' @param sliced_data_list A list of lists, where each list contains `unique_pt_id`, `seqnum`, and a data frame (`data`) with sliced data.
#' @param window_day_col A string specifying the name of the column indicating the presence of blood culture days within a specified window.
#' @param aim An integer specifying the aim criteria for qualifying antimicrobial treatments: 1 for original criteria or 2 for extended criteria.
#' @param abx_days An integer specifying the number of consecutive antimicrobial days.
#' @param creat_hi_lo_ratio The ratio of high to low creatinine levels to define renal dysfunction (default is 2, as per the ASE toolkit).
#' @param creat_hi_cutoff The cutoff value for high creatinine levels (default is 44 µmol/L, as per the Risk, Injury, Failure, Loss, and End-stage (RIFLE) renal disease classification system guidelines).
#' @param tbili_hi_cutoff The cutoff value for high bilirubin levels (default is 34.2 µmol/L, as per the ASE toolkit).
#' @param tbili_hi_lo_ratio The ratio of high to low bilirubin levels to define liver dysfunction (default is 2, as per the ASE toolkit).
#' @param lact_hi_cutoff The cutoff value for high lactate levels (default is 2 mmol/L, as per the ASE toolkit).
#' @param plt_lo_cutoff The cutoff value for low platelet counts (default is 100 10^9/L, as per the ASE toolkit).
#' @param plt_lo_hi_ratio The ratio of low to high platelet counts to define hematologic dysfunction (default is 0.5, as per the ASE toolkit).
#' @return A list of updated slices with the transformations applied.
#' @examples
#' # Example sliced_data_list
#' sliced_data_list <- list(
#'   list(unique_pt_id = 1, seqnum = 1111, data = data.frame(day = 1:3, bcx_daily = c(0, 1, 0), vasop_daily = c(0, 1, 0), imv_daily = c(0, 0, 1), lact_daily_hi = c(1.5, 2.5, 1.8), tbili_daily_hi = c(30, 35, 40), tbili_baseline = c(20, 20, 20), creat_daily_hi = c(40, 45, 50), creat_baseline = c(20, 20, 20), plt_daily_lo = c(150, 80, 90), plt_baseline = c(200, 200, 200), esrd_icd = c(0, 0, 0), new_abx_start = c(0, 1, 0), abx_daily = c(0, 1, 1), death = c(0, 0, 0), transfer_acute = c(0, 0, 0), ALL_DAYS = c(3, 3, 3), window_day = c(1, 1, 0))),
#'   list(unique_pt_id = 2, seqnum = 2222, data = data.frame(day = 1:3, bcx_daily = c(0, 0, 1), vasop_daily = c(0, 0, 0), imv_daily = c(0, 1, 0), lact_daily_hi = c(1.9, 2.1, 3.0), tbili_daily_hi = c(32, 36, 38), tbili_baseline = c(25, 25, 25), creat_daily_hi = c(35, 60, 55), creat_baseline = c(25, 25, 25), plt_daily_lo = c(110, 70, 50), plt_baseline = c(180, 180, 180), esrd_icd = c(0, 0, 0), new_abx_start = c(0, 1, 0), abx_daily = c(0, 1, 1), death = c(0, 1, 0), transfer_acute = c(0, 0, 1), ALL_DAYS = c(3, 3, 3), window_day = c(1, 1, 0)))
#' )
#' apply_all_transformations(sliced_data_list, "window_day", 1)
#' @import dplyr
#' @import purrr
#' @import future
#' @import furrr
#' @export
apply_all_transformations <- function(sliced_data_list,
                                      window_day_col,
                                      aim,
                                      abx_days = 4,
                                      creat_hi_lo_ratio = 2,
                                      creat_hi_cutoff = 44,
                                      tbili_hi_cutoff = 34.2,
                                      tbili_hi_lo_ratio = 2,
                                      lact_hi_cutoff = 2,
                                      plt_lo_cutoff = 100,
                                      plt_lo_hi_ratio = 0.5) {
  plan(multisession)  # Set up parallel processing

  updated_sliced_data_list <- future_map(sliced_data_list, function(slice_info) {
    slice_info$data <- add_aod_daily(slice_info$data,
                                     window_day_col,
                                     creat_hi_lo_ratio,
                                     creat_hi_cutoff,
                                     tbili_hi_cutoff,
                                     tbili_hi_lo_ratio,
                                     lact_hi_cutoff,
                                     plt_lo_cutoff,
                                     plt_lo_hi_ratio)
    slice_info$data <- qualifying_abx_duration(slice_info$data, window_day_col, aim, abx_days)
    slice_info$data <- define_sepsis(slice_info$data)
    return(slice_info)
  }, .options = furrr_options(seed = TRUE))

  return(updated_sliced_data_list)
}

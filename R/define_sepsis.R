#' Define Sepsis
#'
#' This function defines sepsis based on the onset of acute organ dysfunction (AOD) and qualifying antimicrobial treatments. It creates indicators for community-acquired and hospital-acquired sepsis.
#'
#' @param data A data frame containing patient data with columns `unique_pt_id`, `seqnum`, `day`, `aod_any_daily_comm`, `aod_any_daily_hosp`, `bcx_daily`, `abx_window_startday`, and `abx_qualifying_ep`.
#' @return A data frame with new columns indicating the onset days for community and hospital sepsis and indicators for sepsis onset types.
#' @examples
#' # Example data frame
#' data <- data.frame(
#'   unique_pt_id = c(1, 1, 1, 2, 2, 2),
#'   seqnum = c(1111, 1111, 1111, 2222, 2222, 2222),
#'   day = c(1, 2, 3, 1, 2, 3),
#'   aod_any_daily_comm = c(0, 1, 0, 0, 0, 1),
#'   aod_any_daily_hosp = c(0, 0, 1, 0, 1, 0),
#'   bcx_daily = c(1, 0, 0, 0, 1, 0),
#'   abx_window_startday = c(1, 1, 1, 2, 2, 2),
#'   abx_qualifying_ep = c(1, 1, 1, 0, 1, 1)
#' )
#' define_sepsis(data)
#' @import dplyr
#' @import purrr
#' @import future
#' @import furrr
#' @export
define_sepsis <- function(data) {
  data <- data %>%
    arrange(unique_pt_id, seqnum, day) %>%
    group_by(unique_pt_id, seqnum) %>%
    mutate(
      first_aod_any_daily_comm = ifelse(any(aod_any_daily_comm==1), min(day[aod_any_daily_comm==1]), NA),
      first_aod_any_daily_hosp = ifelse(any(aod_any_daily_hosp==1), min(day[aod_any_daily_hosp==1]), NA),
      first_bcx_day = ifelse(any(bcx_daily==1), min(day[bcx_daily==1]), NA),
      onset_day_comm = pmin(abx_window_startday, first_bcx_day, first_aod_any_daily_comm, na.rm = TRUE),
      onset_day_hosp = pmin(abx_window_startday, first_bcx_day, first_aod_any_daily_hosp, na.rm = TRUE),
      sepsis_com_v2 = ifelse(onset_day_comm <=2 & (abx_qualifying_ep==1) & (aod_any_daily_comm==1), 1, 0),
      sepsis_hosp_v2 = ifelse(onset_day_hosp>2 & (abx_qualifying_ep==1) & (aod_any_daily_hosp==1), 1, 0)
    ) %>%
    ungroup()

  return(data)
}

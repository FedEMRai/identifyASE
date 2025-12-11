#' Add Window Day Variable
#'
#' This function adds a new column to the data indicating whether there is any blood culture day within a specified window around each day.
#'
#' @param data A data frame containing patient data with columns `unique_pt_id`, `seqnum`, `day`, and `bcx_daily`.
#' @param window_day_col A string specifying the name of the new column to be added.
#' @param window An integer specifying the number of days before and after each day to check for blood culture days.
#' @return A data frame with the new column indicating the presence of blood culture days within the specified window.
#' @examples
#' data <- data.frame(
#'   unique_pt_id = c(1, 1, 1, 2, 2, 2),
#'   seqnum = c(1111, 1111, 1111, 2222, 2222, 2222),
#'   day = c(1, 2, 3, 1, 2, 3),
#'   bcx_daily = c(0, 1, 0, 0, 1, 0)
#' )
#' add_window_day(data, "window_day", 1)
#' @import dplyr
#' @import purrr
#' @import future
#' @import furrr
#' @export
add_window_day <- function(data, window_day_col, window) {
  data <- data %>%
    arrange(unique_pt_id, seqnum, day) %>%
    group_by(unique_pt_id, seqnum) %>%
    mutate(!!window_day_col := map_dbl(day, ~ {
      if (any(bcx_daily[day >= (.x - window) & day <= (.x + window)] == 1, na.rm = TRUE)) {
        return(1)
      } else {
        return(0)
      }
    })) %>%
    ungroup()

  return(data)
}

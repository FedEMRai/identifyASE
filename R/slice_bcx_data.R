# Load necessary libraries
pack <- c("dplyr", "purrr", "future", "furrr")

new.packages <- pack[!(pack %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pack, require, character.only = TRUE)
options(scipen = 999)

#' Slice Data Around Blood Culture Collections
#'
#' This function slices daily data around blood culture collections for each patient, creating a list of data frames
#' containing the specified window of days before and after each blood culture collection.
#'
#' @param data A data frame containing patient data with columns `unique_pt_id`, `seqnum`, `day`, and `bcx_daily`.
#' @param slide_day_before An integer specifying the number of days before each blood culture collection to include in the slice (default is 2, per the ASE toolkit).
#' @param slide_day_after An integer specifying the number of days after each blood culture collection to include in the slice (default is 6 to track sequential antimicrobial usage if a new antimicrobial starts on the final BCWP day).
#' @return A list of lists, where each list contains `unique_pt_id`, `seqnum`, and a data frame (`data`) with the sliced data.
#' @examples
#' data <- data.frame(
#'   unique_pt_id = c(1, 1, 1, 2, 2, 2),
#'   seqnum = c(1111, 1111, 1111, 2222, 2222, 2222),
#'   day = c(1, 2, 3, 1, 2, 3),
#'   bcx_daily = c(0, 1, 0, 0, 1, 0)
#' )
#' slice_bcx_data(data)
#' @import dplyr
#' @import purrr
#' @import future
#' @import furrr
#' @export
slice_bcx_data <- function(data,
                           slide_day_before=2,  # this equals to # of window days
                           slide_day_after=6    # this equals to # of window days+4: 6 or 7
) {

  # Filter out patients without any blood culture days
  filtered_data <- data %>%
    group_by(unique_pt_id, seqnum) %>%
    filter(any(bcx_daily == 1)) %>%
    ungroup() %>%
    arrange(unique_pt_id, seqnum, day)

  # Split the data by patient
  patient_data_list <- filtered_data %>%
    group_split(unique_pt_id, seqnum)

  # Initialize an empty list to store the slices with IDs
  sliced_data_list <- list()

  # Process each patient's data
  for (patient_data in patient_data_list) {
    unique_pt_id <- unique(patient_data$unique_pt_id)
    seqnum <- unique(patient_data$seqnum)
    bcx_days <- which(patient_data$bcx_daily == 1)
    for (bcx_day in bcx_days) {
      start_day <- max(1, bcx_day - slide_day_before)
      end_day <- min(nrow(patient_data), bcx_day + slide_day_after)
      slice <- patient_data[start_day:end_day, ]
      slice_info <- list(unique_pt_id = unique_pt_id, seqnum = seqnum, data = slice)
      sliced_data_list <- append(sliced_data_list, list(slice_info))
    }
  }

  return(sliced_data_list)
}

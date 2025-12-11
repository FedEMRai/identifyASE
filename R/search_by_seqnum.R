# Load necessary libraries
pack <- c("dplyr", "purrr", "future", "furrr")

new.packages <- pack[!(pack %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pack, require, character.only = TRUE)
options(scipen = 999)

#' Search for Data by Sequence Number
#'
#' This function searches through a list of sliced data and returns the data slices that match a given sequence number.
#'
#' @param sliced_data_list A list of lists, where each list contains `unique_pt_id`, `seqnum`, and a data frame (`data`) with sliced data.
#' @param target_seqnum An integer specifying the sequence number to search for.
#' @return A list of data frames that match the specified sequence number. If no matches are found, an empty list is returned.
#' @examples
#' # Example sliced_data_list
#' sliced_data_list <- list(
#'   list(unique_pt_id = 1, seqnum = 1111, data = data.frame(day = 1:3, bcx_daily = c(0, 1, 0))),
#'   list(unique_pt_id = 2, seqnum = 2222, data = data.frame(day = 1:3, bcx_daily = c(0, 1, 0)))
#' )
#' search_by_seqnum(sliced_data_list, 1111)
#' @import dplyr
#' @import purrr
#' @import future
#' @import furrr
#' @export
search_by_seqnum <- function(sliced_data_list, target_seqnum) {
  results <- lapply(sliced_data_list, function(slice_info) {
    if (slice_info$seqnum == target_seqnum) {
      return(slice_info$data)
    } else {
      return(NULL)
    }
  })

  # Filter out NULL values from the results
  results <- results[!sapply(results, is.null)]

  return(results)
}

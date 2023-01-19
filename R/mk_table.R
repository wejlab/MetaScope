# Helper function to create a taxonomy table
mk_table <- function(intable, taxon_ranks) {
  this_t <- as.data.frame(intable)
  if (!identical(this_t$rank, character(0))) {
    t_n <- nrow(this_t)
    if (this_t$rank[t_n] == "no rank" && this_t$rank[t_n - 1] == "species") {
      this_t$rank[t_n] <- "strain"
    }
    final_tab <- this_t %>%
      dplyr::filter(.data$rank %in% taxon_ranks) %>%
      dplyr::right_join(dplyr::tibble(rank = taxon_ranks),
                        by = "rank") %>%
      dplyr::arrange(factor(.data$rank, levels = taxon_ranks)) %>%
      dplyr::select(.data$name)
    final_tab[seq_along(taxon_ranks), ] %>%
      magrittr::set_names(taxon_ranks) %>% return()
  }
}

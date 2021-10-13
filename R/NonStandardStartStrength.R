
experimentally_determined_start_strengths <- function(){
  path=system.file("extdata/nonstandard_start_tis_efficiency.tsv", package = "startstrength")
  read.csv(path, sep = "\t", header = TRUE) %>%
    dplyr::tibble() %>%
    dplyr::mutate(StartCodon=substr(TIS.Sequence, 5, 7))
}

plot_distributions_experimentally_determined_start_strengths <- function(){
  experimentally_determined_start_strengths() %>%
    ggplot2::ggplot(ggplot2::aes(TIS.Efficiency)) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(~ StartCodon, scales = "free_y") +
    ggplot2::scale_x_log10()
}

#Input sequence must be AAAA[UUG]A


#' Start Strength
#'
#' Use sequence context to identify the EXPERIMENTALLY determined tranlsation initiation start efficiency
#' Built using data from DOI: 10.1093/nar/gkx1114 . Please cite if you find this tool useful
#'
#' @param input_sequence  sequences of non-standard start sites to test. All input sequences must be 8bp long with start codon in pos 5bp-7bp, for example AAAA(UUG)A. No brackets should be included. (character)
#'
#' @return Translation initiation start efficiency of each input sequence (dataframe)
#' @export
#'
#' @examples
#' start_strength_nonstandard("TACTTTGG")
start_strength_nonstandard <- function(input_sequence){
  assertthat::assert_that(all(nchar(input_sequence) == 8), msg = "Please make sure all input sequences are 8bp long with start codon in pos 5bp-7bp, for example AAAA[UUG]A (no brackets should be included)")

  #uppercase sequence and turn any T's to U's
  input_sequence <- toupper(gsub(pattern = "T", replacement = "U", x = input_sequence))

  assertthat::assert_that(all((strsplit(input_sequence, split = "") %>% unlist() %>% unique()) %in% c("A","U","C","G")), msg = "Found non A,U,C,G,T character. Please remove then run again")

  head(input_sequence, n=3) %>%
    sapply(FUN = function(x) {
      message("Start codons assumed to be at pos 5-7bp (indicated in square brackets): ")
      message(substr(x, 1, 4), "[", substr(x, 5, 7),"]", substr(x, 8, 9))
      })
  message("\n")
  start_codons = substr(input_sequence, 5, 7)

  start_strengths_database_df <- experimentally_determined_start_strengths()
  supported_start_codons <- unique(start_strengths_database_df$StartCodon)


  assertthat::assert_that(all(start_codons %in% supported_start_codons), msg = paste0("Unsupported start codon. We support only ", paste0(supported_start_codons, collapse = ", ")))

  dplyr::tibble(
    StartCodon=start_codons,
    Sequence=input_sequence,
    TISefficiency=start_strengths_database_df$TIS.Efficiency[match(input_sequence, start_strengths_database_df$TIS.Sequence)]
    )

}

#' Start Strength
#'
#' Use sequence context to identify the EXPERIMENTALLY determined tranlsation initiation start efficiency
#' Built using data from DOI: 10.1093/nar/gkx1114 . Please cite if you find this tool useful
#'
#' @param input_sequence  sequences of non-standard start sites to test. All input sequences must be 8bp long with start codon in pos 5bp-7bp, for example AAAA(UUG)A. No brackets should be included (character)
#'
#' @return Plot showing how Translation initiation efficiency in given contexts relates to the distribution of all possible sequence contexts
#' @export
#'
#' @examples
#' plot_start_start_strength_nonstandard("TACTTTGG")
plot_start_start_strength_nonstandard <- function(input_sequence) {
  start_strengths_database_df <- experimentally_determined_start_strengths()

  results <- start_strength_nonstandard(input_sequence)
  print(results)

  relevant_start_codons <- unique(results$StartCodon)


  start_strengths_database_df %>%
    dplyr::filter(StartCodon %in% relevant_start_codons) %>%
    ggplot2::ggplot(ggplot2::aes(TIS.Efficiency)) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(~ StartCodon, scales = "free_y") +
    ggplot2::geom_vline(data = results, ggplot2::aes(xintercept=TISefficiency, color=Sequence)) +
    ggplot2::scale_x_log10()
}

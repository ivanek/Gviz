#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param use.defaults `logical`. Derive the subsetting range from the track's
#' own defaults rather than using `from` and `to` verbatim.
#' @param filter A named `list` of additional filters to be applied in the
#' Biomart query.
#' @param range A [`GRanges`][GenomicRanges::GRanges-class] object with the
#' genomic coordinates of the gene model items.
#' @param from,to Numeric scalar, giving the range of genomic coordinates to
#' limit the tracks in. Note that `from` cannot be larger than `to`.
#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.
#'
#' @param start An integer scalar with the genomic start coordinate for the
#' gene model range.
#' @param end An integer scalar with the genomic end coordinate for the gene
#' model range.
#' @param biomart An optional [`Mart`][biomaRt::Mart-class] object providing
#' access to the EBI Biomart webservice. By default the appropriate Ensembl
#' data source is selected based on the provided genome and chromosome.
#' @param strand Character scalar, the strand for which to fetch gene
#' information from Biomart. One of `+`, `-`, or `+-`.
#' @param chromosome The chromosome on which the track's genomic ranges are
#' defined. A valid UCSC chromosome identifier. Please note that at this stage
#' only syntactic checking takes place, i.e., the argument value needs to be a
#' single integer, numeric character or a character of the form `chrx`, where
#' `x` may be any possible string. The user has to make sure that the respective
#' chromosome is indeed defined for the track's genome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point. If no mapping from genome to Biomart Ensembl data
#' source is possible, the `biomart` argument needs to be provided by the user.
#' @param stacking The stacking type for overlapping items of the track. One of
#' `c(hide, dense, squish, pack, full)`. Currently, only hide (do not show the
#' track items), squish (make best use of the available space) and dense (no
#' stacking at all) are implemented.
#' @param filters A list of additional filters to be applied in the Biomart
#' query. See [biomaRt::getBM()] for details.
#' @param featureMap Named character vector or list to map between the fields
#' in the Biomart data base and the features as they are used to construct the
#' track. If multiple values are provided in a single list item, the package
#' will use the first one that is defined in the selected Biomart.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param symbol,transcript,gene,entrez Character vector giving one or several
#' gene symbols, Ensembl transcript identifiers, Ensembl gene identifiers, or
#' ENTREZ gene identifiers, respectively. The genomic locus of their gene model
#' will be fetched from Biomart instead of providing explicit start and end
#' coordinates.
#' @param ... Additional items which will all be interpreted as further
#' display parameters. See [`settings`] and the "Display Parameters"
#' section below for details.

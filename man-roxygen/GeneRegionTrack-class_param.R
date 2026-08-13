#' @param GdObject Object of class [`GdObject`][GdObject-class].
#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param value Value to be set.
#' @param defaults A `list` of the constructor's default arguments, used to
#' fill in whatever has not been provided in `args`.
#' @param object Object inheriting from
#' [`GeneRegionTrack`][GeneRegionTrack-class].
#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.
#' @param reference A `character` scalar with the path to the referenced file.
#' @param mapping A named `list` mapping the columns of the imported data to
#' the metadata columns of the track's [`GRanges`][GenomicRanges::GRanges-class]
#' object.
#' @param args A `list` of the arguments the object has been constructed with.
#' @param type A `character` scalar with the name of the track class, used to
#' label the streamed reference in the object's `show` method.
#'
#' @param start,end An integer scalar with the genomic start or end coordinate
#' for the gene model range. If those are missing, the default value will
#' automatically be the smallest (or largest) value, respectively, in `rstarts`
#' and `rends` for the currently active chromosome. When building a
#' [`GeneRegionTrack`][GeneRegionTrack-class] from a
#' [`TxDb`][GenomicFeatures::TxDb-class] object, these arguments can be used to
#' subset the desired annotation data by genomic coordinates. Please note that
#' in this case the `chromosome` parameter must also be set.
#' @param rstarts An integer vector of the start coordinates for the actual
#' gene model items, i.e., for the individual exons. The relationship between
#' exons is handled via the `gene` and `transcript` factors. Alternatively,
#' this can be a vector of comma-separated lists of integer coordinates, one
#' vector item for each transcript, and each comma-separated element being the
#' start location of a single exon within that transcript. Those lists will be
#' exploded upon object instantiation and all other annotation arguments will be
#' recycled accordingly to regenerate the exon/transcript/gene relationship
#' structure. This implies the appropriate number of items in all annotation and
#' coordinate arguments.
#' @param rends An integer vector of the end coordinates for the actual gene
#' model items. Both `rstarts` and `rends` have to be of equal length.
#' @param rwidths An integer vector of widths for the actual gene model items.
#' This can be used instead of either `rstarts` or `rends` to specify the range
#' coordinates.
#' @param feature Factor (or other vector that can be coerced into one), giving
#' the feature types for the individual track exons. When plotting the track
#' to the device, if a display parameter with the same name as the value of
#' `feature` is set, this will be used as the track item's fill color.
#' Additionally, the feature type defines whether an element in the
#' [`GeneRegionTrack`][GeneRegionTrack-class] is considered to be coding or
#' non-coding. The 'Details' section as well as the section about the
#' `thinBoxFeature` display parameter further below has more information on
#' this. See also [`grouping`] for details.
#' @param exon Character vector of exon identifiers. Its values will be used
#' as the identifier tag when plotting to the device if the display parameter
#' `showExonId=TRUE`.
#' @param strand Character vector, the strand information for the individual
#' track exons. It may be provided in the form `+` for the Watson strand, `-`
#' for the Crick strand or `*` for either one of the two. Please note that all
#' items within a single gene or transcript model need to be on the same strand,
#' and erroneous entries will result in an error.
#' @param transcript Factor (or other vector that can be coerced into one),
#' giving the transcript memberships for the individual track exons. All items
#' with the same transcript identifier will be visually connected when plotting
#' to the device. See [`grouping`] for details. Will be used as labels when
#' `showId=TRUE` and `geneSymbol=FALSE`.
#' @param gene Factor (or other vector that can be coerced into one), giving
#' the gene memberships for the individual track exons.
#' @param symbol A factor with human-readable gene name aliases which will be
#' used as labels when `showId=TRUE` and `geneSymbol=TRUE`.
#' @param chromosome The chromosome on which the track's genomic ranges are
#' defined. A valid UCSC chromosome identifier if
#' `options(ucscChromosomeNames=TRUE)`. Please note that in this case only
#' syntactic checking takes place, i.e., the argument value needs to be an
#' integer, numeric character or a character of the form `chrx`, where `x` may
#' be any possible string. The user has to make sure that the respective
#' chromosome is indeed defined for the track's genome. If not provided here,
#' the constructor will try to build the chromosome information based on the
#' available inputs, and as a last resort will fall back to the value `chrNA`.
#' Please note that by definition all objects in the `Gviz` package can only
#' have a single active chromosome at a time (although internally the
#' information for more than one chromosome may be present), and the user has to
#' call the `chromosome<-` replacement method in order to change to a different
#' active chromosome. When creating a
#' [`GeneRegionTrack`][GeneRegionTrack-class] from a
#' [`TxDb`][GenomicFeatures::TxDb-class] object, the value of this parameter can
#' be used to subset the data to fetch only transcripts from a single
#' chromosome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point. If not provided here, the constructor will try to
#' extract this information from the provided inputs, and eventually will fall
#' back to the default value of `NA`.
#' @param stacking The stacking type for overlapping items of the track. One of
#' `c(hide, dense, squish, pack, full)`. Currently, only hide (do not show the
#' track items), squish (make best use of the available space) and dense (no
#' stacking at all) are implemented.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param importFunction A user-defined function to be used to import the data
#' from a file. This only applies when the `range` argument is a character
#' string with the path to the input data file. The function needs to accept an
#' argument `x` containing the file path and has to return a proper
#' [`GRanges`][GenomicRanges::GRanges-class] object with all the necessary
#' metadata columns set. A set of default import functions is already
#' implemented in the package for a number of different file types, and one of
#' these defaults will be picked automatically based on the extension of the
#' input file name. If the extension cannot be mapped to any of the existing
#' import functions, an error is raised asking for a user-defined import
#' function via this argument. Currently the following file types can be
#' imported with the default functions: `gff`, `gff1`, `gff2`, `gff3`, `gtf`.
#' @param stream A logical flag indicating that the user-provided import
#' function can deal with indexed files and knows how to process the additional
#' `selection` argument when accessing the data on disk. This causes the
#' constructor to return a
#' [`ReferenceGeneRegionTrack`][ReferenceGeneRegionTrack-class] object which
#' will grab the necessary data on the fly during each plotting operation.
#' @param ... Additional items which will all be interpreted as further
#' display parameters. See [`settings`] and the "Display Parameters"
#' section below for details.

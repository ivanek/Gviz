#' @slot dp Object of `DisplayPars-class`, the display settings controlling the
#' look and feel of a track. See settings for details on setting graphical
#' parameters for tracks.
#' @slot name Object of class `character`, a human-readable name for the track
#' that will be used in the track's annotation panel if necessary.
#' @slot imageMap Object of `ImageMap-class`, containing optional information
#' for an HTML image map. This will be created by the `drawGD` methods when the
#' track is plotted to a device and is usually not set by the user.

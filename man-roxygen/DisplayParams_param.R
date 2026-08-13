#' @param class A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param name Name of the retrieved parameter.
#' @param value Value to be set.
#' @param object Object of class [`DisplayPars`][DisplayPars-class].
#' @param interactive `logical`. Emit the message explaining that `setPar` no
#' longer supports pass-by-reference semantics and that its result has to be
#' reassigned.
#' @param asIs `logical`. Return the queried parameters as a list. When
#' `FALSE`, the default, the result of a single-parameter query is unlisted
#' for convenience.
#' @param hideInternal `logical`. Omit the internal display parameters, i.e.,
#' those whose names are prefixed with `.__`.
#' @param recursive `logical`. For composite tracks, also set the display
#' parameters on each of the contained sub-tracks.
#' @param ... All named arguments are stored in the object's environment as
#' individual parameters, regardless of their type.

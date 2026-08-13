#' @param GdObject Object of class [`GdObject`][GdObject-class].
#' @param name Name of the retrieved parameter.
#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param value Value to be set.
#' @param recursive `logical`. For composite tracks, also set the display
#' parameters on each of the contained sub-tracks.
#' @param object Object inheriting from
#' [`ReferenceTrack`][ReferenceTrack-class].
#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.

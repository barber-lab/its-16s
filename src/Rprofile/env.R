#!/usr/bin/env R

#
# .Rprofile for interaktive analysis with make env
#

message("Calling .Rprofile for interactive analysis")

# Isolate new libraries
.libPaths(c("/analysis/R-site-library", .libPaths()))


# Auto start project
setHook("rstudio.sessionInit", function(newSession) {
  if (newSession && is.null(rstudioapi::getActiveProject())) {
    rstudioapi::openProject(".")
  }
}, action = "append")

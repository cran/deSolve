## internal helper functions for printing solver return code messages
## these functions are not exported

## print combined messages (message and numeric output)
printmessage <-function(message1, state, message2 = NULL) {
  if (is.null(message1)) {
    cat("\n", paste(formatC(1:length(message1), "##", width = 2), message1,
              signif(state, digits = getOption("digits")), "\n"), "\n")
  } else {
    cat("\n", paste(formatC(1:length(message1), "##", width = 2), message1,
              signif(state, digits = getOption("digits")), message2, "\n"), "\n")
  
  }
}

## print short messages
printM <- function(message) cat(message, "\n")



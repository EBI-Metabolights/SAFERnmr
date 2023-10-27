# thanks, GPT!

nestedListToHTML <- function(lst, indent = 0, newline = TRUE) {
  html_string <- ""
  for (key in names(lst)) {
    html_string <- paste0(html_string, strrep(" ", indent), "<details>\n")
    html_string <- paste0(html_string, strrep(" ", indent), "<summary>", key, "</summary>\n")
    html_string <- paste0(html_string, strrep(" ", indent), "<blockquote>\n")

    if (is.list(lst[[key]])) {
      html_string <- paste0(html_string, nestedListToHTML(lst[[key]], indent+1, newline))
    } else {
      for (element in lst[[key]]) {
        html_string <- paste0(html_string, strrep(" ", indent+1), element, "<br>\n")
      }
    }
    html_string <- paste0(html_string, strrep(" ", indent), "</blockquote>\n")
    html_string <- paste0(html_string, strrep(" ", indent), "</details>\n")
    if (newline) {
      html_string <- paste0(html_string, "\n")
    }
  }
  return(html_string)
}

nestedListToHTML <- function(lst, indent = 0, newline = TRUE) {
  html_string <- ""
  for (key in names(lst)) {
    html_string <- paste0(html_string, strrep(" ", indent), "<details>\n")
    html_string <- paste0(html_string, strrep(" ", indent), "<summary>", key, "</summary>\n")
    html_string <- paste0(html_string, strrep(" ", indent), "<blockquote>\n")

    if (is.list(lst[[key]])) {
      html_string <- paste0(html_string, nestedListToHTML(lst[[key]], indent + 1, newline))
    } else {
      for (element in lst[[key]]) {
        html_string <- paste0(html_string, strrep(" ", indent + 1), element, "<br>\n")
      }
    }

    html_string <- paste0(html_string, strrep(" ", indent + 1), "<details>\n")
    html_string <- paste0(html_string, strrep(" ", indent + 1), "<summary>description</summary>\n")
    html_string <- paste0(html_string, strrep(" ", indent + 1), "<blockquote>\n")
    html_string <- paste0(html_string, strrep(" ", indent + 2), "description<br>\n")
    html_string <- paste0(html_string, strrep(" ", indent + 1), "</blockquote>\n")
    html_string <- paste0(html_string, strrep(" ", indent + 1), "</details>\n")

    html_string <- paste0(html_string, strrep(" ", indent), "</blockquote>\n")
    html_string <- paste0(html_string, strrep(" ", indent), "</details>\n")
    if (newline) {
      html_string <- paste0(html_string, "\n")
    }
  }
  return(html_string)
}



# txt <- nestedListToHTML(pars)
# Convert the nested list to HTML with proper indentation and newline characters
html_output <- nestedListToHTML(pars, newline = TRUE)

# Writing the HTML output to a file
writeLines(html_output, "/Users/mjudge/Documents/GitHub/SAFERnmr/pars.html")

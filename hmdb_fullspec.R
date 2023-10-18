library(XML)
library(base64enc)
library(xml2)
library(rvest)
hmdb_fullspec <- function(nmrml.file){

  nmrml.file <- '/Users/mjudge/Downloads/500_1H_for_JSviewer_HMDB0000159.nmrML'
  f <- read_html(nmrml.file)
  
  spec <- html_element(f, xpath = '//spectrumDataArray')
  
  p <- htmlTreeParse(f)
  ptree <- p$children$html
  
  
  lines <- readLines(nmrml.file)
    has.data <- grepl(pattern = 'spectrumDataArray', x = lines, perl = TRUE) %>% which
    data.line <- lines[has.data] %>% .[1]
    <spectrumDataArray *[a-zA-Z]*=\\\"[a-zA-Z0-9]*\\\" *"
    a <- 'abcdef'

    # sda <- regexpr(data.line, pattern = '<spectrumDataArray.*</spectrumDataArray>')
    sda <- regexpr(data.line, pattern = '<spectrumDataArray')
    sda <- regexpr(data.line, pattern = '<spectrumDataArray')
    sda <- regexpr(data.line, pattern = '</spectrumDataArray>')
    sda <- regexpr(data.line, pattern = '<spectrumDataArray')
    
    a <- stringr::str_extract(data.line, pattern = '<.*>')

  fullspec <- list(data = NA,
                   ppm = NA)
  
  return(fullspec)
}


}


# Install and load the xml2 package if you haven't already
install.packages("xml2")
library(xml2)

# Load the XML document
xml_file <- "/Users/mjudge/Downloads/500_1H_for_JSviewer_HMDB0000159.nmrML"  # Replace with the path to your NMRML file
doc <- read_xml(fixed.file)

# Define the XML namespaces used in the NMRML document
ns <- c(nmrml = "http://nmrml.org/schema")

# Extract all spectrumDataArray elements
spectrum_data_elements <- xml_find_all(doc, "//nmrml:spectrumDataArray", namespaces = ns)

# Iterate through each spectrumDataArray and extract binary data
for (element in spectrum_data_elements) {
  # Extract the text content (base64-encoded binary data)
  base64_data <- xml_text(element)
  
  # Decode the base64 data to binary format
  binary_data <- base64decode::base64decode(base64_data)
  
  # Now, 'binary_data' contains the binary spectral data
  # You can process or save it as needed
  # For example, you can write it to a binary file
  binary_file_path <- paste0("spectrum_", xml_attr(element, "id"), ".bin")
  writeBin(binary_data, binary_file_path)
  
  # Print some information for reference
  cat("Spectrum ID:", xml_attr(element, "id"), "\n")
  cat("Encoded Length:", xml_attr(element, "encodedLength"), "\n")
  cat("Byte Format:", xml_attr(element, "byteFormat"), "\n")
  cat("Binary Data extracted and saved to", binary_file_path, "\n\n")
}
Make sure to replace "your_file.xml" with the actual path to your NMRML file. This script will extract the binary spectral data from each <spectrumDataArray> element, decode it from base64, and save it to binary files with informative names.

You'll also need to install the base64decode package if you haven't already for decoding the base64 data. You can install it using install.packages("base64decode").









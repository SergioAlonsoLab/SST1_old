# Functions to read and write IMPPC sample barcodes


bcode.read <- function(barcode,what=c("all",
                                      "collection",
                                      "patient",
                                      "type",
                                      "sample")) {
  what <- match.arg(what)
  switch(what,
         all = barcode,
         collection = substr(barcode,1,3),
         patient = substr(barcode,5,9),
         type = substr(barcode,11,12),
         sample = substr(barcode,14,16))
  
}

bcode.read(imppc$Tumor.barcode,"s")

bcode.patient <- function(barcodes) {
  substr(barcodes,5,9)
  
}

bcode.patient(imppc$Tumor.barcode)

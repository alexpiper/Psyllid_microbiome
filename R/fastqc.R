fastqc_installer <- function (url, destdir = "bin") {
  if (missing(url)) {
    download_page <- xml2::read_html("http://www.bioinformatics.babraham.ac.uk/projects/download.html")
    link_hrefs <- download_page %>% rvest::html_nodes("a") %>% 
      rvest::html_attr("href")
    fastqc_href <- grep("fastqc/fastqc.*.zip", link_hrefs, 
                        perl = TRUE) %>% link_hrefs[.] %>% .[1]
    url <- paste0("http://www.bioinformatics.babraham.ac.uk/projects/", 
                  fastqc_href)
  }
  if(!dir.exists(destdir)){dir.create(destdir)} # Create first fastqc
  if(dir.exists(paste0(destdir,"/fastqc"))){unlink(paste0(destdir,"/fastqc"),recursive=TRUE)} # Remove old version
  
  destfile <- file.path(destdir, basename(url))
  if(exists(destfile)){file.remove(destfile)} # Remove old zip file
  
  utils::download.file(url, destfile = destfile)
  utils::untar(destfile,exdir = destdir)  ## check contents
  file.remove(destfile)
}

##Problem - this only seems to work when the files are in the same folder as fastqc
fastqc_local <- function (fq.dir, qc.dir = NULL, threads = 2, fastqc.path = "FastQC") {
  
  if (is.null(qc.dir))(qc.dir <- file.path(fq.dir, "FASTQC"))
    dir.create(qc.dir)
    
    current.wd <- getwd()
    setwd(fastqc.path)
    args <- paste0("-Xmx250m -cp .;.sam-1.103.jar;./jbzip2-0.9.jar uk.ac.babraham.FastQC.FastQCApplication ","../", fastqFs, collapse = " ")
    
    print(paste0(args," /n")) 
    
    system2("java", args = args)
    setwd(current.wd) 
}



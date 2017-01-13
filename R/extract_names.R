extract_nums <- function(lst){
    unlist(sapply(lst, 
                  function(x)
                      ifelse(length(x) > 2,
                             x[[3]]$name,
                             x[[2]])))
}

extract_all_nums<- function(lst){
    len <- max(which(substr(names(lst), 1, 1) == "T"))
    sapply(lst[1:len], function(x) extract_nums(x))
}

### naming in RVine routine
naming <-  function(numb){ 
    bef <- paste(as.character(numb[2]),
                 as.character(numb[1]),
                 sep=",",
                 collapse="")
    aft <- if(length(numb) > 2){
        gsub(" ",
             ",",
             do.call(paste, as.list(as.character(numb[3:length(numb)])))) 
    }  else ""
    sep <- if(length(numb) > 2) " ; " else ""
    paste(bef, aft, sep=sep, collapse="")
}

split_num <- function(x){
    unlist(lapply(strsplit(x, " ; "), strsplit, split = ","))
}
split_name <- function(x){
    unlist(lapply(strsplit(x, " ; ", fixed = TRUE), strsplit, split = ","))
}

# delete from directory

dir_from <- "/gscratch/home/jiang14/mhp/output"

setwd(dir_from)

filenames <- list.files()

for (i in 1:length(filenames)) {
  temp <- unlist(strsplit(filenames[i], split = '.', fixed = TRUE))
  l <- length(temp)
  if (temp[l] == "out") {
    file.remove(filenames[i])
  }


}

dir_from <- "/gscratch/csde/jiang14/output"

setwd(dir_from)

filenames <- list.files()

for (i in 1:length(filenames)) {
  temp <- unlist(strsplit(filenames[i], split = '.', fixed = TRUE))
  l <- length(temp)
  if (temp[l] == "out") {
    file.remove(filenames[i])
  }


}



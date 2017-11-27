## Function taken from VineCopula package

allfams <- c(0:10,
             13, 14, 16:20,
             23, 24, 26:30, 33, 34, 36:40,
             104, 114, 124, 134, 204, 214, 224, 234)
tawns <- which(allfams > 100)
onepar <- setdiff(which(allfams %% 10 %in% c(1, 3, 4, 5, 6)), tawns)
twopar <- seq_along(allfams)[-c(1, onepar)]

with_rotations <- function(nums) {
  unique(unlist(lapply(nums, get_rotations)))
}

get_rotations <- function(fam) {
  sgn <- sign(fam)  # indicator for negative selection
  fam <- sgn * fam  # ensure that fam is positive from here on

  if (fam %in% c(0, 1, 2, 5)) {
    # no roations for independence, gaussian, student and frank copulas
    out <- fam
  } else if (fam %in% c(3, 13, 23, 33)) {
    out <- c(3, 13, 23, 33)
  } else if(fam %in% c(4, 14, 24, 34)) {
    out <- c(4, 14, 24, 34)
  } else if(fam %in% c(6, 16, 26, 36)) {
    out <- c(6, 16, 26, 36)
  } else if(fam %in% c(7, 17, 27, 37)) {
    out <- c(7, 17, 27, 37)
  } else if(fam %in% c(8, 18, 28, 38)) {
    out <- c(8, 18, 28, 38)
  } else if(fam %in% c(9, 19, 29, 39)) {
    out <- c(9, 19, 29, 39)
  } else if(fam %in% c(10, 20, 30, 40)) {
    out <- c(10, 20, 30, 40)
  } else if(fam %in% c(104, 114, 124, 134)) {
    out <- c(104, 114, 124, 134)
  } else if(fam %in% c(204, 214, 224, 234)) {
    out <- c(204, 214, 224, 234)
  }

  # adjust for negative selection
  sgn * out
}

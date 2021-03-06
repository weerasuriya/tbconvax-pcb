normalize.contact.matrices <- function(C, popv, make.sym = F) {
  # function to normalize them so that
  # the matrix with all the contacts is normalized so that its dominant eigenvalue is 1
  # and other matrices keep their contributions
  if (make.sym) {
    Csym <- lapply(C, function(x, popv) (x + t(x) * ((popv) %*% t(1 / popv))) / 2, popv) # make sure contacts are reciprocal
  } else {
    Csym <- C # if make.sym = F leave it as is
  }
  eig1 <- Re(eigen(Csym["all"]$all)$value[1]) # save dominant eigenvalue of the matrix with all contacts
  # divide all matrices by the real part of the dominant matrix of all of the contacts
  # that way all the rescaled matrices still sum to C_all = C_work + C_home + C_school + C_other
  Cnorm <- lapply(Csym, function(x, eig1) x / eig1, eig1)
  return(Cnorm)
}

# Round number to any number
# see: https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x

rounder <- function(x,y) {
  if(y >= 0) { x + (y - x %% y)}
  else { x - (x %% abs(y))}
}
find.h1 <- function(x, eps=1e-4,h.max=100) {
  h.min=1e-6
  while ((h.max-h.min)>eps) {
    h.cur <- (h.min+h.max)/2
    cur.kde <- density(x, bw=h.cur)
    cur.d1 <- (cur.kde$y[2:512]-cur.kde$y[1:511]) /
        (cur.kde$x[2:512]-cur.kde$x[1:511])
    cur.d0 <- min(which(cur.d1<0))
    if (sum(cur.d1[cur.d0:511]>0)==0) {
      h.max <- h.cur
    }
    else {
      h.min <- h.cur
    }
  }
  return((h.min+h.max)/2)
}

rm(list = ls())
while(length(dev.list()) > 0){dev.off()}

# L2 inner product of f and g on a given interval
inner <- function(f, g, interval, mesh = 0.01){
  x <- seq(interval[1], interval[2], by = mesh)
  (interval[2] - interval[1]) * mean(f(x) * g(x))
}

# proj_f(g), projection of g onto f
fn.proj <- function(to.project.onto, to.be.projected, interval, mesh = 0.01){
  E <- new.env()
  E$a <- inner(f = to.project.onto, 
               g = to.be.projected, 
               interval, mesh) / inner(f = to.project.onto, 
                                       g = to.project.onto, 
                                       interval, mesh)
  E$f <- to.project.onto
  function(x){
    E$a * E$f(x)
  }
}

# orthogonalize one function against another
orthize <- function(to.orthogonalize.against, to.be.orthogonalized, interval, mesh = 0.01){
  E <- new.env()
  E$u <- to.orthogonalize.against
  E$v <- to.be.orthogonalized
  E$interval <- interval
  E$mesh <- mesh
  function(x){
    E$v(x) - fn.proj(to.project.onto = E$u,
                     to.be.projected = E$v,
                     interval = E$interval,
                     mesh = E$mesh)(x)
  }
}

# turn a basis of functions into an orthogonal basis
gram.schmidt <- function(fn.list, interval, mesh = 0.01){
  ortho.list <- vector("list", length = length(fn.list))
  for(k in 1:length(fn.list)){
    if(k == 1){
      ortho.list[[k]] <- fn.list[[k]] # first vector in orth. basis = first vector in basis
    }else{
      v <- fn.list[[k]] # fn from basis that will be orthogonalized against all the k-1 orthogonal functions weve already generated
      u <- ortho.list[1:(k - 1)] # the k-1 orth. functions we've already generated
      u.k <- vector("list", length = k) # to store the sequence of functions generated as we orth.ize v against each function in u, one-by-one
      u.k[[1]] <- v # the first function in this sequence is v, not yet orth.ized at all
      for(j in 2:k){
        u.k[[j]] <- orthize(to.orthogonalize.against = u[[j - 1]], # each previously generated orth. function appears here, in turn
                            to.be.orthogonalized = u.k[[j - 1]], # the previous element of the same list we're adding to
                            interval = interval,
                            mesh = mesh)
      }
      ortho.list[[k]] <- u.k[[k]] # the last element of u.k is what v looks like after orth.ization against all the other orth. functions we've generated 
    }
  }
  ortho.list
}

f <- function(x){
  x * exp(x - x^2)
}

g <- function(x){
  1 / (1 + x^2)
}

h <- function(x){
  sin(abs(x)^1.4) + sin(pi * x)
}

interval <- c(-3, 3)

fn.list <- list(f, g, h)

curve(f, from = interval[1], to = interval[2], lwd = 2, n = 1001, asp = 1)
curve(g, from = interval[1], to = interval[2], lwd = 2, n = 1001, asp = 1)
curve(h, from = interval[1], to = interval[2], lwd = 2, n = 1001, asp = 1)

f.g <- fn.proj(to.project.onto = f, to.be.projected = g, interval = interval)
f.h <- fn.proj(to.project.onto = f, to.be.projected = h, interval = interval)
g.h <- fn.proj(to.project.onto = g, to.be.projected = h, interval = interval)
curve(f.g, from = interval[1], to = interval[2], lwd = 2, n = 1001, asp = 1)
curve(f.h, from = interval[1], to = interval[2], lwd = 2, n = 1001, asp = 1)
curve(g.h, from = interval[1], to = interval[2], lwd = 2, n = 1001, asp = 1)


ortho.list <- gram.schmidt(fn.list = fn.list,
                           interval = interval)

ol1 <- ortho.list[[1]]
ol2 <- ortho.list[[2]]
ol3 <- ortho.list[[3]]
curve(ol1, from = interval[1], to = interval[2], lwd = 2, n = 1001, asp = 1)
curve(ol2, from = interval[1], to = interval[2], lwd = 2, n = 1001, asp = 1)
curve(ol3, from = interval[1], to = interval[2], lwd = 2, n = 1001, asp = 1)

# should all be near 0
inner(f = ol1, g = ol2, interval = interval)
inner(f = ol1, g = ol3, interval = interval)
inner(f = ol2, g = ol3, interval = interval)

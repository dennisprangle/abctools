nn.ent <-function (th, k = 4){

th <- as.matrix(th)
p <- ncol(th)
n <- nrow(th)

if (p == 1) {
	t <- .C(C_nnone, as.double(t(th)), as.integer(n),
        as.integer(k), D = as.double(0))
}
else {
	t <- .C(C_nnk, as.double(t(th)), as.integer(n),
	as.integer(p), as.integer(k), D = as.double(0))
}

ent <- t$D

ent
}


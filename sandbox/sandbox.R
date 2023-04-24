# ********************************************
# Cholesky tests
# April 2023
# ********************************************
rm(list=ls())

n <- 5; N <- 20
m <- matrix(rnorm(n * N), ncol=n)
R <- cor(m)
Ri <- solve(R)

# simple Cholesky, R = L %*% t(L) = L %*% U
U <- chol(R); L <- t(U)
all.equal(R, t(U) %*% U) # [1] TRUE

# Cholesky on Ri
U.Ri <- chol(Ri)
all.equal(Ri, t(U.Ri) %*% U.Ri) # [1] TRUE

# On L
Li <- solve(L)
all.equal(Ri, t(Li) %*% Li) # [1] TRUE

# Getting Ri from U from R, R^-1 = U^-1 %*% t(U^-1)
Ui <- solve(U) # is still upper triangular
all.equal(Ri, Ui %*% t(Ui)) # [1] TRUE
all.equal(Ri, t(Ui) %*% Ui) # Not true

Ui. <- backsolve(U, diag( dim( U )[1] ))
all.equal(Ui, Ui.) # [1] TRUE

# transform m: t(x) %*% R^-1 %*% x = t(x) %*% 
z1 <- m %*% Ui
m2 <- z1 %*% U
all.equal(m, m2) # [1] TRUE

# ********************************************
# mcmc2anc test
# April 2023
# ********************************************
rm(list=ls())

# prepare data:
data(carnivores)
M <- carnivores$C.proc
rownames(M) <- rownames(carnivores$M)
back.ages = as.numeric(unlist(strsplit(carnivores$tree$tip.label, "\\^"))[seq(from=2, to=2*19, by=2)])
back.ages = max(back.ages) - back.ages

# reconstruct:
recM1 <- mcmc2anc(carnivores$tree, M, carnivores$mcmc, "t_", "r_g1_", back.ages)
recM2 <- mcmc2anc(carnivores$tree, M, carnivores$mcmc, "t_", "r_g1_", back.ages, carnivores$var.foxes, carnivores$R.sh)
all.equal(recM1, recM2)
plot(recM1, recM2); abline(0, 1)
# This test confirmed that the reconstruction on the MCMC does not depend on the transform
# This test does not work anymore as the option to give R.sh and c to the code is removed

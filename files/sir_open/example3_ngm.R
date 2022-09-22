
# Construct matrix F containing expressions for all completely new infections
# entering infected compartments
F1 <- quote(lambda_1 * S_1)
F2 <- quote(lambda_2 * S_2)
F3 <- quote(lambda_3 * S_3)

# Construct a matrix V- containing expressions for all losses out of each 
# infected compartment
Vm1 <- quote(gamma * I_1)
Vm2 <- quote(gamma * I_2)
Vm3 <- quote(gamma * I_3)

# Contract a matri V+ containing expressions for all gains into each infected
# compartment that does not represent new infections, but transfers among
# infectious classes
Vp1 <- 0
Vp2 <- 0
Vp3 <- 0

# Substract V+ form V-
V1 <- substitute(a -b, list(a = Vm1, b = Vp1))
V2 <- substitute(a -b, list(a = Vm2, b = Vp2))
V3 <- substitute(a -b, list(a = Vm3, b = Vp3))

# Generate partial derivatives for the two Jacobian matrices f and v, which are
# derivatives for F and V with respect to the n infectious state variables
f11 <- D(F1, "I_1")
f21 <- D(F2, "I_2")
f31 <- D(F3, "I_3")

v11 <- D(V1, "I_1")
v21 <- D(V2, "I_2")
v31 <- D(V3, "I_3")

# Evaluate the matrices at the disease free equilibrium
lambda_1 <- quote((b * 7) + (b * 5) + (b * 1))

lambda_2 <- quote((b * 2) + (b * 9) + (b * 1))

lambda_3 <- quote((b * 1) + (b * 3) + (b * 2))

paras <- list(b = 0.02, gamma = 1 / 14, lambda_1 = lambda_1,
              lambda_2 = lambda_2, lambda_3 = lambda_3)

# Derive the greatest eigenvalue of fv-1|dfe
f_1 <- with(paras,
            matrix(c(eval(lambda_1), 0, 0,
                     0, eval(lambda_2), 0,
                     0, 0, eval(lambda_3)),
                   nrow = 3, byrow = TRUE))

v_1 <- with(paras,
            matrix(c(eval(v11), 0, 0,
                     0, eval(v21), 0,
                     0, 0, eval(v31)
            ),
            nrow = 3, byrow = TRUE))

max(eigen(f_1 %*% solve(v_1))$values)

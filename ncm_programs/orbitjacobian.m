syms u v p q real
r = (u^2+v^2)^(1/2)
F = [p; q; -u/r^3; -v/r^3]
y = [u ; v; p; q]
J = simple(jacobian(F,y))
r5J = simple(r^5*J)
lambda = simple(eig(J))
rplambda = simple(r^(3/2)*lambda)

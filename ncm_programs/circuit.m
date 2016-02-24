% Circuit parameters

r = round(100*rand(8,1))
v0 = round(100*randn)

% Kirchoff's voltage law

A = [1 -1  0  0
     0  1 -1  0
    -1  0  1  0
     0 -1  0  0
     0  0 -1  1
     1  0  0  0
     0  0  0 -1
    -1  0  0  1]

% Symbolically

symr = sym('[r12 r13 r14 r23 r34 r25 r35 r45]');
R = A'*diag(symr)*A

% Numerically 

R = A'*diag(r)*A;
b = [0 0 0 v0]'
i = R\[0 0 0 v0]'

% Kirchoff's current law

B = [1 -1  0  0
     1  0 -1  0
     1  0  0 -1
     0  1 -1  0
     0  0  1 -1
     0  1  0  0
     0  0  1  0
     0  0  0  1]

% Symbolically

symg = sym('[g12 g13 g14 g23 g34 g25 g35 g45]');
G = B'*diag(symg)*B

% Numerically 

g = 1./r;
g35 = g(7);
G = B'*diag(g)*B
c = [0 0 g35*v0 0]'
v = G\c

% Check consistency

d = [0 0 0 0 0 0 v0 0]';
[(B*v-d)./(A*i) r]

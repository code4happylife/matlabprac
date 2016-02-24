function x = pagerank2(G,p)
% PAGERANK2  Google's PageRank
% x = pagerank2(G,p) uses inverse iteration to compute the page
% rank for a connectivity matrix G with a damping factory p,
% (default is .85).
% See also: SURFER, PAGERANK1, PAGERANK3.

if nargin < 2, p = .85; end

% Eliminate any self-referential links

G = G - diag(diag(G));
  
% c = out-degree

[n,n] = size(G);
c = sum(G,1);

% Form the Markov transition matrix.

k = find(c~=0);
D = sparse(k,k,1./c(k),n,n);
e = ones(n,1);
z = ((1-p)*(c~=0) + (c==0))/n;
A = p*G*D + e*z;

% One step of inverse iteration to find Markov vector, A*x = x.

warning off MATLAB:nearlySingularMatrix
I = eye(n,n);
x = (I - A)\e;
x = x/sum(x);
warning on MATLAB:nearlySingularMatrix

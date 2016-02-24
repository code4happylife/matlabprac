function x = pagerank3(G,p)
% PAGERANK3  Google's PageRank
% x = pagerank1(G,p) uses the power method to compute the page
% rank for a connectivity matrix G with a damping factory p,
% (default is .85).
% See also: SURFER, PAGERANK1, PAGERANK2.

if nargin < 2, p = .85; end

% Eliminate any self-referential links

G = G - diag(diag(G));
  
% c = out-degree

[n,n] = size(G);
c = sum(G,1);

% Form the components of the Markov transition matrix.

k = find(c~=0);
D = sparse(k,k,1./c(k),n,n);
G = p*G*D;
e = ones(n,1);
z = ((1-p)*(c~=0) + (c==0))/n;

% Power method to find Markov vector, A*x = x.

x = e/n;
xs = zeros(n,1);
while norm(x-xs,inf)> 1.e-6*norm(x,inf)
   xs = x;
   x = G*x + e*(z*x);
end

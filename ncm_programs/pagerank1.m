function x = pagerank1(G,p)
% PAGERANK1  Google's PageRank
% x = pagerank1(G,p) uses direct solution of the sparse page
% rank equations to compute the page rank for a connectivity
% matrix G with a damping factory p, (default is .85).
%
% x = simple(pagerank1(G,sym('p')) returns a symbolic result.
%
% See also: SURFER, PAGERANK2, PAGERANK3.

if nargin < 2, p = .85; end

% Eliminate any self-referential links

G = G - diag(diag(G));
  
% c = out-degree

[n,n] = size(G);
c = sum(G,1);

% Scale column sums to be 1 (or 0 where there are no out links).

k = find(c~=0);
D = sparse(k,k,1./c(k),n,n);

% Solve (I - p*G*D)*x = e

e = ones(n,1);
I = speye(n,n);
x = (I - p*G*D)\e;

% Normalize so that sum(x) == 1.

x = x/sum(x);

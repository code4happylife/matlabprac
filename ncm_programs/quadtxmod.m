function [Q,fcount] = quadtx(F,a,b,tol,varargin)
%QUADTX  Evaluate definite integral numerically.
%   Q = QUADTX(F,A,B) approximates the integral of F(x) from A to B
%   to within a tolerance of 1.e-6.  F is a function handle or an
%   anonymous function that defines F(x).
%
%   Q = QUADTX(F,A,B,tol) uses the given tolerance instead of 1.e-6.
%
%   Arguments beyond the first four, Q = QUADTX(F,a,b,tol,p1,p2,...),
%   are passed on to the integrand, F(x,p1,p2,..).
%
%   [Q,fcount] = QUADTX(F,...) also counts the number of evaluations
%   of F(x).
%
%   See also QUAD, QUADL, DBLQUAD, QUADGUI.

% Default tolerance
if nargin < 4 | isempty(tol)
   tol = 1.e-6;
end

% Initialization
c = (a + b)/2;
fa = F(a,varargin{:});
fc = F(c,varargin{:});
fb = F(b,varargin{:});
fcount = 3;

% Fudge endpoints to avoid infinities.
if ~isfinite(fa)
   warning('Modifying endpoint')
   fa = F(a+eps*(b-a),varargin{:});
   fcount = fcount+1;
end
if ~isfinite(fb)
   warning('Modifying endpoint')
   fb = F(b-eps*(b-a),varargin{:});
   fcount = fcount+1;
end

% Recursive call 
[Q,fcount] = quadtxstep(F, a, b, tol, fa, fc, fb, fcount, varargin{:});
maxcnt = 10000;
if fcount >= maxcnt
   warning('Maximum function count exceeded.  Singularity likely.')
end

% ---------------------------------------------------------

function [Q,fcount] = quadtxstep(F,a,b,tol,fa,fc,fb,fcount,varargin)

% Recursive subfunction used by quadtx.

maxcnt = 10000;
h = b - a; 
c = (a + b)/2;
fd = F((a+c)/2,varargin{:});
fe = F((c+b)/2,varargin{:});
Q1 = h/6 * (fa + 4*fc + fb);
Q2 = h/12 * (fa + 4*fd + 2*fc + 4*fe + fb);
fcount = fcount + 2;
if (abs(Q2 - Q1) <= tol) | (fcount >= maxcnt)
   Q  = Q2 + (Q2 - Q1)/15;
else
   [Qa,fcount] = quadtxstep(F, a, c, tol, fa, fd, fc, fcount, varargin{:});
   [Qb,fcount] = quadtxstep(F, c, b, tol, fc, fe, fb, fcount, varargin{:});
   Q  = Qa + Qb;
end

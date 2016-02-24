function membraneshow(indices,m,n);
if nargin < 1, indices = 1:138; end
if nargin < 2, m = 30; end
if nargin < 3, n = 20; end
shg
clf
p = 0;
for k = indices(:)'
   p = p+1;
   subplot(2,2,p);
   [L,lam] = membranetx(k,m,n);
   L(m+1:2*m+1,m+1:2*m+1) = NaN;
   L = rot90(L);
   contourf(L,-1:.2:1);
   title(sprintf('%3d',k))
   text(m-8,-4,sprintf('%7.4f',lam))
   axis square
   axis off
   drawnow
   if p == 4 & k ~= indices(end)
      pause
      clf
      p = 0;
   end
end

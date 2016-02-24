function pennypois
% PENNYPOIS
% 2-D Poisson equation, source = penny, boundary values = 0.

load penny
P = flipud(P);

% Lighted surf plot of data

clf
shg
surf(P);
daspect([1,1,128])
colormap(copper)
shading interp
material metal
lighting gouraud
view(2)
axis tight
axis off
set(gca,'zlimmode','auto','climmode','manual');
light('pos',[1,2,2000],'style','inf');
toggle('next')

% Contour plot of data

contour(P,20)
axis square
colorbar
toggle('next')

% Contour plot of solution to Poisson equation

n = 128;
h = 1/(n+1);
e = ones(n,1);
D = spdiags([e -2*e e],[-1 0 1],n,n);
I = eye(n,n);
A = (1/h^2)*(kron(D,I) + kron(I,D));
u = reshape(A\P(:),n,n);
contour(u,20)
axis square
colormap(copper)
colorbar
toggle('next')

% Contour plot of Laplacian of solution

f = (4/h^2)*del2(u);
contour(f,20)
colorbar
axis square
toggle('close')
close(gcf)


% ------------------------

function toggle(str)
uic = uicontrol('units','norm','pos',[.02 .02 .10 .04], ...
   'string',str,'style','toggle');
while get(uic,'value') == 0
   drawnow
end

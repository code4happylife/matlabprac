%A4  
%   A4 paper.
sigma = sqrt(2);
x = [0 sigma sigma 0 0];
y = [0 0 1 1 0];
u = [sigma/2 sigma/2];
v = [0 1];
plot(x,y,'b',u,v,'b--')
text(sigma/2-.05,1.05,'sqrt(2)')
text(sigma/4-.05,-.05,'sqrt(2)/2');
text(3*sigma/4-.05,-.05,'sqrt(2)/2')
text(-.05,.50,'1')
title('A4 paper')
axis equal
axis off
set(gcf,'color','white')

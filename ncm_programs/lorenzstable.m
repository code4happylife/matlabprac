function lorenzstable
%LORENZSTABLE  Find rho so that critical point becomes stable.

   shg
   clf reset
   p = get(gcf,'pos');
   set(gcf,'color','black','doublebuff','on', ...
      'pos',[p(1) p(2)-(p(3)-p(4))/2 p(3) p(3)])
   stop = uicontrol('style','toggle','string','stop', ...
      'units','norm','pos',[.04 .02 .10 .04],'value',0);

   % Find rho so that two eigenvalues of the Jacobian lie on the
   % imaginary axis and the third lies on the negative real axis.

   rho = fzerotx(@lorenzeig,[20 28])

   sigma = 10;
   beta = 8/3;
   eta = sqrt(beta*(rho-1));
   A = [ -beta    0     eta
            0  -sigma   sigma 
         -eta   rho    -1  ];
   
   % The critical points are the null vectors of A.
   % The initial value of y(t) is near one of the critical points.
   
   yc = [rho-1; eta; eta];
   y0 = yc + [0; 0; 1];
   
   % Integrate forever, or until the stop button is toggled.
   
   tspan = [0 Inf];
   opts = odeset('reltol',1.e-6,'outputfcn',@lorenzplot,'refine',1);
   ode23(@lorenzeqn, tspan, y0, opts, A);



% ------------------------------

function ydot = lorenzeqn(t,y,A)
%LORENZEQN  Equation of the Lorenz chaotic attractor.
%   ydot = lorenzeqn(t,y,A).
%   The differential equation is written in almost linear form.
%      ydot = A*y
%   where
%      A = [ -beta    0     y(2)
%               0  -sigma   sigma 
%            -y(2)   rho    -1  ];

A(1,3) = y(2);
A(3,1) = -y(2);
ydot = A*y;


% ------------------------------

function fin = lorenzplot(t,y,job,A)
%LORENZPLOT   Plot the orbit of the Lorenz chaotic attractor.

if isequal(job,'init')

   % Initialize axis and comet, R = axis settings, L = length of comet.

   rho = A(3,2);
   R = [0 30  0 30 0 30];
   L = 150;
   set(gca,'color','black','pos',[.03 .05 .93 .95])
   axis(R);
   axis off

   comet = line(y(1),y(2),y(3),'linestyle','none','marker','.', ...
      'color','y','erasemode','none','markersize',6);

   beta = -A(1,1);
   eta = sqrt(beta*(rho-1));
   yc = [rho-1; eta; eta];
   line(yc(1),yc(2),yc(3),'linestyle','none','marker','o','color','g')
   line(yc(1),-yc(2),-yc(3),'linestyle','none','marker','o','color','g')

   ax = [R(2) R(1) R(1) R(1) R(1)];
   ay = [R(3) R(3) R(4) R(3) R(3)];
   az = [R(5) R(5) R(5) R(5) R(6)];
   p = .9;
   q = 1-p;
   grey = [.4 .4 .4];
   line(ax,ay,az,'color',grey);
   text(p*R(1)+q*R(2),R(3),p*R(5),sprintf('%3.0f',R(1)),'color',grey)
   text(q*R(1)+p*R(2),R(3),p*R(5),sprintf('%3.0f',R(2)),'color',grey)
   text(R(1),p*R(3)+q*R(4),p*R(5),sprintf('%3.0f',R(3)),'color',grey)
   text(R(1),q*R(3)+p*R(4),p*R(5),sprintf('%3.0f',R(4)),'color',grey)
   text(R(1),R(3),p*R(5)+q*R(6),sprintf('%3.0f',R(5)),'color',grey)
   text(R(1),R(3),q*R(5)+p*R(6),sprintf('%3.0f',R(6)),'color',grey)
   fin = 0;

   cameratoolbar('setmode','orbit')

elseif isequal(job,'done')

   stop = findobj('string','stop');
   set(stop,'string','close','callback','close(gcf)')
   fin = 1;

else

   % Update comet
   h = get(gca,'children');
   comet = h(end);
   set(comet,'xdata',y(1,end),'ydata',y(2,end),'zdata',y(3,end));
   drawnow;

   % Pause and restart

   stop = findobj('string','stop');
   fin = get(stop,'value');
end


% ------------------------------


function lambda = lorenzeig(rho)
beta = 8/3;
sigma = 10;
eta = sqrt(beta*(rho-1));
A = [ -beta    0     eta; 0 -sigma sigma; -eta rho -1];
y = [rho-1; eta; eta];
J = A + [0 y(3) 0; 0 0 0; 0 -y(1) 0];
lambda = eig(J);
lambda = real(lambda(1));

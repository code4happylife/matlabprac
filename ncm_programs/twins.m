function twins
% TWINS  Plot Tom's and Ben's weight.

%       Date       Tom     Ben
W = [10 27 2001    5 10    4  8  
     11 19 2001    7  4    5 11
     12 03 2001    8 12    6  4 
     12 20 2001   10 14    8  7
     01 09 2002   12 13   10  3
     01 23 2002   14  8   12  0
     03 06 2002   16 10   13 10];

t = datenum(W(:,[3 1 2]));
tom = W(:,4:5)*[1; 1/16];
ben = W(:,6:7)*[1; 1/16];
u = (min(t):max(t))';
h = plot(t,[tom ben],'o',u,[pchiptx(t,tom,u) pchiptx(t,ben,u)],'-');
set(h(3),'color',get(h(1),'color'))
set(h(4),'color',get(h(2),'color'))
axis([min(t)-14 max(t)+14 3 18])
datetick('x',12,'keeplimits')
legend({'Tom','Ben'},4)
title('Twins'' weights')
ylabel('Pounds')

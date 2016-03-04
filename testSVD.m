B = gallery(3)
% [A,B,C,...] = gallery(matname,P1,P2,...) returns the test matrices
% specified by the quoted string matname. The matname input is the name of
% a matrix family selected from the table below. P1,P2,... are input
% parameters required by the individual matrix family. The number of
% optional parameters P1,P2,... used in the calling syntax varies from
% matrix to matrix. The exact calling syntaxes are detailed in the
% individual matrix descriptions below.
% 
% [A,B,C,...] = gallery(matname,P1,P2,...,classname) produces a matrix of
% class classname. The classname input is a quoted string that must be
% either 'single' or 'double' (unless matname is 'integerdata', in which
% case 'int8', 'int16', 'int32', 'uint8', 'uint16', and 'uint32' are also
% allowed). If classname is not specified, then the class of the matrix is
% determined from those arguments among P1,P2,... that do not specify
% dimensions or select an option. If any of these arguments is of class
% single then the matrix is single; otherwise the matrix is double.
% 
% gallery(3) is a badly conditioned 3-by-3 matrix and gallery(5) is an
% interesting eigenvalue problem.
% 
% The gallery holds over fifty different test matrix functions useful for
% testing algorithms and other purposes.
svd(sym(B))
[U,S,V]=svd(B)

% testSVD
% 
% B =
% 
%   -149   -50  -154
%    537   180   546
%    -27    -9   -25
% 
%  
% ans =
%  
%                                                                                                                                                                                                          (447196886221/(9*(27^(1/2)*1875941161589056874332268^(1/2)*(i/27) + 299052995063664025/27)^(1/3)) + ((27^(1/2)*1875941161589056874332268^(1/2)*i)/27 + 299052995063664025/27)^(1/3) + 668737/3)^(1/2)
%  (668737/3 - ((27^(1/2)*1875941161589056874332268^(1/2)*i)/27 + 299052995063664025/27)^(1/3)/2 - (3^(1/2)*(447196886221/(9*(27^(1/2)*1875941161589056874332268^(1/2)*(i/27) + 299052995063664025/27)^(1/3)) - ((27^(1/2)*1875941161589056874332268^(1/2)*i)/27 + 299052995063664025/27)^(1/3))*i)/2 - 447196886221/(18*(27^(1/2)*1875941161589056874332268^(1/2)*(i/27) + 299052995063664025/27)^(1/3)))^(1/2)
%  (668737/3 - ((27^(1/2)*1875941161589056874332268^(1/2)*i)/27 + 299052995063664025/27)^(1/3)/2 + (3^(1/2)*(447196886221/(9*(27^(1/2)*1875941161589056874332268^(1/2)*(i/27) + 299052995063664025/27)^(1/3)) - ((27^(1/2)*1875941161589056874332268^(1/2)*i)/27 + 299052995063664025/27)^(1/3))*i)/2 - 447196886221/(18*(27^(1/2)*1875941161589056874332268^(1/2)*(i/27) + 299052995063664025/27)^(1/3)))^(1/2)
%  
% 
% U =
% 
%    -0.2691   -0.6798    0.6822
%     0.9620   -0.1557    0.2243
%    -0.0463    0.7167    0.6959
% 
% 
% S =
% 
%   817.7597         0         0
%          0    2.4750         0
%          0         0    0.0030
% 
% 
% V =
% 
%     0.6823   -0.6671    0.2990
%     0.2287   -0.1937   -0.9540
%     0.6944    0.7193    0.0204
C=gallery(3)
[X,lambda] = eig(C);
condest(X)
lambda=eig(C)
kappa=condeig(C)
format long
delta =1.e-6;
lambda=eig(C+delta*randn(3,3))
lambda - (1:3)'
delta*condeig(C)
% c = condest(A) computes a lower bound c for the 1-norm condition number
% of a square matrix A.
% 
% ans =
% 
%    1.2002e+03
D=gallery(5)
lambda=eig(D)
plot(eig(gallery(5)))
e=eig(D)
plot(real(e),imag(e),'r*',0,0,'ko')
axis(.1*[-1 1 -1 1])
axis square
e=eig(D+eps*randn(5,5) .*D)
% e=eig(D+eps*randn(5,5) .*D);plot(real(e),imag(e),'r*',0,0,'ko');axis(.1*[-1 1 -1 1]);axis square
% e=eig(D+eps*randn(5,5) .*D);plot(real(e),imag(e),'r*',0,0,'ko');axis(.1*[-1 1 -1 1]);axis square
% e=eig(D+eps*randn(5,5) .*D);plot(real(e),imag(e),'r*',0,0,'ko');axis(.1*[-1 1 -1 1]);axis square
E = gallery(5)
format long e
svd(E)
% ans =
% 
%      1.010353607103610e+05
%      1.679457384066870e+00
%      1.462838728085645e+00
%      1.080169069985621e+00
%      1.957628001581530e-14
while 1
    clc
    svd(E+eps*randn(5,5) .*E)
    pause(.25)
end
A=gallery(3)
[T,B]=schur(A)

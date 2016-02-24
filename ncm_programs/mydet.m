function d = mydet(A)
% DETERM  My version of det.
[L,U,P,sig] = lutxsig(A);
d = sig*prod(diag(U));

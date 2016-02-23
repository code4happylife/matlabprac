function y=mkron(varargin)
% y=mkron(x1,x2,x3,...) or y=mkron(x1,a)
% --------------------------------------
% Returns a Kronecker product of multiple input matrices x1,x2,x3... Needs
% at least 2 input arguments. Alternatively, one can provide one matrix x1
% and a scalar a to compute the repeated Kronecker product of x1 with
% itself.
%
% y         =   matrix, Kronecker product of x1,x2,x3..
%
% x1,2,3,.. =   matrix, matrix of arbitrary size.
%
% a         =   scalar, only used when given with 1 matrix input argument.


n=length(varargin);

if n<2
    error('Need at least 2 input arguments.')
end

% check whether second argument is a scalar, if so, then repeatedly apply
% kronecker product onto first argument
if isscalar(varargin{2})%当参数为[1,1]的标量时，isscalar返回逻辑1（true），否则返回逻辑0（false）
    n=varargin{2};%varargin is a MATLAB function;It is an input variable in a function definition statement
    y=varargin{1};%that allows the function to accept any number of input arguments.
    for i=2:n
        y=kron(y,varargin{1});%kron -Kronecker tensor product.
    end                       %This MATLAB function returns the Kronecker tensor product of X and Y.
else
    y=varargin{1};
    for i=2:n
        y=kron(y,varargin{i});
    end
end

end
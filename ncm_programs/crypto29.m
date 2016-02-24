function y = crypto29(x)

% Use a two-character Hill cipher with arithmetic modulo 29.

% Convert ASCII text to integers in range 0:96.
x(x=='.') = '}';
x(x==',') = '{';
x(x < 'a') = '|';
x(x > '}') = '|';
x = mod(real(x-'a'),29);

% Reshape into a matrix with 2 rows and floor(length(x)/2) columns.
n = 2*floor(length(x)/2);
X = reshape(x(1:n),2,n/2);

% Encode with matrix multiplication modulo 29.
A = [18 5; 5 11];
Y = mod(A*X,29);

% Reshape into a single row.
y = reshape(Y,1,n);

% If length(x) is odd, encode the last character. 
if length(x) > n
   y(n+1) = mod(28*x(n+1),29);
end

% Convert to ASCII characters.
y = char(y+'a');
y(y=='}') = '.';
y(y=='{') = ',';
y(y=='|') = ' ';

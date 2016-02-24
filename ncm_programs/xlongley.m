%     See: National Institute of Standards and Technology
%          Standard Reference Data Sets
%          Longley data set
%              http://www.itl.nist.gov/div898/strd/lls/data/Longley.shtml
%
%     Reference: Longley, J. W. (1967). 
%     An Appraisal of Least Squares Programs for the Electronic
%     Computer from the Viewpoint of the User. 
%     Journal of the American Statistical Association, 62, pp. 819-841.  
%
%     This classic dataset of labor statistics was one of the first used 
%     to test the accuracy of least squares computations. The response 
%     variable (y) is the Total Derived Employment and the predictor 
%     variables are GNP Implicit Price Deflator with Year 1954 = 100 (x1),
%     Gross National Product (x2), Unemployment (x3), Size of Armed Forces
%     (x4), Non-Institutional Population Age 14 & Over (x5), and Year (x6). 
%
%     Model:  
% 
%       y \approx beta_0 + sum_k=1^6 beta_k x_k
%

   load longley.dat
   [n,p] = size(longley);
   y = longley(:,1);
   X = longley(:,2:p);

%  Regression coefficients.
%  Include a column of ones to get beta_0.

   format long e
   e = ones(n,1);
   beta = [e X]\y

%  Plot y with residual error bars.

   figure(1)
   r = y - (beta(1) + X*beta(2:7));
   t = X(:,6);
   errorbar(t,y,abs(r))

%  Correlation coefficients.
%  Many of these values are close to 1, indicating that the columns
%  of X are highly correlated with each other
%  The correlation coefficients are computed by normalizing X
%  so that its columns have mean 0 and norm 1.

   format short
   corrcoeff = corrcoef(X)

%  Normalize y

   w = y - mean(y);
   w = w/norm(w);

%  Normalize X

   mu = mean(X);
   Z = X - mu(e,:);
   nrm = sqrt(diag(Z'*Z));
   Z = Z/diag(nrm);

%  Note that, after normalization, the response variable and the
%  six predictor variables are all highly correlated.

   figure(2)
   plot(t,Z,'-',t,w,'-ko')
   legend({'IPD','GNP','Unemployment','Armed Forces', ...
           'Population','Year','Total Derived Employment'},4)


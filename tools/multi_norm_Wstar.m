function [Wstar, pval] = multi_norm_Wstar(X)
% =================================================================
% This function calculates Villasenor Alva and Gonzalez Estrada's  W* test statistic for MVN.
%
% syntax : [Wstar, pval]= multi_norm_Wstar(X)
% inputs:
%   X: multivariate data, n x p matrix, n:sample size, p:dimensionality.
%
% outputs:
%   Wstar: test statistic
%   pval: p-value
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% References:
% Villasenor Alva J.,Gonz ?alez Estrada E. "A generalization of shapiro-wilk¡¦s test for multivariate normality."
% Communications in Statistics - Theory and Methods. 2009;38:1870¡V1883.
%
% =================================================================
[n,p] = size(X);
% Standardize input data X into a matrix with identity covariance matrix
if rank(cov(X)) < p % cov(X) is nonsingular                        
    Z=X;% ill-conditioned, no standardization
else
    Z=(chol(cov(X))'\X')';
end

Wstar=mean(W_stat_1992(Z));

y=log(n);
mun=-1.5861-0.31082*y-0.083751*y^2+0.0038915*y^3;
sn2=(exp(-0.4803-0.082676*y+0.0030302*y^2))^2;
s12=log((p-1+exp(sn2))/p);
mu1=mun+ 0.5*sn2 - 0.5*s12;
% cv=1 - exp(mu1 + sqrt(s12)*norminv(1-alfa))
pval= 1- normcdf((log(1-Wstar)-mu1)/sqrt(s12));

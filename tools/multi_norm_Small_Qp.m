function [Qp, P] = multi_norm_Small_Qp(X)
% =================================================================
% This function calculates Small's Omnibus statistics for Multivariate Normality test.
%
% syntax :[Q, P] = multi_norm_Small_Qp(X)
% Inputs :
%   X: multivariate data, n x p matrix, n:sample size, p:dimensionality.
%
% outputs :
%   Q: test statistics by Small, modified by Doornik and Hansen
%   P: p-value
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% References:
% Small NJH (1980). "Marginal skewness and kurtosis in testing multivariate normality." 
%           Applied Statistics, 29, 85-87. 
% D'Agostino RB (1970). "Transformation to Normality of the Null Distribution of g1." 
%           Biometrika, Vol 57, No.3 pp. 679-681.
% Doornik, JA, Hansen H. (2008). "An Omnibus Test for Univariate and Multivariate Normality." 
%     Oxford Bulletin of Economics and Statistics, 70, Supplement (2008), p.927-939.
%
% =================================================================

[n, p] = size(X);
X0=X - repmat(mean(X),n,1);
b1=mean(X0.^3)./mean(X0.^2).^(1.5);% skewness
% b1=skewness(X); % 1 x p, need statistic toolbox
b2=mean(X0.^4)./mean(X0.^2).^2;%kurtosis
% b2=kurtosis(X);% 1 x p, need statistic toolbox

U1 = corr(X).^3;
Z=b12Z(b1, n); % convert skewness to Z
q1=Z/U1*Z';  % statistic
        
U2 = corr(X).^4;
Z=b22Z_DH(b1,b2,n);% convert kurtosis to Z
% Z=b22Z(b2, n);% convert kurtosis to Z
q2=Z/U2*Z';% statistic
       
Qp=q1+q2;
P=1-chi2cdf(Qp,2*p);




function [Q, P] = multi_norm_Small(X)
% =================================================================
% This function calculates Small's Q1, Q2,Q3 statistics for Multivariate Normality test.
%
% syntax :[Q, P] = multi_norm_Small(X)
% Inputs :
%   X: multivariate data, n x p matrix, n:sample size, p:dimensionality.
%
% outputs :
%   Q: Q=[Q1 Q2 Q3], test statistics
%   P: P=[p1 p2 p3], p-values for Q1,Q2,Q3
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
% Ansombe FJ, Glynn WJ(1983). "Distribution of the Kurtosis Statistics b2 for Normal Samples." 
%           Biometrika, Vol 70, No.1 pp. 227-234.
%
% =================================================================

[n, p] = size(X);
X0=X - repmat(mean(X),n,1);
b1=mean(X0.^3)./mean(X0.^2).^(1.5);% skewness
% b1=skewness(X); % 1 x p, need statistic toolbox
b2=mean(X0.^4)./mean(X0.^2).^2;%kurtosis
% b2=kurtosis(X);% 1 x p, need statistic toolbox
% for Q1
        U1 = corr(X).^3;
        Z=b12Z(b1, n); % convert skewness to Z
        q1=Z/U1*Z';  % statistic
        p1=1-chi2cdf(q1,p);

% for Q2        
        U2 = corr(X).^4;
        Z=b22Z(b2, n);% convert kurtosis to Z
        q2=Z/U2*Z';% statistic
        p2=1-chi2cdf(q2,p);

% for Q3        
        q3=q1+q2;
        p3=1-chi2cdf(q3,2*p);

Q=[q1 q2 q3];
P=[p1 p2 p3];



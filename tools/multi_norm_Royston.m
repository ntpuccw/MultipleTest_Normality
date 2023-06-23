function [H,P_Royston]=multi_norm_Royston(X)
% This function computes Royston's Multivariate Normal test statistic.
%
% syntax : [H,P_Royston]=multi_norm_Royston(X)
%
% inputs :
%   X: multivariate data, n x p matrix, n:sample size, p:dimensionality.
% 
% outputs :
%   H: statistic
%   P_Royston: p-value.
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% References:
% Royston JP (1983). "Some techniques for assessing multivariate normality based on the Shapiro-Wilk $W$." 
%           Appl. Statist., 32, 121-133.	
% Royston JP (1992). "Approximating the Shapiro-Wilk W test for non-normality." 
%       Stat. Comput. 2,117-119.
% =================================================================

    [n,p]=size(X);
    c_lambda=5;
    c_mu=0.715;
    c_nu=polyval([-0.0018034 0.015124 0 0.21364], log(n));

    U=triu(corr(X),1);%compute correlation matrix
    rho_hat=U(U>0);%p x 1, retrieve correlation coefficients r_ij for i ~= j

    cij_bar=rho_hat.^c_lambda.*(1-c_mu/c_nu*(1-rho_hat).^c_mu);%p x 1
    c_bar=sum(cij_bar)/(p^2-p);%1 x 1
    e=p./(1+(p-1)*c_bar);
    W=W_stat_1992(X);%1 x p
    P=wtest1992(W,n);%1 x p
    K=norminv(P/2).^2;
    G=mean(K,2);
    H=G.*e;  
    P_Royston = 1-chi2cdf(H,e); % p-value
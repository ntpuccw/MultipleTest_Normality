function X=multi_Khintchine(p,n)
%=================================================================
% This function generates multivariate Khintchine random numbers.
%
% syntax : X=multi_Khintchine(p,n)
% inputs :
%   p:  p-variate
%   n:  sample size
%
% outputs :
%   X: n x p multivariate Khintchine matrix
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% Reference: "Multivariate Statistical Simulation," by Mark E. Johnson
%=================================================================

    Z=sqrt(chi2rnd(3,n,1));
    Z=repmat(Z,1,p);
    U=unifrnd(0,1,n,p);
    X=Z.*(2*U-1);    

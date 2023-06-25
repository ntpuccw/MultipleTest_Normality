function X=multi_PearsonII(p,m,n)
% =================================================================
% This function generate random numbers from multivariate Pearson Type II.
% 
% syntax : X=multi_Pearson(p,m,n)
%
% inputs :
%   p:  p-variates
%   m:  Pearson Type II parameter
%   n:  sample size
%
% outputs :
%   X: multivariate data
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% Reference: "Multivariate Statistical Simulation," by Mark E. Johnson
% =================================================================

    mu=zeros(1,p);
%     Sigma=(2*m+p+2)*eye(p);
%     B=chol(Sigma,'lower');% Sigma should be positive definite, used for
%     general Sigma.
    B=sqrt(2*m+p+2)*eye(p);% for diagonal Sigma.
    R=sqrt(betarnd(p/2,m+1,n,1));% n x 1
    %--- generate random numbers on the unit sphere
%     mu_for_U=zeros(1,p);
%     Sigma_for_U=eye(p);
%     Y=mvnrnd(mu_for_U,Sigma_for_U,n);% n x p
    Y = randn(n,p);% n x p, independent standard normal data
    U=Y./repmat(sqrt(sum(Y.^2,2)),1,p); % n x p
    %------- generate by the Cambanis representation -------------
    X=repmat(R,1,p).*(U*B')+repmat(mu,n,1);% n x p, p.110, 116
return
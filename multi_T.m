function X=multi_T(p,nu,rho,n)
% This function generates random numbers from General Multivariate T
% distribution with nu degrees of freedom. 
% The Pearson Type VII with parameter m is one of its special case where
% m=(p+nu)/2 and nu*Sigma_pearson=Sigma_t
% The Multivariate Cauchy distribution is a special case of it with nu=1.
%
% syntax : X=multi_T(p,nu,rho,n)
% 
% inputs :
%   p: p-variates
%   nu: degree of freedom
%   rho: the correlation coefficient of Sigma
%   n:sample size
%
% outputs :
%   X:multivariate data
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% Reference: "Multivariate Statistical Simulation," by Mark E. Johnson

    mu=zeros(1,p);
    Sigma=eye(p);Sigma(Sigma==0)=rho;%equal rho       
    R=chol(Sigma);
    Z=repmat(mu,n,1) + randn(n,p)*R;
%     Z=mvnrnd(zeros(1,p), Sigma, n);% n x p
    S=sqrt(nu./repmat(chi2rnd(nu,n,1),1,p)); % n x p
    X=S.*Z + repmat(mu,n,1);

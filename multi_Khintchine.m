function X=multi_Khintchine(p,n, marginal)
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
    if marginal == "normal"
% Khintchine distribution with normal marginals
        Z=sqrt(chi2rnd(3,n,1)); %page 156 Johnson, independent Uniform
%         Z =sqrt(3*gamma(1.5)/gamma(2.5)*gamrnd(1.5,1,n,1)) % the same,
%         p.36 Johnson
        Z=repmat(Z,1,p); % Identical generator
        U=unifrnd(0,1,n,p);
        X=Z.*(2*U-1);    
        return
    end
    if marginal == "exponential"
        % Khintchine distribution with General Exponential Power marginals 
        a = 1.5;
        all_p = [2 3 4 5 7 10];
        % the values of tau for each p is calculated in
        % determine_KHN_tau.mlx
        tau_for_p = [0.399 0.34 0.3 0.271 0.231 0.194]; % for p=[2 3 4 5 7 10]
        tau = tau_for_p(find(p==all_p,1));
        Z=sqrt(3*gamma(a)/gamma(a+2*tau))*gamrnd(a,1,n,1).^tau; %p.36 Johnson
        Z=repmat(Z,1,p); % Identical generator
        U=unifrnd(0,1,n,p);
        X=Z.*(2*U-1); 
        
% Khintchine distribution with Exponential marginals     
%         Z=sqrt(chi2rnd(3,n,p));
%         U=unifrnd(0,1,n,1);
%         U=repmat(U,1,p);
%         X=Z.*U;
        return
    end


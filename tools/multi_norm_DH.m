function [Ep, pval]= multi_norm_DH(X)
% =================================================================
% This function calculates DH's multivariate omnibus statistic Ep.
%
% syntax : [Ep, pval]= multi_norm_DH(X)
% inputs:
%   X: multivariate data, n x p matrix, n:sample size, p:dimensionality.
%
% outputs:
%   Ep: test statistic
%   pval: p-value
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% References:
% Doornik, JA, Hansen H. (2008). "An Omnibus Test for Univariate and Multivariate Normality." 
%     Oxford Bulletin of Economics and Statistics, 70, Supplement (2008), p.927-939.
%
% =================================================================

[n,p] = size(X);
% C=corr(X);% need statistics toolbox
% also by
% S = cov(X,1);
% V=diag(1./sqrt(diag(S)));%V^(-1/2)
% C=V*S*V;
C=corrcoef(X);
rp=rank(C);
% need to cope with rank deficiency suggested by Doornik, JA, Hansen H. (2008)
if rp < p
    [H,Lambda] = eig(C);
    [tmp,I]=sort(diag(Lambda));
    X=X*H(:,I(p-rp+1:end));
    C=corr(X);
    p=rp;
end
% end of rank deficiency

% standardization such that Y*Y'/n=Ip
S = cov(X,1);
[H,Lambda] = eig(C);
L=diag(1./sqrt(diag(Lambda)));% Lambda^(-1/2)
V=diag(1./sqrt(diag(S)));%V^(-1/2)
Y=H*L*H'*V*((X - repmat(mean(X),n,1)))';% p x n
% end of standardization

Y0=Y' - repmat(mean(Y'),n,1);
B1=mean(Y0.^3)./mean(Y0.^2).^(1.5);% skewness
% B1=skewness(Y');% 1 x p,  need statistics toolbox
B2=mean(Y0.^4)./mean(Y0.^2).^2;%kurtosis
% B2=kurtosis(Y');% 1 x p, need statistics toolbox
Z1=b12Z(B1,n);
Z2=b22Z_DH(B1,B2,n);
Ep=Z1*Z1'+Z2*Z2';
pval = 1 - chi2cdf(Ep, 2*p);

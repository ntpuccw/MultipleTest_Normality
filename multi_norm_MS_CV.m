function b1p = multi_norm_MS_CV(X, c)
% Mardia Skewness for MVN
n = size(X, 1);
dX = (X - repmat(mean(X),n,1));% (x_i - x_bar) for all i
if c == 1  % normalized by n
    S = cov(X,1);
else   % normalized by n-1
    S = cov(X);
end
% G = dX*inv(S)*dX';  
G = dX*(S\dX');  
b1p = sum(sum(G.^3))/n^2;  % 
% K = ((p+1)*(n+1)*(n+3))/(n*(((n+1)*(p+1))-6));  % correction for n>=50
% v = p*(p+1)*(p+2)/6; 
% % if n>50
% %     g = (n*b1p)/6;  % MS test statistics ~ chi2(v)
% % else
%     g = (n*b1p*K)/6;  % for small sample
% % end
% pval = 1 - chi2cdf(g, v);  % p-value of skewness to the right

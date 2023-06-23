function b2p = multi_norm_MK_CV(X, c)
% Mardia Kurtosis for MVN
% This program is specially designed for generating statistics only for
% collecting Critical Values.

% [n, p] = size(X);
n = size(X, 1);
dX = (X - repmat(mean(X),n,1));% (x_i - x_bar) for all i
if c == 1  % normalized by n
    S = cov(X,1);
else   % normalized by n-1
    S = cov(X);
end
% G = dX*inv(S)*dX';  
G = dX*(S\dX');  
b2p=trace(G.^2)/n;  % 
% g2 = (b2p-((p*(p+2)*(n-1))/(n+1)))/(sqrt((8*p*(p+2))/n));  % ~ N(0,1)
% pval = 2*(1-normcdf(abs(g2)));  % p-value of kurtosis Âù§À
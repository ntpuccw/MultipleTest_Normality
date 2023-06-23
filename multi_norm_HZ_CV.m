function [HZ_h141, HZ_hS, HZ_hL] = multi_norm_HZ_CV(X,c)
% Calculate statistic of Henze-Zirkler's Multivariate Normality Test.

if c == 1  %covariance matrix normalizes by (n) [=default]
    S = cov(X,1);
else   %covariance matrix normalizes by (n-1)
    S = cov(X);
end

[n,p] = size(X);

difT = (X - repmat(mean(X),n,1)); 

% Dj = diag(difT*inv(S)*difT');  %squared-Mahalanobis' distances
Dj = diag(difT*(S\difT'));
% Y = X*inv(S)*X';
Y = X*(S\X');

Djk = - 2*Y' + diag(Y')*ones(1,n) + ones(n,1)*diag(Y')'; %this procedure was
                       %taken into account in order to '..avoiding loops and
                       %it (the file) runs much faster' we thank to Johan
                       %(J.D.) for it valuabe comment (15/12/08) to improve 
                       %it.

% b = 1/(sqrt(2))*((2*p + 1)/4)^(1/(p + 4))*(n^(1/(p + 4))); % (2.5) smoothing parameter
h = 1.41;
HZ_h141 = HZbyh(Djk,Dj,S,n,p,h);
h = 0.448 + 0.026*p;
HZ_hS = HZbyh(Djk,Dj,S,n,p,h);
h = 0.928 + 0.049*p;
HZ_hL = HZbyh(Djk,Dj,S,n,p,h);

function HZ = HZbyh(Djk,Dj,S,n,p,h)
    b=1/(h*sqrt(2));

    if (rank(S) == p)    
        HZ = n * (1/(n^2) * sum(sum(exp( - (b^2)/2 * Djk))) - 2 *...
            ((1 + (b^2))^( - p/2)) * (1/n) * (sum(exp( - ((b^2)/(2 *...
            (1 + (b^2)))) * Dj))) + ((1 + (2 * (b^2)))^( - p/2))); % equation (2.1) (2.6)
    else
        HZ = n*4;
    end

return
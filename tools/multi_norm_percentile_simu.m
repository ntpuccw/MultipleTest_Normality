function wmin_m = multi_norm_percentile_simu(X,m,q)
% =================================================================
% This function generates the empirical q*100 th percentile of the Wmin_m(q)
% statistic from a set of W-tests, which contains m statistics W(X*c_i),
% i=1,2,...,m. c_i are m uniformly scattered points in p-dimensional uint
% sphere.
% limit: the dimension of the input data is restricted to between 2 and 10
% because of the availability of the empirical critical values. The sample
% size is suggested to be bewteen [10 200]. 
%
% syntax : wmin_m=multi_norm_percentile(X,p,m,q)
%
% inputs :
%   X: multivariate data with dimension in [2 10].
%   m: how many uniformly scattered points, 10000 is suggested to be the
%   minimum.
%   q: the percentile e.g. q=0.05 for the 5th percentile. q can be a
%   vector.
%
% outputs :
%   wmin_m : the wmin_m(q) test statistic from m W-tests. 
%   pavl: p-value(s) computed by interpolation from the attached empirical critical values.
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% References:
% Wang CC, Hwang YT (2010). "A new functional statistic for multivariate normality." 
%       Statistics and Computing, Statistics and Computing, 21(4), 501-509.
% =================================================================
%

[n,p]=size(X);
% if ( p < 2 || p > 10)
%     error('Invalid data dimension.');
% end
% Standardize input data X into a matrix with identity covariance matrix
if rank(cov(X)) < p % cov(X) is nonsingular                        
    Z=X;% ill-conditioned, no standardization
else
    Z=(chol(cov(X))'\X')';
end

%-- The Box-Muller method
%     mu=zeros(1,p);
%     Sigma=eye(p);
%     Y=mvnrnd(mu,Sigma,m);% m x p
    Y = randn(m,p);% Independent multivariate normal.
%     YL=sort((Y./repmat(sqrt(sum(Y.^2,2)),1,p))*Z',2); % linear combinations
    YL=(Y./repmat(sqrt(sum(Y.^2,2)),1,p))*Z'; % linear combinations, m x n
    % for the sake of insufficient memory,use loop.
    if n*m < 1e07
        R=W_stat_1992(YL');% (1 x m) W test stats   
    else
        R=zeros(1,m);
        for i=1:m
            R(i)=W_stat_1992(YL(i,:));
        end
    end
    Rr=sort(R);
    wmin_m=Rr(ceil(m*q));
    
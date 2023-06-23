function Z=b22Z_DH(b1,b2,N)
% =================================================================
% This function transforms the kurtosis b2 from a gamma distribution to chi2, then to standard Normal Z. 
%
% syntax : Z=b22Z(b1,b2,N)
% inputs :
%       b1:  the skewness derived by the conventional defintion for sqrt(b1).
%       b2:  the kurtosis derived by the conventional defintion for b2. It can be a value or a vector.
%       N:   sample size from which b2 derived. N>=8 is required to have
%            better approximation.
% outputs :
%       Z: sample(s) from N(0,1) and have the same size as the input b2.
%
% Written by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% Reference:
% Doornik, JA, Hansen H. (2008). "An Omnibus Test for Univariate and Multivariate Normality." 
%     Oxford Bulletin of Economics and Statistics, 70, Supplement (2008), p.927-939.
% =================================================================
delta=(N-3)*(N+1)*(N^2+15*N-4);
a=(N-2)*(N+5)*(N+7)*(N^2+27*N-70)/(6*delta);
c=(N-7)*(N+5)*(N+7)*(N^2+2*N-5)/(6*delta);
k=(N+5)*(N+7)*(N^3+37*N^2+11*N-313)/(12*delta);
alfa=a+b1.^2*c;
X=(b2-1-b1.^2)*2*k;
Z=(nthroot(X./(2*alfa),3)-1+1./(9*alfa)).*sqrt(9*alfa);

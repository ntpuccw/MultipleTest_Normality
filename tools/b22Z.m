function Z=b22Z(b2,N)
% =================================================================
% This function transforms the null distribution of b2 to standard Normal Z. 
%
% syntax : Z=b22Z(b2,N)
% inputs :
%       b2:  the kurtosis derived by the conventional defintion for b2. It can be a value or a vector.
%       N:   sample size from which b2 derived. N>=20 is required to have
%            better approximation.
% outputs :
%       Z: sample(s) from N(0,1) and have the same size as the input b2.
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% Reference:
% Ansombe FJ, Glynn WJ (1983). "Distribution of the Kurtosis Statistics b2 for Normal Samples." 
%                   Biometrika, Vol 70, No.1 pp. 227-234.
% =================================================================

Eb2=3*(N-1)/(N+1);
Vb2=24*N*(N-2)*(N-3)/((N+1)^2*(N+3)*(N+5));
SBb2=6*(N^2-5*N+2)/((N+7)*(N+9))*sqrt(6*(N+3)*(N+5)/(N*(N-2)*(N-3)));
A=6+8/SBb2*(2/SBb2 + sqrt((1+4/SBb2^2)));
x=(b2-Eb2)/sqrt(Vb2);
Z=((1-2/(9*A)) - nthroot((1-2/A)./(1 + x*sqrt(2/(A-4))),3))/sqrt(2/(9*A));
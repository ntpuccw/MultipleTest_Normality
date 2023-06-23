function Z=b12Z(g1,N)
% =================================================================
% This function transforms the null distribution of sqrt(b1) to standard Normal Z. 
%
% syntax : Z=b12Z(g1,N)
% inputs :
%       g1: the skewness derived by the conventional defintion for
%           sqrt(b1) . g1 can be a value or a vector.
%       N:  sample size from which g1 derived. N>=8 is required to have
%           better approximation.
%
% outputs :
%       Z: sample(s) from N(0,1) and have the same size as the input g1.
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% Reference:
% D'Agostino RB (1970). "Transformation to Normality of the Null Distribution of g1." 
%       Biometrika, Vol 57, No.3 pp. 679-681
%
% =================================================================

Y=g1*sqrt((N+1)*(N+3)/(6*(N-2)));
beta2=3*(N^2+27*N-70)*(N+1)*(N+3)/((N-2)*(N+5)*(N+7)*(N+9));
W2=sqrt(2*beta2-2)-1;
delta=1./sqrt(log(sqrt(W2)));
alfa=sqrt(2/(W2-1));
Z=delta*asinh(Y/alfa);

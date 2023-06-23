function p=wtest1992(W,n)
% This function transforms the W-statistic into Z (output p-value) based on Royston(1992) 
% Note: The less W is, the less normality. Small W causes large Z.
%
% syntax : p=wtest1992(W,n)
%
% inputs :
%   W : W test statistic
%   n : sample size
%
% outputs:
%
%   p : p-value
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% References:
% Royston JP (1992). "Approximating the Shapiro-Wilk W test for non-normality." 
%       Stat. Comput. 2,117-119.
% =================================================================

if n>=12
    w=log(1-W);
    mu=-1.5861-0.31082*log(n)-0.083751*log(n)^2+0.0038915*log(n)^3;
    sigma=exp(-0.4803-0.082676*log(n)+0.0030302*log(n)^2);
    Z=(w-mu)/sigma;
    p=1-normcdf(Z);
elseif n>=4 && n<=11
    g=-2.273+0.459*n;
	w=-log(g-log(1-W));
    mu=0.544-0.39978*n+0.025054*n^2-0.0006714*n^3;
    sigma=exp(1.3822-0.77857*n+0.062767*n^2-0.0020322*n^3);
    Z=(w-mu)/sigma;
    p=1-normcdf(Z);
end
% special concern for n=3
if n==3
    p=6/pi*(asin(sqrt(W))-asin(sqrt(0.75)));% p=F(W), the cumulative probability of W statistic
end
function W=W_stat_1992(y)
%
% This function computes W statistic for data with sample size from  3 to 1000 based on
% Royston (1992).
%
% syntax : W=W_stat_1992(y)
% 
% inputs :
%   y: data, can be a vector or a matrix (column-wide)
% 
% outputs :
%   W: W test statistic
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% References:	
% Royston JP (1992). "Approximating the Shapiro-Wilk W test for non-normality." 
%       Stat. Comput. 2,117-119.
% =================================================================

n=size(y,1);
if n==1
    y=y';%convert to column vector
    n=length(y);
end
y=sort(y);
 i=1:n;
 m=norminv((i-3/8)/(n+1/4)); % approximate mean
 c=1/sqrt(m*m')*m;
 x=1/sqrt(n);
 %--- compute coefficients a1,a2,...,an
 a=zeros(1,n);
 a(end)=c(end)+0.221157*x-0.147981*x^2-2.071190*x^3+4.434685*x^4-2.706056*x^5;
 a(end-1)=c(end-1)+0.042981*x-0.293762*x^2-1.752461*x^3+5.682633*x^4-3.582663*x^5;
 if n<=5
     phi=(m*m'-2*m(end)^2)/(1-2*a(end)^2);
     a(2:end-1)=m(2:end-1)/sqrt(phi);
     a(1)=-a(end);
 else
     phi=(m*m'-2*m(end)^2-2*m(end-1)^2)/(1-2*a(end)^2-2*a(end-1)^2);
     a(3:end-2)=m(3:end-2)/sqrt(phi);
     a(1:2)=-[a(end) a(end-1)];
 end
 %-- In Royston' paper, n=3 is excluded.
 if n==3%set theoretical value for a at n=3
     a=[-0.7071 0 0.7071];
 end
 %-- compute W-statistic
 W=(a*y).^2./((n-1)*var(y));% W-statistic

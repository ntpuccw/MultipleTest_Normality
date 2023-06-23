clear
alfa=0.2;
n=100;
p=2;

y=log(n);
mun=-1.5861-0.31082*y-0.083751*y^2+0.0038915*y^3;
sn2=(exp(-0.4803-0.082676*y+0.0030302*y^2))^2;
s12=log((p-1+exp(sn2))/p);
mu1=mun+ 0.5*sn2 - 0.5*s12;
cv=1 - exp(mu1 + sqrt(s12)*norminv(1-alfa))
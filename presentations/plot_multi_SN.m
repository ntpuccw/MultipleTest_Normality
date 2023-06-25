% This program draws a contour plot for the skewed normal fuction.
clear
x1 = -3:.2:3; x2 = -3:.2:3;
[X1,X2] = meshgrid(x1,x2);
a=[2 2];% shape factor
rho=0.5;
mu=zeros(1,2);
Sigma=eye(2);Sigma(Sigma==0)=rho;
F = dmsn([X1(:) X2(:)], mu, Sigma, a);%compute in a vector format 
F = reshape(F,length(x2),length(x1));%reshape into a matrix
contour(X1,X2,F,30)
grid
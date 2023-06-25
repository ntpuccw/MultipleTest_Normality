% This program make a function call to generate multivariate T
% random numbers. Then the locations of these random numbers are plotted
% along with the corresponding distribution function.
clear
% --- Draw graph for p=2 for confirmation
p=2;% dimension
nu=3;% Multivariate T parameter
n=1000;% sample size
rho=0;
Sigma=eye(p);Sigma(Sigma==0)=rho;
X=multi_T(p,nu,rho,n);
m=(nu+p)/2;
D=diag(X*(Sigma\X'));% n x n
% z=gamma(m)/gamma(nu/2)/pi^(p/2)/sqrt(det(Sigma))*((1+D).^(-m)); %n x 1
z=gamma(m)/gamma(nu/2)/(pi*nu)^(p/2)/sqrt(det(Sigma))*((1+D/nu).^(-m)); %n x 1
xlin = linspace(min(X(:,1)),max(X(:,1)),33);
ylin = linspace(min(X(:,2)),max(X(:,2)),33);
% xlin = linspace(-2,2,33);
% ylin = linspace(-2,2,33);
[XX,YY] = meshgrid(xlin,ylin);
Z = griddata(X(:,1),X(:,2),z,XX,YY,'cubic');
meshc(XX,YY,Z) %interpolated
axis tight; hold on
plot3(X(:,1),X(:,2),z,'.','MarkerSize',15) 
hold off
title(strcat('Multivariate T Distribution with \nu=' ,num2str(nu)))
% This program make a function call to generate multivariate Pearson Type
% II random numbers. Then the locations of these random numbers are plotted
% along with the corresponding distribution function.

% --- Draw graph for p=2 
p=2;% dimension
m=0.5;% Pearson Type parameter
n=1000;% sample size
mu=zeros(1,p);
Sigma=(2*m+4)*eye(p);
X=multi_PearsonII(p,m,n);
D=diag(X*(Sigma\X'));% n x n
z=gamma(p/2+m+1)/gamma(m+1)/pi^(p/2)/sqrt(det(Sigma))*((1-D).^m); %n x 1
% xlin = linspace(min(X(:,1)),max(X(:,1)),33);
% ylin = linspace(min(X(:,2)),max(X(:,2)),33);
xlin = linspace(-3,3,33);
ylin = linspace(-3,3,33);
[XX,YY] = meshgrid(xlin,ylin);
Z = griddata(X(:,1),X(:,2),z,XX,YY,'cubic');
meshc(XX,YY,Z) %interpolated
axis tight; hold on
plot3(X(:,1),X(:,2),z,'.','MarkerSize',15) 
hold off
% title(strcat('Pearson Type II(m=',num2str(m),')'))
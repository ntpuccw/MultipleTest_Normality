% This program generates mixed MVN random numbers. 
% The locations of these random numbers are plotted
% along with the corresponding distribution function.

% --- Draw graph for p=2 
p=2;%dimension
n=1000;% sample size
pi1=0.1;pi2=1-pi1;% proportion
mu1=zeros(1,p);
Sigma1=eye(p);%equal rho            
mu2=0*ones(1,p);
Sigma2=Sigma1;
Sigma2(Sigma2==0)=0.9;Sigma2(Sigma2==1)=1;
N1=binornd(n,pi1);
N2=n-N1;
X1=mvnrnd(mu1, Sigma1, N1);
X2=mvnrnd(mu2, Sigma2, N2);
X=[X1;X2];% random data
%-- make a plot
z=pi1*mvnpdf(X,mu1,Sigma1)+pi2*mvnpdf(X,mu2,Sigma2); %n x 1
% xlin = linspace(min(X(:,1)),max(X(:,1)),33);
% ylin = linspace(min(X(:,2)),max(X(:,2)),33);
xlin = linspace(-3,5,33);
ylin = linspace(-3,5,33);
[XX,YY] = meshgrid(xlin,ylin);
Z = griddata(X(:,1),X(:,2),z,XX,YY,'cubic');
meshc(XX,YY,Z) %interpolated
axis tight; hold on
plot3(X(:,1),X(:,2),z,'.','MarkerSize',15) 
hold off
% title('mixed MVN')
zlim([-0.1, 0.35])
% This program make a function call to generate Khintchine GEP distribution
% random numbers. The locations of these random numbers are plotted
% along with the corresponding distribution function (use MVN instead)

% --- Draw graph for p=2 for confirmation

p=2;
X=multi_Khintchine(p,2000, "exponential");
x1=X(:,1); x2=X(:,2);
mu = mean(X);
Sigma = cov(X);
z=mvnpdf(X,mu,Sigma); % use MVN instead; can't derive pdf of Khintchine GEP

xlin = linspace(-5,5,100);
ylin = linspace(-5,5,100);
[XX,YY] = meshgrid(xlin,ylin);
Z = griddata(x1,x2,z,XX,YY,'natural');
meshc(XX,YY,Z) %interpolated
axis tight; hold on
plot3(x1,x2,z,'.','MarkerSize',15) 
zlim([-0.05, 0.15])
hold off



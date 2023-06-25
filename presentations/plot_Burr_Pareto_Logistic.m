% This program make a function call to generate Burr-Pareto-Logistic distribution
% random numbers. The locations of these random numbers are plotted
% along with the corresponding distribution function.

% --- Draw graph for p=2 for confirmation
a = 1; % alpha = 1
[X,Y]=meshgrid(-3:0.1:3);
Z=(a+1)/a*normpdf(X).*normpdf(Y)./(normcdf(X).*normcdf(Y)).^(1/a+1)./ ...
    (1./normcdf(X).^(1/a) + 1./normcdf(Y).^(1/a) -1 ).^(a+2);
mesh(X,Y,Z)
figure
contour(X,Y,Z,100)
% disp('Hit any key to continue')
% pause
figure
p=2; n = 1000;
X = multi_BurrParetoLogistic(p, n, 'normal');
x1=X(:,1); x2=X(:,2);
z=(a+1)/a*normpdf(x1).*normpdf(x2)./(normcdf(x1).*normcdf(x2)).^(1/a+1)./ ...
    (1./normcdf(x1).^(1/a) + 1./normcdf(x2).^(1/a) -1 ).^(a+2);
xlin = linspace(-3,3,60);
ylin = linspace(-3,3,60);
[XX,YY] = meshgrid(xlin,ylin);
Z = griddata(x1,x2,z,XX,YY,'cubic');
meshc(XX,YY,Z) %interpolated
axis tight; hold on
plot3(x1,x2,z,'.','MarkerSize',15) 
zlim([-0.05, 0.2])
hold off



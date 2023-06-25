% This program make a function call to generate Khintchine distribution
% random numbers. The locations of these random numbers are plotted
% along with the corresponding distribution function.

% --- Draw graph for p=2 for confirmation

[X,Y]=meshgrid(-3:0.1:3);
Z=0.5-0.5*normcdf(max(abs(X),abs(Y)),0,1); % p. 157 (8.14) 
mesh(X,Y,Z)
disp('Hit any key to continue')
pause
p=2;
X=multi_Khintchine(p,1000, "normal");
x1=X(:,1); x2=X(:,2);
z=0.5-0.5*normcdf(max(abs(x1),abs(x2)),0,1);
xlin = linspace(-3,3,33);
ylin = linspace(-3,3,33);
[XX,YY] = meshgrid(xlin,ylin);
Z = griddata(x1,x2,z,XX,YY,'cubic');
meshc(XX,YY,Z) %interpolated
axis tight; hold on
plot3(x1,x2,z,'.','MarkerSize',15) 
hold off



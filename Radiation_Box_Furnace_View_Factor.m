% EJR 2016
% Determine diffuse view factor 
%  across 5x4x3 box furnace

numberOfPhotons = 30000;

X = 5; % Box Furnace Dimensions
Y = 4;
Z = 3; 

x1 = rand(numberOfPhotons,1).*X; % Start position
y1 = rand(numberOfPhotons,1).*Y;
 R = rand(numberOfPhotons,1); % (Random number, uniformly distributed from 0-1)
theta   = acos(sqrt(R));  % Emission angle! (To produce Uniform likelihood over unit hemisphere surface)
                          % To see why this gives uniform flux to a hemisphere, refer to: 
                          % Chen (2007) http://dx.doi.org/10.1103/PhysRevE.75.051504
                          % Or Rees (2015) http://dx.doi.org/10.1016/j.bpj.2015.09.023 
   phi  = 2*pi*rand(numberOfPhotons,1);

x3 = x1 + Z.*tan(theta).*cos(phi);
y3 = y1 + Z.*tan(theta).*sin(phi);

viewList13 = (x3>0).*(y3>0).*(x3<X).*(y3<Y);
viewFactor = sum(viewList13)/numberOfPhotons


viewFactorN = cumsum(viewList13)./(1:numberOfPhotons)' ;

figure(1)
plot(viewFactorN, 'lineWidth', 2)
xlabel('Number of "photons"', 'fontSize', 18)
ylabel('View Factor', 'fontSize', 18)
set(gca, 'fontSize', 18)
ylim([0 0.5])
grid on


%% 
myCol = zeros(length(viewList13), 3);
myCol(:,1) = (viewList13 == 1);
myCol(:,3) = (viewList13 == 0);

figure(3)
plot3([x1(1), x3(1)], [y1(1), y3(1)], [0 3], 'color' , myCol(1,:))
hold on
  for lp = 2:30
     plot3([x1(lp), x3(lp)], [y1(lp), y3(lp)], [0 3], 'color' , myCol(lp,:)) 
  end
  xlim([-2 5])
  ylim([-2 6])
  grid on
hold off
xlabel('y', 'fontSize', 16);
ylabel('x', 'fontSize', 16);
zlabel('z', 'fontSize', 16);

hold on
  plot3([0,0,X,X,0],[0,Y,Y,0,0], [0,0,0,0,0], 'k', 'lineWidth', 3);
  plot3([0,0,X,X,0],[0,Y,Y,0,0], [Z,Z,Z,Z,Z], 'r--', 'lineWidth', 2);
hold off
%%
% Exact method:

X = 4/3;
Y = 5/3;

F12 = (2/(pi*X*Y))*( (log( ((1+X^2)*(1+Y^2)/(1+X^2+Y^2))^0.5 )) ...
       + X*sqrt(1+Y^2)*atan(X/(sqrt(1+Y^2))) ... 
       + Y*sqrt(1+X^2)*atan(Y/(sqrt(1+X^2))) ...
       - X*atan(X) - Y*atan(Y) )
 

%%
% 2D opposing strips H = 3;
H=3;

F12 = sqrt(1+H^2) - H

numberOfPhotons = 30000;
x1 = rand(numberOfPhotons,1);

R = rand(numberOfPhotons,1);
   theta = acos(R);
   
 x2 = x1 + H*tan(R);
 
 viewList12 = (x2<1).*(x2>0);
 
mean(viewList12)
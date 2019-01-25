% EJR 2016
% Determine the view factor 
%  across a 5x4x3 box furnace
%  assuming a 'Lambertian' (diffuse) surface

numberOfPhotons = 30000;

X = 5; % Box Furnace Dimensions
Y = 4;
Z = 3; 

x1 = rand(numberOfPhotons,1).*X; % Start position, x-coordinate
y1 = rand(numberOfPhotons,1).*Y;
 R = rand(numberOfPhotons,1); % (Random number, uniformly distributed from 0-1)
theta   = acos(sqrt(R));  % Emission angle! (To produce Uniform likelihood over unit hemisphere surface)
                          % To see why this gives uniform flux to a hemisphere, refer to: 
                          % Chen (2007) http://dx.doi.org/10.1103/PhysRevE.75.051504
                          % Or Rees (2015) http://dx.doi.org/10.1016/j.bpj.2015.09.023 
													% 'theta' is the polar angle of emission, in spherical polar
													% coordinates.
													% Uniform likelihood of arrival over hemisphere
													% is the 'Lambertian' surface emission property.
   phi  = 2*pi*rand(numberOfPhotons,1); % Azimuthal angle of emission, uniformly distributed. 

x3 = x1 + Z.*tan(theta).*cos(phi); % x-position of arrival at opposite surface
y3 = y1 + Z.*tan(theta).*sin(phi);

viewList13 = (x3>0).*(y3>0).*(x3<X).*(y3<Y); % Logical test. Which photons hit the target surface
viewFactor = sum(viewList13)/numberOfPhotons % The proportion hitting the target surface is the view factor

% Proportion hitting target, as number of simulated photons increases:
viewFactorN = cumsum(viewList13)./(1:numberOfPhotons)' ; 

figure(1)
plot(viewFactorN, 'lineWidth', 2)
xlabel('Number of "photons"', 'fontSize', 18)
ylabel('View Factor', 'fontSize', 18)
set(gca, 'fontSize', 18)
ylim([0 0.5])
grid on


%% Plot a 3D visualisation of photon tracks
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
	title('Hits (red) and misses (blue)')
hold off
xlabel('y', 'fontSize', 16);
ylabel('x', 'fontSize', 16);
zlabel('z', 'fontSize', 16);

hold on
  plot3([0,0,X,X,0],[0,Y,Y,0,0], [0,0,0,0,0], 'k', 'lineWidth', 3);
  plot3([0,0,X,X,0],[0,Y,Y,0,0], [Z,Z,Z,Z,Z], 'r--', 'lineWidth', 2);
hold off
%% Alternative: compute the result of the (cosine rule) integral
% Exact method:

X = 4/3;
Y = 5/3;

F12 = (2/(pi*X*Y))*( (log( ((1+X^2)*(1+Y^2)/(1+X^2+Y^2))^0.5 )) ...
       + X*sqrt(1+Y^2)*atan(X/(sqrt(1+Y^2))) ... 
       + Y*sqrt(1+X^2)*atan(Y/(sqrt(1+X^2))) ...
       - X*atan(X) - Y*atan(Y) )
 

%% Get an extra result for the 2D case (surfaces infinitely long in the y-direction)
% 2D opposing strips H = 3, where H = Z/X
H=3;

% Algebraic result:
F12 = sqrt(1+H^2) - H

% Simulation result:
numberOfPhotons = 30000;
x1 = rand(numberOfPhotons,1);

R = rand(numberOfPhotons,1);
   theta = acos(R);
   
 x2 = x1 + H*tan(R);
 
 viewList12 = (x2<1).*(x2>0);
 
mean(viewList12)
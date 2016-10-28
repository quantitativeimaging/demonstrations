%% Example 14: Demonstration of Lagrangian minimisation principle
%   "What is the highest point within a given radius of a house on the side
%   of a hill?" 
%  - contrained optimisation if radius prevents reaching summit
%  - uncontstrained optimisation if the summit is inside the range

% Consider a hill of height (100 - x^2 - y^2)
xgv = -7:0.5:7;
ygv = -7:0.5:7;

[X,Y] = meshgrid(xgv,ygv);

Z = 100 - (X.^2 + Y.^2);

figure
mesh(X,Y,Z)

set(gca, 'fontSize', 18)
xlabel('x', 'fontSize', 18);
ylabel('y', 'fontSize', 18);
title('z=100 - (x^2 + y^2)', 'fontSize', 18);
caxis([50 100])

% 
hold on
scatter3(6,0, 70, 120, 'ko', 'filled') % Plot position of house

%
theta = 0:0.1:(2*pi+0.1);  % Plot a 4-mile radius around house (4 mile in x-y)
xRad = 4*cos(theta) +6;
yRad = 4*sin(theta);
zRad = 96*ones(length(theta),1);
plot3(xRad,yRad,zRad, 'lineWidth', 3)

% Now draw the contour of Z for the maximum possible z...

theta2 = 0:0.1:(2*pi+0.1);  % Plot a 4-mile radius around house (4 mile in x-y)
xRad2 = 2*cos(theta2);
yRad2 = 2*sin(theta2);
zRad2 = 96*ones(length(theta2),1);
plot3(xRad2,yRad2,zRad, 'r',  'lineWidth', 3)

% And a contour for a non-maximised Z:
theta3 = 0:0.1:(2*pi+0.1);  % Plot a 4-mile radius around house (4 mile in x-y)
xRad3 = 3*cos(theta3);
yRad3 = 3*sin(theta3);
zRad3 = 96*ones(length(theta3),1);
plot3(xRad3,yRad3,zRad, 'm',  'lineWidth', 3)

hold off

% Observe that the highest point on the 4-mile radius is 
% where the constraint lies parallel to the contour of equal height
% If so then we can't increase height by moving along the constraint
%   Whereas if not then we can
% This gives the principle of Lagrangian maximisation (or minimisation)
% Specifically, that contour of equal function value (height) lies parallel
%  to the line of the constraint at a stationary point.
% This fact is then expressed in terms of the normal to those lines being 
% parallel - i.e grad(f) is parallel to grad(g)


%%
%% Example 16: Carrots and Oatmeal
% http://stiglerdiet.com/blog/2012/Jan/09/carrots-oatmeal-operations-research/

% Plot graphs:
figure(10)
scatter([-10 60], [-1 2 ])
grid on
set(gca, 'fontSize', 18)
xlabel('Carrots (lbs)', 'fontSize', 18);
ylabel('Oatmeal (lbs)', 'fontSize', 18);
title('Possible combinations of food', 'fontSize', 18);
grid off
box on
  xlim([0 50])
  ylim([-0.05 1.6])
set(gca,'xTick',0:10:50)
set(gca,'yTick',0:0.5:1.5)
set(gcf,'color','w')
set(gca,'XMinorTick','on','YMinorTick','on')

xGrid = 0:0.05:50;
yGrid = 0:0.01:2;
[X,Y] = meshgrid(xGrid, yGrid);
Z = 0.5.*X + 0.5*Y;
Z((0.839.*X + 25.*Y) > 40) = 0;   % Too much fat
Z((172.*X + 1732.*Y) < 2000) = 0; % Not enough cals
Z((19.*X + 15.*Y) < 60) = 0; % Not enough Vit C

figure(11)
surf(X,Y,Z)
shading interp
view([0 90])
colormap(jet)
set(gca, 'fontSize', 18)
xlabel('Carrots (pounds)', 'fontSize', 18);
ylabel('Oatmeal (pounds)', 'fontSize', 18);
zlabel('Cost, GBP', 'fontSize', 18);
title('Linear programming', 'fontSize', 18);
%grid off
  xlim([-0.5 50])
  ylim([-0.1 1.6])
set(gca,'xTick',0:10:50)
set(gca,'yTick',0:0.5:1.5)
set(gcf,'color','w')
cmap = colormap(jet(256));
cmap(1,:) = [1,1,1];
colormap(cmap);
%colorbar()
% clabel('cost')

% Use a Matlab function to solve this linear programming problem:
% x = [carrots; oatmeal]

f = [0.5; 0.5]; % cost of carrots; oatmeal
A = [0.839  25   % Fat
    -172   -1732 % -ve calories...
    -19    -15]; % -ve vit C...
b = [40; -2000; -60 ];
lb = [0;0];

[x,fval,exitflag,output,lambda] = linprog(f,A,b,[],[],lb);

Carrots = x(1)
Oatmeal = x(2)

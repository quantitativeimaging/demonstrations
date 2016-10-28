%% Example 11: page 28 example
xgv = -1:0.1:5;
ygv = 1:0.1:6;

[X1,X2] = meshgrid(xgv,ygv);

Z = X1.^2 + X2.^2 - 4*X1 - 8*X2;

figure(7)
mesh(X1,X2,Z)
  view([30 60])

set(gca, 'fontSize', 18)
xlabel('X_1', 'fontSize', 18);
ylabel('X_2', 'fontSize', 18);
title('z = X_1^2 + X_2^2 - 4*X_1 - 8*X_2', 'fontSize', 18);

% Find minimum:
% Option 1: Exhaustive search of evaluated positions
[value, idx] = min(Z(:));
Z(idx)
X1(idx)
X2(idx)

% Option 2: Try a library function:
myFunction = @(r) r(1).^2 + r(2).^2 - 4*r(1) - 8*r(2); % Define function
r0 = [0,0]; % Initial guess of coords, "r" is an array for [x1, x2]

[x,fval] = fminsearch(myFunction, r0)
hold on
 scatter3(x(1), x(2), fval+0.5, 120, 'mo', 'filled')
hold off
% legend('mesh', 'fminsearch result')

%% Example 12: A harder problem 
% This is a maximisation of a function with oscillating ridges...
% Many methods will struggle with this type of problem...
xgv = -6:0.1:6;
ygv = -6:0.1:6;

[X1,X2] = meshgrid(xgv,ygv);

Z = (cos(X1+X2) + sin(2*X1) + sin(X2) + cos(4*X1 + 4*X2)).*...
     exp(-(X1.^2+X2.^2)/16 );
figure(7)
mesh(X1,X2,Z)
  view([30 60])
xlabel('X', 'fontSize', 18);
ylabel('Y', 'fontSize', 18);
title('z = (cos(x+y)+sin(2x)+sin(2y)+cos(4x+4y))*exp(-(x^2+y^2)/16)', ...
    'fontSize', 18);

  
% Define function:
myFunction = @(r)-1*(cos(r(1)+r(2))+sin(2*r(1))+sin(2*r(2))+ ...
                     cos(4*r(1)+4*r(2))).*exp(-(r(1).^2 + r(2).^2)/16); 
r0 = [0,0]; % Initial guess of coords, x(1) = x1, x(2) = x2

[x,fval] = fminsearch(myFunction, r0)
hold on
 scatter3(x(1), x(2), -fval, 120, 'mo', 'filled')
hold off
% legend('mesh', 'fminsearch result')

%%
% Example 13: Open tank volume:
Lgv = 0:0.02:3; % Plausible range of l
Wgv = 0:0.02:3;

[L,W] = meshgrid(Lgv,Wgv);

V = ( L.*W.*(10 - L.*W) )./(2*(L+W));

figure(8)
mesh(L,W,V)

set(gca, 'fontSize', 18)
xlabel('L', 'fontSize', 18);
ylabel('W', 'fontSize', 18);
zlabel('V', 'fontSize', 18);
% title('MATLAB meshgrid for open-top tank volume', 'fontSize', 18);

zlim([2 3.5])  % Make plot scale nice.
caxis([2 3.5])

[maxValue,index] = max(V(:)) % Find best Volume and co-ordinates

bestLength = L(index)
bestWidth  = W(index)
bestheight = (10 - bestLength*bestWidth)/(2* (bestLength+bestWidth))

% Alternatively, use fminsearch (of -f(x,y) to find maximum):
volFunction = @(r)  -(r(1).*r(2).*(10 - r(1).*r(2)))./(2*(r(1)+r(2)));
r0 = [2,2];

[r,fval] = fminsearch(volFunction, r0);
r
l = r(1)
w = r(2)
Vol   = -fval % minus sign  is for using fminsearch to find maximum...

hold on
 scatter3(l,w,Vol,  120, 'mo', 'filled')
hold off

%
%% Check Hessian with symbolic algebra:
syms x y
V = x*y*(10-x*y)/(2*(x+y));

H = hessian(V, [x y])

x = l;
y = w;
% Then run H at the command line
H_values = eval(H)
det(H_values)
% Positive determinant, with negative second deriv.s implies maximum. 

% Note: long-hand, Hessian should be a symmetric matrix with:
num  = 2*(l + w)^2*(-2*w^2*l-2*w^3) - w^2*(10-l^2-2*l*w)*4*(l+w);
num2 = 2*(l+w)^2*(20*w - 2*w*l^2 - 6*l*w^2) - w^2*(10-l^2-2*l*w)*4*(l+w);
den  = 4*(l+w)^4;
f11 = num/den 
f12 = num2/den


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
title('MATLAB meshgrid for a hill', 'fontSize', 18);
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
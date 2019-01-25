% Matlab examples for Engineering Maths, CET IB, 2016
% Eric Rees
% 


% To use these scripts:
%  1. Start Matlab
%  2. Enter these scripts on the console
%      one line at a time
%  3. Or open this file in the editor 
%      and run entire sections (Ctrl + Enter)


% EXAMPLE 1: Trapezium rule errors
%  Integrate exp(x) dx from 0 to 1
% 
clear    % Clear workspace

n = 5;   % Number of steps
h = 1/n; % Interval

x = 0:h:1;
f = exp(x);

Itrapezium = 0.5*h*(f(1)+2*sum(f(2:end-1))+f(end));
Iexact     = exp(1) - 1;
difference = Itrapezium - Iexact;


% EXAMPLE 1b: Plot graph for trapezium rule errors for several values of 

listOfN          = zeros(100,1); % Pre-allocate an empty list 
listOfDifference = zeros(100,1);
for loop = 1:100
  n = loop;
  h = 1/n;   % Interval
  x = 0:h:1;
  f = exp(x);
  Itrapezium = 0.5*h* (f(1) + 2*sum(f(2:end-1)) + f(end));
  Iexact     = exp(1) - 1;
  difference = Itrapezium - Iexact;
  listOfN(loop)          = n;
  listOfDifference(loop) = difference;
end

figure(1)
 scatter(listOfN(1:10), listOfDifference(1:10));
 set(gca, 'fontSize', 16)
 xlabel('n', 'fontSize', 18);
 ylabel('Error', 'fontSize', 18);
 title('Trapezium Rule Error', 'fontSize', 18);
 set(gcf,'Position',[100,100,400,300]);
 set(gcf,'color','w')
figure(2)
 scatter(log10(listOfN), log10(listOfDifference));
 set(gca, 'fontSize', 16)
 xlabel('log_{10} (number of steps)', 'fontSize', 18);
 ylabel('log_{10} (Total error)', 'fontSize', 18);
 title('Trapezium Rule Error', 'fontSize', 18);
 set(gcf,'Position',[100,100,500,400]);
 grid on
 legend('Trapezium rule')
 
 
%% 
%  EXAMPLE 2: Simpson's rule integration
%  Also for  exp(x) dx from 0 to 1 
clear

n = 2;     % Consider n intervals. n is even.
h = 1/n;   % The width of each interval is h
                   
x = 0:h:1;    % x coordinates 
y = exp(x);   % y coordinates

Isimpson = (1/3)*h*(y(1) + 4*sum( y(2:2:end-1) ) + ...
                    2*sum( y(3:2:end-2) ) +y(end) );
Iexact = exp(1) - 1;
difference = Isimpson - Iexact;

% EXAMPLE 2b:
%  Compare Simpson's rule error with Trapezium rule:
%  For several n-values (must be even):
listOfN          = zeros(50,1);
listOfDiffTrap   = zeros(50,1);
listOfDiffSimp   = zeros(50,1);
for loop = 1:50
  n = 2*loop;
  h = 1/n;     % Interval
  x = 0:h:1;
  y = exp(x);
  Itrapezium = 0.5*h* (y(1) + 2*sum(y(2:end-1)) + y(end));
  Isimpson   = (1/3)*h*(y(1) + 4*sum( y(2:2:end-1) ) + ...
                        2*sum( y(3:2:end-2) ) +y(end) );
  
  Iexact     = exp(1) - 1;
  
  listOfN(loop)          = n;
  listOfDiffTrap(loop) = Itrapezium - Iexact;
  listOfDiffSimp(loop) = abs( Isimpson - Iexact );
end

figure(3)
 scatter(log10(listOfN), log10(listOfDiffTrap), 'bo');
 hold on
  scatter(log10(listOfN), log10(listOfDiffSimp), 'rx');
 hold off
 set(gca, 'fontSize', 16)
 xlabel('log_{10} n', 'fontSize', 18);
 ylabel('log_{10} (Total Error)', 'fontSize', 18);
 title('', 'fontSize', 18);
 set(gcf,'Position',[100,100,400,300]);
 set(gcf,'color','w')
 legend('Trapezium','Simpson','fontSize',14)
 grid on
ylim([-10 0])


%%
% Example 3: Noise Tolerance of Trapezium and Simpson's Rules
%            integrate  exp dx  from x = 0 to 1 with noise

clear
n = 10;
h=1/n;
xData = 0:h:1; 

nRepeats = 1000;   % Run this 1000 times 

Isimpson = zeros(length(nRepeats),1); % Pre-allocate memory
Itrap =    zeros(length(nRepeats),1); % 

for loop = 1:nRepeats
  yData = exp(xData) + 0.5*randn(1, length(xData)); 
  Isimpson(loop) = (1/3) * h * (yData(1) + 4*sum(yData(2:2:end-1)) + ...
                          2*sum(yData(3:2:end-2)) +yData(end) );
  Itrap(loop)    = 0.5*h*(yData(1)+2*sum(yData(2:end-1))+yData(end));
end

Iexact = exp(1) - 1;

errorOfSimpson = abs( Isimpson - Iexact) ;
errorOfTrap    = abs( Itrap    - Iexact) ;
mean(errorOfSimpson)
mean(errorOfTrap)

% Bar chart for mean absolute errors:
figure(4)
 bar( [mean(errorOfSimpson) , mean(errorOfTrap)] )
 set(gca,'XTickLabel',{'Simpson', 'Trapezium'})
 set(gca, 'fontSize', 16)
 title('Mean absolute error of noisy I = e^x dx', 'fontSize', 18);
colormap summer
ylabel('Mean Absolute Error', 'fontSize', 16);

% Scatterplot for absolute errors (hardly enlightening)
figure(5)
 scatter((1:length(errorOfSimpson)), errorOfSimpson, 'rx')
 hold on
 scatter((1:length(errorOfTrap)), errorOfTrap, 'bo')
 hold off
legend('Simpson','Trapezium')
set(gca, 'fontSize', 16)
title('Total error, with noise', 'fontSize', 18);
xlabel('Test number', 'fontSize', 16);
ylabel('Absolute error', 'fontSize', 16);


%% 
%  Example 4: Euler and Modified Euler for ODEs
% Exponential decay:
% 
% Note:
%   Interesting cases are (h = 0.01, 0.1, 0.15, 0.2, 0.203)
%    Note the axis limits, xlim() and ylim(), may hide 
%    the divergence to infinity for large h here

steps = 200;    % Number of steps
h = 0.01;       % Step size

t = 0;    % Initial time
y = 1;    % Initial amount

listT = zeros(steps,1);
listY = zeros(steps,1);
for loop = 1:steps
   listT(loop) = (loop-1)*h;
   listY(loop) = y;
     
    f = -10*y;    % Slope
    y = y + h*f; 
    t = t+h; 
end
figure(6)
 plot(listT, listY, 'lineWidth', 2);

 set(gca, 'fontSize', 16)
 xlabel('time', 'fontSize', 18);
 ylabel('N', 'fontSize', 18);
 title('Exponential decay simulation', 'fontSize', 18);
 set(gcf,'Position',[100,100,400,300]);
 set(gcf,'color','w')
 xlim([0 1])
 
 
%% 
%  Example 5:  Euler and its errors

n = 5;    % Set step size and number of steps
h = 0.1;
x = 0;    % Set initial values
y = 1;

listX = zeros(n,1);
listY = zeros(n,1);

for loop = 1:n
   listX(loop) = x;  % Store results
   listY(loop) = y;
    
    k1 = x+y;         % Slope at start of  step
    
    y = y + h*( k1 );
    x = x+h;
end

yExact = 2*exp(x) - x - 1;
yError = yExact - y;

figure(7)

plot(listX, 2*exp(listX) - listX - 1, 'r','lineWidth', 4 )
hold on
   plot(listX, listY, 'lineWidth', 2); % Sim
   scatter(listX, listY, 200, 'bo')
hold off
set(gca, 'fontSize', 14)
xlabel('x', 'fontSize', 14);
ylabel('y', 'fontSize', 14);
% title('Euler', 'fontSize', 14);
legend('exact', 'Euler')
set(gcf,'Position',[100,100,400,300]);
set(gcf,'color','w')


%%
% Example 6: Modified Euler and its errors
n = 5;    % Set step size and number of steps
h = 0.1;
x = 0;    % Set initial values
y = 1;

listX = zeros(n,1);
listY = zeros(n,1);

for loop = 1:n
   listX(loop) = x;  % Store results
   listY(loop) = y;
    
    k1 = x+y;         % Slope at start of  step
    yG = y + h*k1;    % Guess y at end of step
    k2 = (x+h) + yG;  % Guess slope at end of step
    
    y = y + 0.5*h*( k1 + k2 );
    x = x+h;
end

yExact = 2*exp(x) - x - 1;
yError = yExact - y;

figure(8)
 hold on
   plot(listX, 2*exp(listX) - listX - 1, 'r','lineWidth', 10 )
   plot(listX, listY, 'lineWidth', 4); % Sim
   scatter(listX, listY, 200, 'bo')
 hold off
 set(gca, 'fontSize', 16)
 xlabel('X', 'fontSize', 18);
 ylabel('Y', 'fontSize', 18);
 title('Modified Euler', 'fontSize', 18);
 legend('exact', 'simulated', 'simulated')
 set(gcf,'Position',[100,100,400,300]);
 set(gcf,'color','w')
 
 
 %%
 %  Example 7: RK4 method and its errors
 
n = 5;   % n steps
h = 0.1;   % step size (large)
listX = zeros(n,1);
listY = zeros(n,1);

x = 0;   % initial x
y = 1;   % initial y
for loop = 1:n
  listX(loop) = x;
  listY(loop) = y;
  
  k1 = x+y;
  k2 = (x + 0.5*h) + (y + 0.5*h*k1);
  k3 = (x + 0.5*h) + (y + 0.5*h*k2);
  k4 = (x + h)     + (y + h*k3);
  
  averageSlopeEstimate = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
  
  y = y + h*averageSlopeEstimate;
  x = x + h;
end

yExact = 2*exp(x) - x - 1;
yError = y - yExact;

figure(9)
plot(listX, 2*exp(listX) - listX - 1, 'r','lineWidth', 10 )
 hold on
   plot(listX, listY, 'b', 'lineWidth', 3)
   scatter(listX, listY, 200, 'bo')
 hold off
 set(gca, 'fontSize', 16)
 xlabel('X', 'fontSize', 18);
 ylabel('Y', 'fontSize', 18);
 title('RK4', 'fontSize', 18);
 legend('exact', 'simulated', 'simulated')
set(gcf,'Position',[100,100,400,300]);
set(gcf,'color','w')


%%
% Example 8: RK4 comparison with Euler, Mod Euler

n = 7;   % n steps
h = 1;   % step size (1 is 'large')
x = 0;   % initial x
y = 1;   % initial y

listX_euler    = zeros(n+1,1); % Prepare lists to store results
listY_euler    = zeros(n+1,1);
listX_euler(1) = x;
listY_euler(1) = y;

listX_modEuler = zeros(n+1,1);
listY_modEuler = zeros(n+1,1);
listX_modEuler(1) = x;
listY_modEuler(1) = y;

listX_RK4      = zeros(n+1,1);
listY_RK4      = zeros(n+1,1);
listX_RK4(1)   = x;
listY_RK4(1)   = y;

for loop = 1:n
   
   % Euler:
   k1  = (listX_euler(loop) + listY_euler(loop));
   listX_euler(loop+1) = listX_euler(loop) + h;
   listY_euler(loop+1) = listY_euler(loop) + h*k1;
   
   % Modified Euler (or RK2):
   k1  = (listX_modEuler(loop) + listY_modEuler(loop));
   yG  = listY_modEuler(loop) + k1*h;
   k2  = (listX_modEuler(loop) + h   + yG);
   listX_modEuler(loop+1) = listX_modEuler(loop) + h;
   listY_modEuler(loop+1) = listY_modEuler(loop) + (1/2)*h*(k1+k2);
    
   % 4th order Runge-Kutta
   k1 = (listX_RK4(loop) + listY_RK4(loop));
   k2 = (listX_RK4(loop) + 0.5*h) + (listY_RK4(loop) + 0.5*h*k1);
   k3 = (listX_RK4(loop) + 0.5*h) + (listY_RK4(loop) + 0.5*h*k2);
   k4 = (listX_RK4(loop) + h)     + (listY_RK4(loop) + h*k3);
    
   averageSlopeEstimate = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
   listX_RK4(loop+1) = listX_RK4(loop) + h;
   listY_RK4(loop+1) = listY_RK4(loop) + h*averageSlopeEstimate;

end


figure(10)
 xx = 0:0.02:(n*h);     % x positions
 yy = 2*exp(xx) - xx - 1; % Exact y
 plot(xx,yy, 'r','lineWidth', 4 ) % Exact
 hold on
   plot(listX_euler, listY_euler, 'b--','lineWidth', 4); % Sim
   plot(listX_modEuler, listY_modEuler, 'g--','lineWidth', 4); % Sim
   plot(listX_RK4, listY_RK4, 'k--','lineWidth', 4); % Sim
 hold off
 set(gca, 'fontSize', 16)
 xlabel('X', 'fontSize', 18);
 ylabel('Y', 'fontSize', 18);
 title('', 'fontSize', 18);
 legend('Exact','Euler', 'Mod Euler', 'RK4')
% set(gcf,'Position',[100,100,400,300]);
% set(gcf,'color','w')


%%
% Example 9: 
%   Accuracy of RK 4 versus Euler and Mod Euler, for different n values:
%

nValues = [1,2,4,10,100];

for loopNs = 1:length(nValues)
  n = nValues(loopNs);   % n steps
  h = 0.1/n;            % step size for total distance 0.5

  x = 0;   % initial x
  y = 1;   % initial y for RK4
  yEu = 1; % Initial y for Euler
  yME = 1; % Initial y for Modified Euler
  
  for loop = 1:n
    k1 = x+y;                 % RK4 method 
    k2 = (x + 0.5*h) + (y + 0.5*h*k1);
    k3 = (x + 0.5*h) + (y + 0.5*h*k2);
    k4 = (x + h)     + (y + h*k3);
    
    averageSlopeEstimate = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    y = y + h*averageSlopeEstimate;
    
    yEu = yEu + h*(x + yEu); % Euler method
    
    k1ME = x+yME;            % Modified Euler method
    yG   = yME + k1ME*h; 
    k2ME = x+h + yG;
    yME  = yME + 0.5*h*(k1ME + k2ME);
    
    x = x + h; 
  end
  
  yRK4(loopNs)   = y;
  yEuler(loopNs) = yEu;
  yModE(loopNs)  = yME;
end
yExact = 2*exp(x) - x - 1;


log_yErrorRK4   = log10( abs(yExact - yRK4)   );
log_yErrorEuler = log10( abs(yExact - yEuler) );
log_yErrorME    = log10( abs(yExact - yModE)    );

logNvalues = log10(nValues);


figure(11)
plot( logNvalues, log_yErrorEuler, 'b', 'lineWidth', 3)
hold on % Overlay some extra plots
  plot(logNvalues, log_yErrorRK4, 'k', 'lineWidth', 3)
  plot(logNvalues, log_yErrorME,  'r', 'lineWidth', 3)
  legend('Euler', 'RK4', 'Mod Euler', 18, 'Location', 'SouthWest');
  scatter( logNvalues, log_yErrorRK4, 70, 'k', 'filled') % Show points
  scatter( logNvalues, log_yErrorEuler, 70, 'b', 'filled') % Show points
  scatter( logNvalues, log_yErrorME, 70, 'r', 'filled') % Show points
hold off
% Make the plot readable by increasing font sizes etc.
set(gca, 'fontSize', 18)
xlabel('log_{10} N', 'fontSize', 18);
ylabel('log_{10} Error', 'fontSize', 18);
title('Runge-Kutta global error', 'fontSize', 18);
grid on


slopeRK = (log_yErrorRK4(end)- log_yErrorRK4(end-1))/(logNvalues(end) - logNvalues(end-1));
disp(slopeRK);


%%
%  Example 10: Euler for a 2nd order ODE (simple harmonic oscillator)

h = 0.1; % Time interval: try 0.02, 0.2
y = 0;   % Initial y value (position)
v = 1;   % Initial dy/dt   (velocity)
t = 0;
for loop = 1:500 
    acceleration = -y; % Slopes
    velocity     = v;
    
    v = v + h*acceleration;
    y = y + h*velocity;
    t = t+h;
end

% Repeat, plotting results:
clear
n = 500;
h = 0.1;  % Time interval (try 0.02, 0.05, 0.2) 

y = 0;    % Initial conditions
v = 1;
t = 0;

listY = zeros(n,1);
for loop = 1:n 
    acceleration = -y; % Slopes
    velocity     = v;
    
    v = v + h*acceleration;
    y = y + h*velocity;
    t = t+h;
    
    listY(loop) = y;
end

xx = (1:n)*h; % time co-ordinates

figure(12)
scatter(xx, listY) % Plot simulation result
hold on
  plot(xx, sin(xx),'r', 'lineWidth', 2) % Exact
hold off
 legend('Euler, 2nd order ', 'Exact')
 xlim([0 10])
 ylim([-2 4])

 set(gca, 'fontSize', 16)
 xlabel(' x', 'fontSize', 18);
 ylabel(' y', 'fontSize', 18);
 title('', 'fontSize', 18);
 set(gcf,'Position',[100,100,400,300]);
set(gcf,'color','w')
 

%%
% Example 11:
% Use Modified Euler for the simultaneous ODEs in Example Sheet 1, Q.6

clear
% (1): INPUT
h = 0.2;            % interval in time, 0.2 seconds
numberOfSteps = 3; 

t = 0;              % initial value of independent variable
x = 0;             % initial value of x
y = 1;             % initial value of y

% Pre-allocate computer memory to store outputs 
% (This is unnecessary if you don't need to recall the iterating values)
% (It is also unnecessary, here, as the computational time is brief anyway)
listT  = zeros( numberOfSteps + 1 , 1 );
listX  = zeros( numberOfSteps + 1 , 1 ); 
listY  = zeros( numberOfSteps + 1 , 1 );

listK1 = zeros( numberOfSteps + 1 , 1 ); % K1 is the initial slope for X
listM1 = zeros( numberOfSteps + 1 , 1 ); % M1 is the initial slope for Y
listK2 = zeros( numberOfSteps + 1 , 1 );
listM2 = zeros( numberOfSteps + 1 , 1 );
listXg = zeros( numberOfSteps + 1 , 1 ); % First guess of next X
listYg = zeros( numberOfSteps + 1 , 1 );

listT(1) = t;
listX(1) = x;
listY(1) = y;

% (2.) PROCESS
%
for n = 1:numberOfSteps
   
    k1 = x*y + t;
    m1 = x - t;
    xg = x + h*k1;
    yg = y + h*m1;
    k2 = xg*yg + (t + h);
    m2 = xg - (t + h);    
    
    t = t + h;
    x = x + 0.5 * h * (k1 + k2);
    y = y + 0.5 * h * (m1 + m2);
    

    listT(n+1) = t;
    listX(n+1) = x;
    listY(n+1) = y;

    listK1(n)  = k1;
    listXg(n)  = xg;
    listK2(n)  = k2;
    listM1(n)  = m1;
    listYg(n)  = yg;
    listM2(n)  = m2;
end

figure(13)
plot(listT, listX, 'lineWidth', 2);
hold on
  plot(listT, listY, 'r','lineWidth', 2 )
hold off
legend('x','y');

%% 
%  Example 12
%    The Lokta-Volterra predator prey model
%    (With populations as a continuous variable)
% 

listT = zeros(n,1);
listX = zeros(n,1);
listY = zeros(n,1);

a = 0.4;
b = 0.02;
c = 0.005;
d = 0.1;

n = 10000;
h = 0.05;

t = 0;   % time
x = 80;  % Baboons
y = 40;  % Cheetahs

for loop = 1:n
    listT(loop) = t;
    listX(loop) = x;
    listY(loop) = y;
    
    f = a*x - b*x*y;
    g = c*x*y - d*y;
    
    x = x + f*h;
    y = y + g*h;
    t = t+h;
end

figure(14)
  plot(listT, listX, 'b', 'lineWidth', 2)
  hold on
    plot(listT, listY, 'r', 'lineWidth', 2);
  hold off
set(gca, 'fontSize', 18)
xlabel('Time', 'fontSize', 18);
ylabel('Number', 'fontSize', 18);
legend('Baboons (prey)','Cheetahs (predator)', 18);
title('A predator prey model', 'fontSize', 18);
grid on
ylim([0 250])
  

figure(15)
  plot(listX, listY, 'b', 'lineWidth', 2)
set(gca, 'fontSize', 18)
xlabel('Baboons', 'fontSize', 18);
ylabel('Cheetahs', 'fontSize', 18);
legend('', 18);
title('A predator prey model', 'fontSize', 18);
grid on

% Modified Euler, same problem:
listT = zeros(n,1);
listX = zeros(n,1);
listY = zeros(n,1);
t = 0;
x = 80;
y = 40;
for loop = 1:n
    listT(loop) = t;
    listX(loop) = x;
    listY(loop) = y;
    
    k1 = a*x - b*x*y;
    m1 = c*x*y - d*y;
    xG = x + k1*h;
    yG = y + m1*h;
    k2 = a*xG - b*xG*yG;
    m2 = c*xG*yG - d*y;
    
    x = x + 0.5*h*(k1+k2);
    y = y + 0.5*h*(m1+m2);
    t = t+h;
end
figure(16)
  plot(listT, listX, 'b', 'lineWidth', 2)
  hold on
    plot(listT, listY, 'r', 'lineWidth', 2);
  hold off
set(gca, 'fontSize', 18)
xlabel('Time', 'fontSize', 18);
ylabel('Number', 'fontSize', 18);
legend('Baboons (prey)','Cheetahs (predator)', 18);
title('Mod Euler', 'fontSize', 18);
ylim([0 150])
grid on


%%
% Example 13, Truncation error
% Demonstrate introduction of wrong solution into y'' = y, y(0)=1, y'(0)=-1
% Given initial conditions, solution should be e^(-x)
% But some e^x can be added due to truncation or rounding error
% This demonstration introduces a (very very bad!) rounding error in order 
% to make the point that numerical errors can introduce wrong parts of the 
% general ODE solution into specific cases. 

n = 100;    % Set step size and number of steps
h = 0.1;
x = 0;    % Set initial values
y = 1;
v = -1;

listX = zeros(n,1);
listY = zeros(n,1);

for loop = 1:n
   listX(loop) = x;  % Store results
   listY(loop) = y;
    
   dvBydx = y;         % Slope at start of  step
    
   y = y + h*( v );
   v = v + h*dvBydx;
   x = x+h;
   
   y = ceil(y*10)/10; % Artificially put in bad numerical rounding error.
	 % Comment out the above line to produce a good simualtion. 
end


figure(17)

plot(listX, listY, 'lineWidth', 2); % Sim
hold on
    scatter(listX, listY, 200, 'bo')
hold off
set(gca, 'fontSize', 14)
xlabel('x', 'fontSize', 14);
ylabel('y', 'fontSize', 14);
title('Euler, large truncation error', 'fontSize', 14);
legend('Euler')
set(gcf,'Position',[100,100,400,300]);
set(gcf,'color','w')
%

%% Engineering Maths Matlab Examples 
% 

% Example 13b: Biexponential decay
% This example demonstrates a stiff system of 2 ODEs
%  in which two processes with different timescales mean a simple Euler 
%  method simulation would need a large number of steps to return an 
%  accurate result.
% This particular example is just Y = Ya + Yb, and could be solved by 
%  adding up the separate exponentual decays Ya and Yb. The problem is that
%  some systems have the type of multi-timescale processes that
%  cannot be separated. And these require more advanced numerical methods
%  such as the adaptive step size...
%

h     = 0.02;
steps = 1000;
k1 = 10;
k2 = 0.1;
t = 0;
y = 2;

listT = zeros(steps,1);
listY = zeros(steps,1);

for loop = 1:steps
  listT(loop) = t;
  listY(loop) = y;
    
  f = -k1*exp(-k1*t) - k2*exp(-k2*t); % Pretend we can't solve this exactly
  y = y + h*f;
  t = t+h;
end

figure(1)
 hold on
   scatter(listT, listY, 'bo')
 hold off
 set(gca, 'fontSize', 16)
 xlabel('t', 'fontSize', 18);
 ylabel('y', 'fontSize', 18);
 title('Biexponential decay', 'fontSize', 18);
 % legend('exact', 'simulated', 'simulated')
set(gcf,'Position',[100,100,400,300]);
set(gcf,'color','w')


%% With adaptive step size for x:
h     = 0.02;
steps = 40;
k1 = 10;
k2 = 0.1;
t = 0;
y = 2;

listT = zeros(steps,1);
listY = zeros(steps,1);

% Try fixed step in y...
for loop = 1:steps
  listT(loop) = t;
  listY(loop) = y;
    
  f = -k1*exp(-k1*t) - k2*exp(-k2*t);
  h = abs(0.05/f);
  y = y + h*f;
  t = t+h;
end

figure(2)
 hold on
   scatter(listT, listY, 'bo')
 hold off
 set(gca, 'fontSize', 16)
 xlabel('t', 'fontSize', 18);
 ylabel('y', 'fontSize', 18);
 title('Biexponential decay', 'fontSize', 18);
 % legend('exact', 'simulated', 'simulated')
set(gcf,'Position',[100,100,400,300]);
set(gcf,'color','w')


%% With adaptive step size (based on a numerical estimate of error)
clear % avoid memory buildup due to loop below
h     = 0.02;
distance = 20;
k1 = 10;
k2 = 0.1;
t = 0;
y = 2;

loop = 1;

targetError = 0.005;

% Try to keep proportional error 1%
while t < distance
  listT(loop) = t; 
  listY(loop) = y;
  
  % 1 step of Euler simulation
  f = -k1*exp(-k1*t) - k2*exp(-k2*t);
  y1step = y + h*f;

  % 2 steps of Euler simulation
  f1 = -k1*exp(-k1*t) - k2*exp(-k2*t);
  yG = y + f1*0.5*h;
  f2 = -k1*exp(-k1*(t+0.5*h)) - k2*exp(-k2*(t+0.5*h));
  y2step = yG + f2*0.5*h;
  
  localErrorEst = abs(y2step - y1step)/2;
  localPropError = localErrorEst/y;
  
  if(localPropError > targetError)
      h = h/2;
      continue;
  elseif(localPropError < targetError/4)
      h = 2*h;
      continue;
  end
  
  t = t+h;
  y = y2step;
  loop = loop+1;
end
  listT(loop) = t;
  listY(loop) = y;
  
figure(3)
   scatter(listT, listY, 'bo')

 set(gca, 'fontSize', 16)
 xlabel('t', 'fontSize', 18);
 ylabel('y', 'fontSize', 18);
 title('Biexponential decay', 'fontSize', 18);
 % legend('exact', 'simulated', 'simulated')
set(gcf,'Position',[100,100,400,300]);
set(gcf,'color','w')
ylim([0 2])
xlim([0 20])


% Page 45 calculation

ej  = 0.8;
sig = 5.67E-8;
h   = 100;
Tg = 800; % Kelvin
Tw = 500; % Kelvin

Temperature = roots([ ej*sig, 0, 0, h, -(h*Tg + sig*ej*Tw^4) ])

% ans =
% 
%   711.8660

%
Tdata = 680:1:740;

yData = h*(Tg - Tdata) - ej*sig*(Tdata.^4 - Tw^4);

figure(1)
plot(Tdata, yData, 'lineWidth', 2)
grid on
xlabel('T_j / K', 'fontSize', 16)
ylabel('f(T_j) / K', 'fontSize', 16)
set(gca, 'fontSize', 14)
%% 
 
 Q = (63 + 2.65 + 55)^(-1) * 5.67E-8 * (294^4 - 77^4)
 
 
 %%

% Solve simultaneous nonlinear equations:
% Q12 -  Qg1         = 0
% Q23 - (Q12 + Qg12) = 0; 
 
e1  = 0.8;
e2  = 0.4;
hg1 = 100;
hg2 = 100;
d1  = 2;
d2  = 10;
sig = 5.67E-8;
Tg = 800;
Tw = 500;

F12 = pi*d1/( (1/e1) + (d1/d2)*(1/e2 - 1) );
F23 = pi*d2* e2;

Eq = @(T) [ F12*sig*( T(1)^4-T(2)^4 ) - pi*d1*hg1*( Tg -T(1) ),  ...
            -F12*sig*( T(1)^4-T(2)^4 ) + ...
            F23*sig*( T(2)^4-Tw^4 )    - ...
            2*pi*d2*hg2*(Tg - T(2)) ];

Temperatures = fsolve(Eq, [800 800])


  
  
% Regression to the mean demonstration
% 26/4/2013

N = 1000;              % Number of samples
trueData = randn(N,1); % Values of samples. Mean = zero, Variance = 1.

% xx are exact obsevations in this example: zero measurement error

% yy are observations with error - e.g. from a measurement error  
%
% Coefficients for generating yy in this sample are chosen to ensure:
% (a) yy and xx have the same mean 
% (b) yy and xx have the same variance
% 
% (a) can be achieved by de-biasing the measurements, or may occur without
% de-biasing if the observation error is a "Zero Mean Error"
% (b) can be achieved by rescaling the measurements
%    - this is legitimate, because we are interested in the change in the 
%      relative position of 2 measurements of the same specimen, and not
%      the absolute values. 
%
% In this example the mean of the data (both true and observed) is zero:
% Regression towards the Mean still occurs if the mean is non-zero
% I.e. the effect is not "regression towards (an arbitary) zero"
% Prove this for yourself with an edited version of this MATLAB script
% (Unless this is obvious to you; which seems unlikely.)

xx = trueData;                      
yy = 0.8*trueData + 0.6*randn(N,1); % Parameters lead to the same variance
                            % Because std deviations add in quadrature
                            % And 0.64 + 0.36 = 1
mean(xx);
mean(yy);
std(xx);
std(yy);  % Means should converge to zero, and stds should converge to 1

p = polyfit(xx,yy,1);  % p returns 2 coefficients fitting r = a_1 * x + a_2

% Observe that the linear fit has a slope of about 0.8, and not 1.
% This is due to Regression Towards The Mean

p2 = polyfit(yy,xx,1); % REMOVE the semicolon for demo code

% Observe that the linear fit is again 0.8, due to regression to the mean

myPlot1=figure(1);
scatter(xx,yy,'+k')

hold on 
simFromXX = floor(min(xx)):ceil(max(xx));
plot(simFromXX,p(1)*simFromXX+p(2),'g','LineWidth',2);

simFromYY = floor(min(yy)):ceil(max(yy));
plot(p2(1)*simFromYY+p2(2), simFromYY,'r','LineWidth',2);

axis equal
plot([-4,4],[-4,4],'color',[0.7 0.7 0.7],'LineWidth',2)

legend('Measurements','Expected Y given X', ...
    'Expected X given Y', 'Y=X', ...
    'Location','SouthEast');
hold off

xlabel('Observation 1', 'fontsize',18)
ylabel('Observation 2 ', 'fontsize',18)
title('Regression To The Mean', 'fontsize',18)

set(gca,'FontSize',18,'fontweight','bold');
set(myPlot1,'Position',[100,100,720,600]); % 720 px wide, 600 high
set(myPlot1,'color','w')

xlim([-4.5 4.5])
ylim([-4.5 4.5])
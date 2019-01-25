% Corrosion Course Matlab Figures
% EJR, 2015, and 2016 update for axis labels.
% 
% Changes from 2014-15 notes
% 1. Use the 'Electrode potential on Y-axis' convention for Evans diagrams
% 2. Add kinetic plot to illustrate ion resistance (salt water vs pure
% water corrosion)
%
%
%% 1 October 2014
%  Plot Pourbaix pH-potential diagram for tracing

figure(1)
 pHdata       = [0, 10];
 E_Hydrogen   = [0, -0.59];
 E_Oxygen     = [1.22, 1.22-0.59];
 E_Cu_CuIons  = [0.34, 0.34];
 E_Cu_CuO     = [0.57, 0.57-0.59];
 
 pH_CuIons_CuO = [7.9/2, 7.9/2];
 E_CuIons_CuO  = [0.34 1.5];
 
  hold on
  plot(pHdata, E_Hydrogen, 'k--', 'lineWidth', 2);
  plot(pHdata, E_Oxygen,   'b--', 'lineWidth', 2);
  
  xlim([0 10])
  ylim([-0.6 1.5])
  grid on
  set(gca,'FontSize',18,'fontweight','bold');
  set(1,'Position',[100,100,800,600]);
  set(gcf,'color','w')
  xlabel('pH','fontSize',18)
  ylabel('E, V_{SHE}','fontSize',18)

  plot(pHdata, E_Cu_CuIons,'r-',  'lineWidth', 4);
  plot(pHdata, E_Cu_CuO,'g--',  'lineWidth', 4);
  plot(pH_CuIons_CuO, E_CuIons_CuO, 'c-',  'lineWidth', 4);
 
 title('pH-potential diagram','fontSize',18)
 
 hold off

%% 1 October 2015
%  Plot Pourbaix pH-potential diagram for Zinc and Iron (simplified)

figure(11)
 pHdata       = [0, 14];
 E_Hydrogen   = [0, -0.059*14];
 E_Oxygen     = [1.22, 1.22-0.059*14];
 E_Zn_Zn2p    = [-0.937, -0.937];
 E_Zn_ZnOH2   = [-0.439 - 8.44*0.059, -0.439 - 10.75*0.059];
 E_Zn_ZnO2_2m = [0.441 - 0.1182*10.75 - 8*0.0295, 0.441 - 0.1182*14 - 8*0.0295];
 E_Zn2p_ZnOH2 = [-0.937, 1.5];
 E_ZnOH2_ZnO2_2m=[-0.439 - 10.75*0.059, 1.5];
 
pH_Zn_Zn2p   = [0, 8.44];
pH_Zn_ZnOH2  = [8.44, 10.75];
pH_Zn_ZnO2_2m = [10.75, 14];
pH_Zn2p_ZnOH2 = [8.44, 8.44];
pH_ZnOH2_ZnO2_2m = [10.75, 10.75];


  plot(pHdata, E_Hydrogen, 'k--', 'lineWidth', 1);
hold on
  plot(pHdata, E_Oxygen,   'b--', 'lineWidth', 1);
  
  xlim([0 14])
  ylim([-1.5 1.5])
  %grid on
  set(gca,'FontSize',18,'fontweight','bold');
  set(gcf,'Position',[100,100,800,600]);
  set(gcf,'color','w')
  xlabel('pH','fontSize',18)
  ylabel('E, V_{SHE}','fontSize',18)

  plot(pH_Zn_Zn2p, E_Zn_Zn2p,'k',  'lineWidth', 2);
  plot(pH_Zn_ZnOH2, E_Zn_ZnOH2,'k',  'lineWidth', 2);
  plot(pH_Zn_ZnO2_2m, E_Zn_ZnO2_2m, 'k',  'lineWidth', 2);
  plot(pH_Zn2p_ZnOH2,  E_Zn2p_ZnOH2, 'k',  'lineWidth', 2);
  plot(pH_ZnOH2_ZnO2_2m, E_ZnOH2_ZnO2_2m, 'k',  'lineWidth', 2);

 title('pH-potential diagram for 1 \mu{}M Zn ions in aqueous solution','fontSize',16)
 
 text(4.5, 0.3, 'Zn^{2+}', 'fontSize', 16)
 text(4.5, -1.2, 'Zn', 'fontSize', 16)
 text(8.5, 0.3, 'Zn(OH)_2', 'fontSize', 16)
 text(12.5, -1, 'ZnO_2^{2-}', 'fontSize', 16)
hold off
 
%% 1 October 2015
%  Plot Pourbaix pH-potential diagram for Iron (simplified)
% - not finished - fill in another way.

figure(12)
 pHdata       = [0, 14];
 E_Hydrogen   = [0, -0.059*14];
 E_Oxygen     = [1.22, 1.22-0.059*14];
 E_Fe2p_Fe3p  = [+0.771, +0.771];
 E_Fe_Fe2p    = [-0.62, -0.62];
 E_ = [0.441 - 0.1182*10.75 - 8*0.0295, 0.441 - 0.1182*14 - 8*0.0295];
 E_ = [-0.937, 1.5];
 E_=[-0.439 - 10.75*0.059, 1.5];
 
pH_Fe2p_Fe3p   = [0, 1.75];
pH_Fe_Fe2p     = [0, 9];
pH_ = [10.75, 14];
pH_ = [8.44, 8.44];
pH_ = [10.75, 10.75];

plot(pHdata, E_Hydrogen, 'k--', 'lineWidth', 1);
hold on
  plot(pHdata, E_Oxygen,   'b--', 'lineWidth', 1);
  
  xlim([0 14])
  ylim([-1.5 1.5])
  %grid on
  set(gca,'FontSize',18,'fontweight','bold');
  set(gcf,'Position',[100,100,800,600]);
  set(gcf,'color','w')
  xlabel('pH','fontSize',18)
  ylabel('E, V_{SHE}','fontSize',18)

  plot(pH_Fe2p_Fe3p, E_Fe2p_Fe3p,'k',  'lineWidth', 2);
  plot(pH_Fe_Fe2p,   E_Fe_Fe2p,'k',  'lineWidth', 2);
%   plot( ,  , 'k',  'lineWidth', 2);
%   plot( ,  , 'k',  'lineWidth', 2);
%   plot( ,  , 'k',  'lineWidth', 2);

 title('pH-potential diagram for 1 \mu{}M Fe ions in aqueous solution','fontSize',16)
 
%  text(4.5, 0.3, 'Zn^{2+}', 'fontSize', 16)
%  text(4.5, -1.2, 'Zn', 'fontSize', 16)
%  text(8.5, 0.3, 'Zn(OH)_2', 'fontSize', 16)
%  text(12.5, -1, 'ZnO_2^{2-}', 'fontSize', 16)
hold off
%% Figure 2 
% Plot the Bulter-Volmer equation with overpotential on x-axis
 
eta = -0.4:0.01:0.4; % volts of overpotential
i0  = 1;             % micro A/cm^2 h2 on platinum
RT  = 8.3*298;
azF = 0.5*1*96500;

I  = i0.*exp(azF.*eta./RT) - i0.*exp(-azF.*eta./RT); % Net current
If = i0.*exp(azF.*eta./RT) ;
Ib = - i0.*exp(-azF.*eta./RT);

figure(2)
 plot(eta, I, 'lineWidth', 4);
 
  xlabel('polarisation / Volts','fontSize',12)
  ylabel('current density / (\mu{}A / cm^2)','fontSize',12)

  grid on
  set(gca,'FontSize',12);
  %set(gca,'FontWeight','bold');
  set(2,'Position',[100,100,800,600]);
  set(gcf,'color','w')
  %    set(gca,'yTick',[-50:25:50])
      
  % xlim([-0.1 0.1])    
 %  ylim([])
 % ylim([-8 8])
      
  hold on
   plot(eta, If, 'r--','lineWidth', 2);
   plot(eta, Ib, 'g', 'lineWidth', 2);
  hold off
  
  xlim([-0.15 0.15]);
  ylim([-15 15]);
  
  legend('Net Current Density','Anodic Current Density',...
         'Cathodic Current Density', 'Location', 'NW')
  
%  Figure 2b (21)
% Demonstration of an Evans diagram with the data from figure 2
% Log current density on x-axis, potential on y-axis. 

LogAbsI  = log10(abs(I)); 
LogAbsIf = log10(abs(If));
LogAbsIb = log10(abs(Ib));

figure(21)
  plot(LogAbsI, eta, 'b', 'lineWidth', 2);
  grid on
  set(gca,'FontSize',18);
  set(21,'Position',[100,100,800,600]);
  set(gcf,'color','w')
  xlabel('log_{10}(absolute current density / (\mu{}A cm^{-2}) )','fontSize',18)
  ylabel('overpotential / Volts','fontSize',18)
  xlim([-0.5 2]);
  ylim([-0.2 0.2]);

  
%
%%  Diffusion limited oxygen cathode polarisation
  beta = 0.059 ;  % Tafel slope, 59 mV/decade
  RT   = 8.3*298;
  zF   = 4*96500; % Oxygen 4-electron process charge transport
  
  iL = 1E-2; % mA/cm2
  i0 = 1E-8;   % mA/cm2
  
  i = 1E-8:1E-4:1;
  
  % Overpotential based on surface concentration
  etaTotal  = -beta.*log10(i./i0) + 2.3.*(RT/zF).*log10(1 - i./iL);
  etaTotal2 = -beta.*log10(i./i0) + 2.3.*(RT/zF).*log10(1 - i./(iL*10));
  etaTotal3 = -beta.*log10(i./i0) + 2.3.*(RT/zF).*log10(1 - i./(iL*100));

  
  etaTotal(i >= iL) = -3
  etaTotal2(i >= iL*10) = -3
  etaTotal3(i >= iL*99) = -3
  % i = [i,100];
  
  figure(3)
  plot(log10(i), 1.22+etaTotal, 'lineWidth', 2)
  
  set(gca,'FontSize',18);
  set(3,'Position',[100,100,800,600]);
  set(gcf,'color','w')

  ylabel('E, V_{SHE}','fontSize',18)
  xlabel('log_{10}(i)','fontSize',18)
  xlim([-8 0.5])
  ylim([0.5 1.25])
  title('Oxygen reduction, i_{0}=10^{-8}, i_{L}=0.01 (A/cm^2)','fontSize',18)
hold on
  plot(log10(i), 1.22+etaTotal2,'r--', 'lineWidth', 2)
  
  plot(log10(i), 1.22+etaTotal3,'g-.', 'lineWidth', 2)
  % title([])
 legend('i_L', '10 i_L', '100 i_L')
hold off

% Mixed potential, oxygen transport:
title([])
hold on
d = 1.5; %vols
plot([-8, -8 + d/0.059 ], [-0.76, -0.76+d],  'k', 'lineWidth', 4)
hold off

ylim([-0.8, 1.25])
xlim([-8 6])
legend('i_L', '10 i_L', '100 i_L', '(M|M^{2+})')


%% 2015 Version: Cu--oxygen
% Second demonstration plot of oxygen diffusion limitation

E0ox = 1.22; % Volts_she
i0ox = 1E-5; % mA/cm^2
bTox = -0.59; % Volts per decade, cathode, negative. Guess 2-e process
E0cu = +0.337; % Volte_she
i0cu = 1E-5; % This is a guess, mA/cm^2
bTcu = 0.59; % 

% Add exact Butler-Volmer plot to the graph.
RT     = 298*8.31;
alfzFcu = 0.5*2*96500;
alfzFox = 0.5*2*96500; % Guess 2-electron process
zFox    = 4*96500;

logiData = -6:0.02:3;

iL_1 = 2;   % Guess, mA/cm^2.
iL_2 = 20;  % High iL
iL_3 = 500; % Effectively diffusion-unlimited

Ecu  = E0cu + (RT/alfzFcu)*asinh((10.^logiData)./(2*i0cu));
Eox  = E0ox - (RT/alfzFox)*asinh((10.^logiData)./(2*i0ox)) ...
       +2.3*(RT/zFox).*log10(1 - (10.^logiData)./iL_1);

Eox2 = E0ox - (RT/alfzFox)*asinh((10.^logiData)./(2*i0ox)) ...
       +2.3*(RT/zFox).*log10(1 - (10.^logiData)./iL_2);
Eox3 = E0ox - (RT/alfzFox)*asinh((10.^logiData)./(2*i0ox)) ...
       +2.3*(RT/zFox).*log10(1 - (10.^logiData)./iL_3);
   
Eox(10.^logiData > iL_1) = -0.5; % Effectively fix to minus off-scale
Eox2(10.^logiData > iL_2) = -0.5;
Eox3(10.^logiData > iL_3) = -0.5;

   
figure(23)
plot(logiData, Ecu, 'r--', 'lineWidth', 4)  % (Fe|Fe2+)
hold on
 plot(logiData, Eox, 'b-',  'lineWidth', 4') % (oxygen cathode)
hold off

xlim([-6 3])
ylim([0.2 1.4])

ylabel('E / Volts SHE','fontSize',16)
xlabel('log_{10}(i / (mA/cm^{2}))','fontSize',16)

legend('i_{a}(Cu|Cu^{2+})', 'i_c(O_2|H_{2}O), i_{L}=2mA/cm^2', ...
         'Location', 'NE')
set(gcf,'Position',[100,100,800,600]);
set(gcf,'color','w')
set(gca,'FontSize',14);
grid on

% Add extra lines for comparison
hold on
  plot(logiData, Eox2, 'k-',  'lineWidth', 2) % (oxygen cathode)
  plot(logiData, Eox3, 'g-.',  'lineWidth', 2') % (oxygen cathode)
hold off

legend('i_{a}(Cu|Cu^{2+})', 'i_c(O_2|H_{2}O), i_{L}=2mA/cm^2', ...
       'i_c(O_2|H_{2}O), i_{L}=20mA/cm^2', ...
       'i_c(O_2|H_{2}O), i_{L}=500mA/cm^2', ...
         'Location', 'NE', 'FontSize', 12)

%% 2015 Version
% Two Electrode Corrosion, no mass transport limitation (i.e. large i_lim)
% Mixed potential diagram, iron-hydrogen
%
% Data
%    E_0 hydr = 0.00  Volts SHE
%    E_0 iron = -0.44 Volts SHE
%    i_0 hydr = 1E-3  mA/cm^2
%    i_0 iron = 1E-5  mA/cm^2
%    Betahydr = 0.118 V/decade at 25 Celsius
%    Betairon = 0.059 V/decade at 25 Celsius

figure(4)

plot([-5, -5+0.44/0.059 ], [-0.44, 0], 'k', 'lineWidth', 4) % iron
hold on
 plot([-3, -3+0.44/0.118 ], [0, -0.44], 'b--', 'lineWidth', 4) % h2
hold off

ylim([-0.5, 0.1])
xlim([-6 4])

grid on
  set(gca,'FontSize',16);
  set(gcf,'Position',[100,100,1000,500]);
  set(gcf,'color','w')

  ylabel('E / Volts SHE','fontSize',18)
  xlabel('log_{10}(i / (mA/cm^{2}))','fontSize',18)
  legend('i_{a}(Fe|Fe^{2+}) Tafel', 'i_c(H_2| H^{+}) Tafel', ...
         'Location', 'SE')

% Add exact Butler-Volmer plot to the graph.
RT     = 298*8.31;
alfzFiron = 0.5*2*96500;
alfzFhydr = 0.5*1*96500;
i0iron = 1E-5;
E0iron = -0.44;
i0hydr = 1E-3;
E0hydr = 0;

logiData = -6:0.1:1;

Eiron  = E0iron + (RT/alfzFiron)*asinh((10.^logiData)./(2*i0iron));
Ehydr  = E0hydr - (RT/alfzFhydr)*asinh((10.^logiData)./(2*i0hydr));

hold on
 plot(logiData, Eiron, 'r--', 'lineWidth', 1)  % (Fe|Fe2+)
 plot(logiData, Ehydr, 'b-',  'lineWidth', 1') % (H+|H2)
hold off

legend('i_{a}(Fe|Fe^{2+}) Tafel', 'i_c(H^{+}|H_2) Tafel', ...
       'i_{a}(Fe|Fe^{2+}) exact', 'i_c(H^{+}|H_2) exact', ...
         'Location', 'SE')
grid on
     
% Add Tafel-equation E_0, i_0 points to graph
hold on
  scatter(-3, 0, 100, 'bo')     % (Fe|Fe2+) E_0, i_0
  scatter(-5, -0.44, 100, 'ko') % (H2|H+)   E_0, i_0 on Fe.
hold off


%% 2015 Version: Cu--oxygen with resistance
% Second demonstration plot of corrosion of copper by dissolved oxygen
% Assumptions:
% Activation limitation as per Tafel law (plotted as full Butler-Volmer)
% Let there be diffusion limitation, but with i_lim larger than the actual
% value of the corrosion current density
% And ohmic overpotential due to solution ion resistivity

E0ox = 1.22 - 0.05; % Volts_SHE   equilibrium potential for oxygen (Nernst correction for, say 1 mM conc)
i0ox = 1E-5; % mA/cm^2     exchange current density for oxygen
bTox = -0.59; % V, Tafel   slope in Volts per decade, cathode, negative. 
E0cu = +0.337; % Volts_she equilibrium potential for (Cu|Cu^2+), say at 1 Molar [Cu^2+] 
i0cu = 1E-4; % mA/cm^2     exchange current density for copper. This is a guess. 
bTcu = 0.59; % 

% Add exact Butler-Volmer plot to the graph.
RT     = 298*8.31;
alfzFcu = 0.5*2*96500;
alfzFox = 0.5*2*96500; % Guess 2-electron process
zFox    = 4*96500;

logiData = -6:0.02:3;

iL_1 = 500;   % mA/cm^2 arbitary diffusion-limited current density 
              %         Note this is a rather high limit for oxygen
              %         But could be possible in a fast flow

len  = 1E-2;  % 1 cm    Characteristic lengthscale for corrosion (guess)
% Ion resitivity for seawater, tap water, and "pure water":
% http://www.corrosion-doctors.org/Corrosion-Kinetics/Ohmic-drop-water.htm
rho1 = 0.3;   % Guess fresh water: 30 ohm cm =  0.3 ohm m
rho2 = 50;    % Guess tap water 5000 ohm cm = 50 ohm m
rho3 = 12000;  % Distilled: 5000 ohm m, "Pure" 1.2E6 ohm cm = 1.2E4 ohm m

EcuR1  = E0cu + (RT/alfzFcu)*asinh((10.^logiData)./(2*i0cu)) ...
       + 0.5*(0.001)*(10.^logiData)*len*rho1;
EoxR1= E0ox - (RT/alfzFox)*asinh((10.^logiData)./(2*i0ox)) ...
       + 2.3*(RT/zFox).*log10(1 - (10.^logiData)./iL_1) ...
       - 0.5*(0.001)*(10.^logiData)*len*rho1;

EcuR2  = E0cu + (RT/alfzFcu)*asinh((10.^logiData)./(2*i0cu)) ...
       + 0.5*(0.001)*(10.^logiData)*len*rho2;
EoxR2= E0ox - (RT/alfzFox)*asinh((10.^logiData)./(2*i0ox)) ...
       + 2.3*(RT/zFox).*log10(1 - (10.^logiData)./iL_1) ...
       - 0.5*(0.001)*(10.^logiData)*len*rho2;

EcuR3  = E0cu + (RT/alfzFcu)*asinh((10.^logiData)./(2*i0cu)) ...
       + 0.5*(0.001)*(10.^logiData)*len*rho3;
EoxR3= E0ox - (RT/alfzFox)*asinh((10.^logiData)./(2*i0ox)) ...
       + 2.3*(RT/zFox).*log10(1 - (10.^logiData)./iL_1) ...
       - 0.5*(0.001)*(10.^logiData)*len*rho3;   

EoxR1(10.^logiData > iL_1) = -0.5; % Effectively fix to minus off-scale
EoxR2(10.^logiData > iL_1) = -0.5;

   
figure(23)
plot(logiData, EcuR1, 'r--', 'lineWidth', 2)  % (  sea water
hold on
 plot(logiData, EcuR2, 'r-',  'lineWidth', 1') %  tap
 plot(logiData, EcuR3, 'r-',  'lineWidth', 3') %  pure
 plot(logiData, EoxR1, 'b--',  'lineWidth', 2') % (oxygen cathode) sea
 plot(logiData, EoxR2, 'b-',  'lineWidth', 1') % (oxygen cathode) tap
 plot(logiData, EoxR3, 'b-',  'lineWidth', 3') % (oxygen cathode) pure
hold off

xlim([-1 3])
ylim([0.5 1.2])

ylabel('E / Volts SHE','fontSize',16)
xlabel('log_{10}(i / (mA/cm^{2}))','fontSize',16)
set(gcf,'Position',[100,100,800,600]);
set(gcf,'color','w')
set(gca,'FontSize',14);
grid on
legend('i_{a}(Cu|Cu^{2+}), sea', ...
       'i_{a}(Cu|Cu^{2+}), tap', ...
       'i_{a}(Cu|Cu^{2+}), pure', ...
       'i_{c}(H_{2}O|O_2), sea', ...
       'i_{c}(H_{2}O|O_2), tap', ...
       'i_{c}(H_{2}O|O_2), pure', ...
         'Location', 'NE')


%% 2015 Evans diagram showing passivation
%
figure(5)
plot([-5,-5+0.46/0.059,-6,-6,4 ], [-0.46,+0.1, +0.1, +0.8,1.3], 'k', 'lineWidth', 4) % iron

hold on
  plot( [-8, 4] ,[0, -12*0.059], 'r--', 'lineWidth', 2);
  plot( [-8, 4] ,[0+0.85, 0.85+ -12*0.059], 'b-.', 'lineWidth', 2);
hold off

ylim([-0.5, 1.6])
xlim([-9 5])

grid on
  set(gca,'FontSize',18);
  set(gcf,'Position',[100,100,800,600]);
  set(gcf,'color','w')

  ylabel('E / V_{SHE}','fontSize',18)
  xlabel('log_{10}(i / (mA/cm^{2}))','fontSize',18)

  legend('i_{a}(Fe|Fe^{2+})', ...
         'i_{c}(weak oxidation)', ...
         'i_{c}(strong oxidation)', ...
          'Location', 'NW')


%% Tafel plots for cathodic inhibitor
figure(7) % Bimetallic diagram
plot([-10, -5,-5+0.46/0.059],[-0.46,-0.46,+0, ], 'k', 'lineWidth', 3) % zinc
 
grid on
  set(gca,'FontSize',18);
  set(gcf,'Position',[100,100,800,600]);
  set(gcf,'color','w')

  ylabel('E, V_{SHE}','fontSize',18)
  xlabel('log_{10}(i / (mA/cm^{2}))','fontSize',18)

  xlim([-10, 4])
  ylim([-0.5 0.1])
  
hold on
  plot([-10,-6, -6+0.5/0.059], [-0,-0,-0.5 ], 'b--.', 'lineWidth', 4) % zinc
  plot([-9, -9+0.5/0.059], [-0,-0.5 ], 'r-.', 'lineWidth', 4) % zinc
hold off
  legend('Iron anode (Fe|Fe^{2+})', '(H_{2}|H^{+}) on iron',...
         '(H_{2}|H^{+}) inhibited')
  
%% Tafel plots for 2-metal corrosion
figure(8)
plot([-9, -5,-5+0.46/0.059], [-0.46,-0.46,+0, ], 'k-', 'lineWidth', 4) % zinc

 
grid on
  set(gca,'FontSize',12);
  set(gcf,'Position',[100,100,500,400]);
  set(gcf,'color','w')
  grid off

  ylabel('E, V_{SHE}','fontSize',16)
  xlabel('log_{10}(i)','fontSize',16)

  xlim([-6.5, 4])
  ylim([-0.5 0.2])
  
hold on
  plot([-9,-6, -6+0.5/0.059], [-0,-0,-0.5 ], 'b-.', 'lineWidth', 2) % zinc
  % plot([-0,-0.5 ],[-9, -9+0.5/0.059], 'm-.', 'lineWidth', 2) % zinc
hold off
  legend('Iron anode i(Fe|Fe^{2+})', 'i(H_{2}|H^{+})')
  
  
%%  SCC
% Create log-log graph of flaw size vs failure stress 
% Add text labels later
figure(9)
 plot(0,0)
 grid on

  set(gca,'FontSize',18);
  set(gcf,'Position',[100,100,800,600]);
  set(gcf,'color','w')

  xlim([0 5])
  ylim([0 5])
  
  xlabel('(log) flaw size ','fontSize',18)
  ylabel('(log) stress','fontSize',18)


%%  Study bimetallic corrosion computationally:
Area_zinc = 0.1;         % Square metres
Area_iron = 1;

E_0_iron     = -0.44;    % Volts SHE (using standard electrode potential)
E_0_zinc     = -0.76;
E_0_hydrogen =  0.00;

i0_ironAnode         = 1E-7; % Exchange current densities (say in A/m^2)
i0_zincAnode         = 1E-7;
i0_H2_cathode_onIron = 1E-6;
i0_H2_cathode_onZinc = 1E-10;
  
potential = -0.76:0.001:0.1;    % Possible mixed surface potentials

iron_overpotential = potential - E_0_iron;
zinc_overpotential = potential - E_0_zinc;
hydr_overpotential = potential - E_0_hydrogen;

% Evaluate currents, here using simplified Butler-Volmer equation
F  = 96500;
twoRT = 8.31*298*2;

% Metals
i_iron_anode = i0_ironAnode * sinh( 2*96500*iron_overpotential/twoRT );
i_zinc_anode = i0_zincAnode * sinh( 2 * F * zinc_overpotential / twoRT );
% Protons (water)
i_H2_cath_on_iron =i0_H2_cathode_onIron*sinh(1*F*hydr_overpotential/twoRT);
i_H2_cath_on_zinc =i0_H2_cathode_onZinc*sinh(1*F*hydr_overpotential/twoRT);

total_current = i_iron_anode * Area_iron      + ...
                i_zinc_anode * Area_zinc      + ...
                i_H2_cath_on_iron * Area_iron + ...
                i_H2_cath_on_zinc * Area_zinc;

figure(8)

plot(potential, log10(abs(total_current)), 'lineWidth', 3 )
 grid on
  set(gca,'FontSize',16);
  set(gcf,'Position',[100,100,800,600]);
  set(gcf,'color','w')

  xlim([-0.8 0.01])
  ylim([-10 2])
  
  xlabel('Mixed Potential, V_{SHE} ','fontSize',24)
  ylabel('Log_{10} (absolute current / A )','fontSize',24)
  
  hold on
   plot(potential, log10(abs(i_zinc_anode * Area_zinc)), 'r','lineWidth', 1 )
   plot(potential, log10(abs(i_iron_anode * Area_iron)), 'g','lineWidth', 1 )
  hold off

figure(81)
plot(potential, abs((total_current).^(1/9)).*(sign(total_current)), 'lineWidth', 3 );
 grid on
  set(gca,'FontSize',18);
  set(gcf,'Position',[100,100,800,600]);
  set(gcf,'color','w')

hold on 
  plot(potential, abs((i_zinc_anode*Area_zinc).^(1/9)).*sign(i_zinc_anode), 'r','lineWidth', 2 )
  plot(potential, abs((i_iron_anode*Area_zinc).^(1/9)).*sign(i_iron_anode), 'g','lineWidth', 2 )
  plot(potential, abs((i_H2_cath_on_iron).^(1/9)).*sign(i_H2_cath_on_iron), 'k','lineWidth', 2 )
  plot(potential, abs((i_H2_cath_on_zinc).^(1/9)).*sign(i_H2_cath_on_zinc), 'k--','lineWidth', 2 )
hold off  

xlim([-0.8 0.01])
ylim([-1 1.5])
  
xlabel('Mixed Potential, V_{SHE} ','fontSize',18)
ylabel('Current^{1/9} / A^{1/9} ','fontSize',18)

legend('Total current', 'Zn|Zn^{2+}', 'Fe|Fe^{2+}', 'H_2|H^{+} on iron', ...
       'H_2|H^{+} on zinc',...
    'Location', 'NW')

%% Plot Evans diagram to illustrate a dominate case of bimetallic corrosion
%
figure(31)

plot([-10, -4, -4 + 0.76/0.1], [-0.76, -0.76, 0],'k', 'lineWidth', 6 )

grid on
xlim([-8 3])
ylim([-0.8 0.4])
ylabel('E, V_{SHE}','fontSize',16)
xlabel('log_{10}(i)','fontSize',16)

hold on
  plot([-10, -3, -3 + 0.44/0.1], [-0.44, -0.44, 0], 'g--', 'lineWidth', 6 );
  plot([-10, -7, -7 + 0.76/0.1], [0, 0, -0.76], 'b-.', 'lineWidth', 4)
  plot([-10, -3, -3 + 0.76/0.1], [0, 0, -0.76], 'm--', 'lineWidth', 4)
hold off

set(gca,'FontSize',14);
set(gcf,'Position',[100,100,800,600]);
legend('i(Zn|Zn^{2+})','i(Fe|Fe^{2+})','i(H^{+}|H_{2}) on Zn',...
         'i(H^{+}|H_{2}) on Fe', 'Location', 'NE')
set(gcf,'color','w')

%%
% High temp corrosion examples

t = 0:0.01:1.5;

mLoss = -0.1*t;
mLin  = +0.5*t;
mPar  = +1*t.^0.5;
mLog  = +0.1*log(t/0.01);

figure(40)
plot(t, mLoss, 'b', 'lineWidth', 2);
hold on
  plot(t, mLin, 'r', 'lineWidth', 1);
  plot(t, mPar, 'k--', 'lineWidth', 2);
  plot(t, mLog, 'k.-', 'lineWidth', 1);
hold off

ylabel('Mass Change / units','fontSize',16)
xlabel('Time / units','fontSize',16)

legend('Mass loss','Linear','Parabolic','Logarithmic','Location', 'NW') 

set(gca,'FontSize',14);


set(gcf,'color','w')

%%
% Hi-T corrosion examples sheet

tData = [5,10,30,60,120,240,360];

m700d = [4.4, 6.4,11.2,15.7,22.2,31.9,38.1];
m900d = [24.0,33.9,57.5,83.2,116.2,165.3,200.1];

figure(32)
plot(tData, m700d, 'b-', 'lineWidth', 2);
hold on 
  plot(tData, m900d, 'r--', 'lineWidth', 2);
hold off
xlabel('time/ min','fontSize',14)
ylabel('\Delta (mass/area) /\mu{}g/mm^2 ','fontSize',14)
set(gca,'FontSize',14);
legend('700 C','900C', 'Location', 'NW')
set(gcf,'color','w')
grid on

figure(33)
plot(log10(tData), log10(m700d), 'b-', 'lineWidth', 2);
hold on 
  plot(log10(tData), log10(m900d), 'r--', 'lineWidth', 2);
hold off
xlabel('log (time/ min)','fontSize',14)
ylabel('log (\Delta (mass/area) /\mu{}g/mm^2) ','fontSize',14)
set(gca,'FontSize',14);
legend('700 C','900C', 'Location', 'NW')
set(gcf,'color','w')
grid on
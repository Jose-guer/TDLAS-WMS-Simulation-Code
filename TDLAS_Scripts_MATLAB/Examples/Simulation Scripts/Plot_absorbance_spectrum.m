%% %%%%%%%%%%%%%%%%%%%% Example 1 %%%%%%%%%%%%%%%%%%%%%%%
% This example simulates the absorption spectrum at a specified T, P, Xi
% where the function plots the absorbance spectrum.
%
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

clc; 
clear;
close all;

species = 'H2O'; % molecule
T = [1200,1800,2200]; % Temperature [K]
P = 1; % Pressure [atm]
Xrad = 0.1; %Radiating species mole fraction
L = 6; % Path length [cm]
dwn = 0.001; % wavenumber spacing [cm^-1] 

% wnL = 7185;  %[cm^-1] - min wavenumber
% wnH = 7186;  %[cm^-1] - max wavenumber

wnL = 6805.5;  %[cm^-1] - min wavenumber
wnH = 6806.5;  %[cm^-1] - max wavenumber

wn = wnL:dwn:wnH;

%This function is meant as a tool for quickly plotting the absorbance
%spectrum. For more versatility see the next example.
SpectrumPlot(T,P,Xrad,L,wn,species) 

%% %%%%%%%%%%%%%%%%%%%% Example 2 %%%%%%%%%%%%%%%%%%%%%%%
% This example simulates the absorption spectrum at a specified T, P, Xi
% where the function returns the wavenumber vector and absorbance spectrum.

clc;

species = 'H2O'; %absorbing molecule, H2O is default
T = 1200; % Temperature [K]
P = 1; % Pressure [atm]
Xrad = 0.1; %Radiating species mole fraction
L = 6; % Path length [cm]

dwn = 0.001; % wavenumber spacing [cm^-1] 
wnL = 7185.2;    %[cm^-1] - min wavenumber
wnH = 7186;  %[cm^-1] - max wavenumber

% wnL = 6805.5;  %[cm^-1] - min wavenumber
% wnH = 6806.5;  %[cm^-1] - max wavenumber

[X, absorption] = HITRAN_Spectrum_Plot(T,P,Xrad,L,dwn,wnL,wnH,species); 

figure; set(gcf,'Position',[403 289 732 595])
plot(X,absorption)
box off
xlabel('wavenumber [$\mathrm{cm^{-1}}$]')
ylabel('$\alpha(\nu)$')
xlim([wnL,wnH])


%% %%%%%%%%%%%%%%%%%%%% Example 3 %%%%%%%%%%%%%%%%%%%%%%%
% This example simulates the absorption spectrum in two ways. The first
% uses the specified T, P, Xi and the second using T, Pi, and wLi. 

species = 'H2O'; % molecule
T = 1800; % Temperature [K]
P = 1; % Pressure [atm]
Xrad = 0.1; %Radiating species mole fraction
Pi = P*Xrad;
L = 6; % Path length [cm]
dwn = 0.001; % wavenumber spacing [cm^-1] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wnL = 7185;  %[cm^-1] - min wavenumber
% wnH = 7186;  %[cm^-1] - max wavenumber
% 
% gamma_self = 0.198; %7185.596
% n_self = 0.53;
% gamma_air = 0.0403;
% n_air = 0.587;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wnL = 6805.5;  %[cm^-1] - min wavenumber
wnH = 6806.5;  %[cm^-1] - max wavenumber

gamma_self = 0.205; %6806.03
n_self = 0.86;
gamma_air = 0.0104;
n_air = -0.164;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wn = wnL:dwn:wnH;

%try to estimate the collisional width to enable comparison
wLi = P * (Xrad*gamma_self * (296/T)^n_self + (1-Xrad)*gamma_air*(296/T)^n_air);

[X1, absor1] = HITRAN_Spectrum_Plot_Alt(T,Pi,wLi,L,dwn,wnL,wnH,species);
[X2, absor2] = HITRAN_Spectrum_Plot(T,P,Xrad,L,dwn,wnL,wnH,species); 


% The differences in the peak amplitude near lineecnter are due to using 
% different self broadening parameters to estimate the doublet-lines 
% collisional width. When specifiying the collisional width, 
% pressure-line-shifting is also not modeled. To improve the comparison,
% pressure line shifting can be turned off by commenting that line of
% code in HITRAN_Spectrum_Plot.m.

figure; set(gcf,'Position',[403 289 732 595]);
hold on
plot(X1,absor1)
plot(X2,absor2)
hold off
box off
legend('$\alpha(\nu,T,P_i,\Delta\nu_c)$','$\alpha(\nu,T,P,X_i)$','box','on')
xlabel('wavenumber [$\mathrm{cm^{-1}}$]')
ylabel('$\alpha(\nu)$')
xlim([wnL,wnH])

%% %%%%%%%%%%%%%%%%%%%% Example 4 %%%%%%%%%%%%%%%%%%%%%%%
% This simulate the Cluster of CO lines near 5 microns commonly used in
% thermometry in the mid-IR.

clc; clear; 
close all;

T = 2000; %[K]   - temperature
L = 2; %[cm]  - path length
dwn = 0.001; %[cm^-1] - wavenumber spacing
XCO = 0.05;

threshold = 10^-4;
wn_units = 1;
S_units = 1;

wnL = 2008;
wnH = 2009;

species = 'CO_HiTemp'; 

figure
hold on

P = 2; %atm
[X_CO, abs_CO] = HITRAN_Spectrum_Plot(T,P,XCO,L,dwn,wnL,wnH,species);
plot(X_CO,abs_CO,'LineWidth',2)

P = 5; %atm
[X_CO, abs_CO] = HITRAN_Spectrum_Plot(T,P,XCO,L,dwn,wnL,wnH,species);
plot(X_CO,abs_CO,'LineWidth',2)

XCO = 1;
[vo_CO,Es_CO,XS_CO,So_CO] = getHITRAN(T,XCO,wnL,wnH,threshold,species,wn_units,S_units);

stem(vo_CO,XS_CO * 1,'k-','LineWidth',2,'Marker','none')
hold off
box off
xlabel('$\nu$ ($\mathrm{cm^{-1}}$)')
ylabel('$\alpha(\nu)$')
legend('P = 2 atm','P = 5 atm','box','off','location','northwest');
ylim([0,0.8])
set(gcf,'Position',[463 266 700 569])
xlim([wnL,wnH])

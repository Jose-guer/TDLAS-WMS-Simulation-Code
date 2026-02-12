%This script is used to simulate synthetic absorption spectra with noise.
%
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

clc; clear; close all;

species = 'H2O'; %absorbing molecule, H2O is default
T = 1800; % Temperature [K]
P = 1; % Pressure [atm]
Xrad = 0.15; %Radiating species mole fraction
L = 10; % Path length [cm]
dwn = 0.001; % wavenumber spacing [cm^-1] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wnL = 7185.3;    %[cm^-1] - min wavenumber
wnH = 7185.8;  %[cm^-1] - max wavenumber
x1 = linspace(wnL,wnH,500);

[X1, absor1] = HITRAN_Spectrum_Plot(T,P,Xrad,L,dwn,wnL,wnH,species); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wnL = 6806.03-0.3;  %[cm^-1] - min wavenumber
wnH = 6806.03+0.2;  %[cm^-1] - max wavenumber
x2 = linspace(wnL,wnH,500);

[X2, absor2] = HITRAN_Spectrum_Plot(T,P,Xrad,L,dwn,wnL,wnH,species); 

%%

y1 = interp1(X1,absor1,x1);
y2 = interp1(X2,absor2,x2);

A = 0.005; % noise
% A = 0; %noise free
noise1 = A*(rand(size(y1))-0.5);
noise2 = A*(rand(size(y2))-0.5);

y1 = y1 + noise1;
y2 = y2 + noise2;

dat.x1 = x1;
dat.absor1 = y1;
dat.x2 = x2;
dat.absor2 = y2;
dat.P = P;
dat.T = T;
dat.X = Xrad;
dat.L = L;

figure; set(gcf,'Position',[403 289 732 595])
plot(x1,y1)
box off
xlabel('wavenumber [$\mathrm{cm^{-1}}$]')
ylabel('$\alpha(\nu)$')


figure; set(gcf,'Position',[403 289 732 595])
plot(x2,y2)
box off
xlabel('wavenumber [$\mathrm{cm^{-1}}$]')
ylabel('$\alpha(\nu)$')


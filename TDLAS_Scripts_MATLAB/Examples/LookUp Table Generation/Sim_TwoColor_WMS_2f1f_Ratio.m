% This script simulates the two-color WMS-2f/1f ratio and WMS-2f/1f signals 
% at a specified pressure and laser specific tuning parameters
% 
% Hardcoded values are taken from [1] and Fig. 5 is reproduced.
% 
% [1] Schwartz, Charles J., et al. "Near-MHz temperature and H2O measurements 
% in post-detonation fireballs of 25 g hemispherical explosives using 
% scanned-wavelength-modulation spectroscopy." 
% Applied Optics 62.6 (2023): 1598-1609.
% 
% Linecenter pressure-shift can be commented out in
% "HITRAN_Spectrum_Plot.m" line 128, when generating scanned-WMS tables. 
% 
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

clc; clear;
close all;

%% ************SIMULATION INPUTS *******************

M = 10^6;
G = 10^9;

fs = 3*G;  %Sampling frequency

%Test gas details
species = 'H2O';

P = 1; % [atm]
T = 800:200:2400;					
H2O = [0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31]';	
L = 12.6; % [cm]

%Background details
Tbg = 296;
H2O_bg = 0;
Lbg = 0; % [cm]

%% Laser 1 - Absorption and WMS-2f/1f Signal
%
wn0 = 7185.596;
fmod = 35.0*M;   %Hz
a = 0.113;       %cm^-1
i0 = 0.160;
i2 = 3.18*10^-3;
psi1 = 6.08;
psi2 = 6.92;

t = -0.5/fmod:1/fs:0.5/fmod;
w = 2*pi*fmod;
wn_sim = wn0 + a*cos(w*t);

delta_wn = 0.001;
wnL = min(wn_sim);
wnH = max(wn_sim);

%%% Background absorption
[X, bg_abs] = HITRAN_Spectrum_Plot(Tbg,P,H2O_bg,Lbg,delta_wn,wnL,wnH,species);
tau_bg_copy = exp(-bg_abs);
tau_bg = interp1(X,tau_bg_copy,wn_sim);

%%% Background WMS signal
[X1f_bg,Y1f_bg,X2f_bg,Y2f_bg,X4f_bg,Y4f_bg] = wms_calibration_free(t,tau_bg,i0,psi1,i2,psi2);

S2f1f_L1 = zeros(length(H2O),length(T));

for i = 1:length(H2O)
    for j = 1:length(T)

        %combustion absorption
        [~, absorption] = HITRAN_Spectrum_Plot(T(j),P,H2O(i),L,delta_wn,wnL,wnH,species);
        tau = exp(-(absorption + bg_abs));
        tau_sim = interp1(X,tau,wn_sim);

        %Absorption Signal
        [X1f,Y1f,X2f,Y2f,X4f,Y4f] = wms_calibration_free(t,tau_sim,i0,psi1,i2,psi2);

        R1f = sqrt(X1f.^2 + Y1f.^2);
        R1f_bg = sqrt(X1f_bg.^2 + Y1f_bg.^2);

        S2f1f_L1(i,j) = sqrt( ((X2f./R1f) - (X2f_bg./R1f_bg)).^2 +  ((Y2f./R1f) - (Y2f_bg./R1f_bg)).^2 );

    end
end


%}

%% Laser 2 - Absorption and WMS-2f/1f Signal

wn0 = 6806.03;
fmod = 45.5*M; %Hz

a = 0.166; %cm^-1
i0 = 0.175;
i2 = 8.394*10^-3;
psi1 = 5.88;
psi2 = 4.22;

t = -0.5/fmod:1/fs:0.5/fmod;
w = 2*pi*fmod;
wn_sim = wn0 + a*cos(w*t);

delta_wn = 0.001;
wnL = min(wn_sim);
wnH = max(wn_sim);

%%% Background absorption
[X, bg_abs] = HITRAN_Spectrum_Plot(Tbg,P,H2O_bg,Lbg,delta_wn,wnL,wnH,species);
tau_bg_copy = exp(-bg_abs);
tau_bg = interp1(X,tau_bg_copy,wn_sim);

%%% Background WMS signal
[X1f_bg,Y1f_bg,X2f_bg,Y2f_bg,X4f_bg,Y4f_bg] = wms_calibration_free(t,tau_bg,i0,psi1,i2,psi2);

S2f1f_L2 = zeros(length(H2O),length(T));

for i = 1:length(H2O)
    for j = 1:length(T)

        %combustion absorption
        [~, absorption] = HITRAN_Spectrum_Plot(T(j),P,H2O(i),L,delta_wn,wnL,wnH,species);
        tau = exp(-(absorption + bg_abs));
        tau_sim = interp1(X,tau,wn_sim);

        %Absorption Signal
        [X1f,Y1f,X2f,Y2f,X4f,Y4f] = wms_calibration_free(t,tau_sim,i0,psi1,i2,psi2);

        R1f = sqrt(X1f.^2 + Y1f.^2);
        R1f_bg = sqrt(X1f_bg.^2 + Y1f_bg.^2);

        S2f1f_L2(i,j) = sqrt( ((X2f./R1f) - (X2f_bg./R1f_bg)).^2 +  ((Y2f./R1f) - (Y2f_bg./R1f_bg)).^2 );

    end
end


%% Plots

legStr = {};
for i = 1:length(H2O)
    legStr(i) = {['$X_{\mathrm{H_2O}}$ = ',num2str(H2O(i))]};
end

figure
plot(T,S2f1f_L2./S2f1f_L1)
xlabel('T [K]')
ylabel('Two-Color WMS-2f/1f Ratio')
box off
legend(legStr,'Location','NorthWest','box','on')
xlim([min(T),max(T)])
set(gcf,'position',[575 277 765 628])

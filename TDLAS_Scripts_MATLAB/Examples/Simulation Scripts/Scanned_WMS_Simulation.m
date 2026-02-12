% This script demonstrates how to simulate scanned-WMS signals using the
% calibration-free WMS model.
% 
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

clc; clear;
close all;

%constants and conversion factors
k = 10^3;
M = 10^6;
G = 10^9;


%% Simulation Inputs

%Radiating Species Details
species = 'H2O';
wn0 = 7185.596; %line center frequency [cm^-1]
% wn0 = 6806.03; %line center frequency [cm^-1]

%Test gas
T = 1800; %[K]
P = 1; %[atm]
Xrad = 0.15;
L = 10; %[cm]

%Background
T_bg = 293;
P_bg = 1;
Xrad_bg = 0.001;
L_bg = 10; %[cm]

%Laser characterization parameters

fm = 35.0*M;   %Hz, modulation frequency
fscan = 500*k; %Hz, scanning frequency
fs = 3*G;      %Hz, sampling frequency

% - modulation -
a = 0.116;    %cm^-1
i0 = 0.1671;
i2 = 2.48E-03;
psi1 = 1.9356*pi;
psi2 = 4.4138*pi;
Idc = 1.2957;

% - scanning -
as = 0.06;
psi_s = 1.411*pi;
i_s = 0.1661;


%% Incident light intensity simulation

wm = 2*pi*fm;
ws = 2*pi*fscan;
t = 0:1/fs:1/fscan;

wn_scan = wn0 + as*cos(ws*t);
I0_scan = 1 + i_s*cos(ws*t + psi_s);

wn_sim = wn_scan + a*cos(wm*t); %cm^-1
I0_sim = Idc*( 1 + i_s*cos(ws*t + psi_s) + i0*cos(wm*t + psi1) + i2*cos(2*wm*t + psi2) );

i0 = i0./I0_scan;
i2 = i2./I0_scan;

%% WMS signal simulation

delta_wn = 0.0005;
wnL = min(wn_sim);
wnH = max(wn_sim);

%background absorption
[X, bg_abs] = HITRAN_Spectrum_Plot(T_bg,P_bg,Xrad_bg,L_bg,delta_wn,wnL,wnH,species);
bg_abs_scan = interp1(X,bg_abs,wn_scan);
bg_abs_sim = interp1(X,bg_abs,wn_sim);

tau_bg_copy = exp(-bg_abs);

%combustion absorption
[X, absorption] = HITRAN_Spectrum_Plot(T,P,Xrad,L,delta_wn,wnL,wnH,species);
absor_scan = interp1(X,absorption,wn_scan);
absor_sim = interp1(X,absorption,wn_sim);
tau = exp(-(absorption+bg_abs));

%Transmitted light intensity signal
It_sim = I0_sim .* exp(-(absor_sim + bg_abs_sim));

S_2f1f = zeros(size(wn_scan));
S_2f = zeros(size(wn_scan));
S_1f = zeros(size(wn_scan));
S_4f = zeros(size(wn_scan));
S_4f2f = zeros(size(wn_scan));
tsim = -0.5/fm:1/fs:0.5/fm;

for j = 1:length(wn_scan)

    wn0 = wn_scan(j);
    wn_sim_new = wn0 + a*cos(wm*tsim);
    tau_sim = interp1(X,tau,wn_sim_new);
    tau_bg = interp1(X,tau_bg_copy,wn_sim_new);

    %Absorption Signal
    [X1f,Y1f,X2f,Y2f,X4f,Y4f] = wms_calibration_free(tsim,tau_sim,i0(j),psi1,i2(j),psi2);

    %%% Background signal
    [X1f_bg,Y1f_bg,X2f_bg,Y2f_bg,X4f_bg,Y4f_bg] = wms_calibration_free(tsim,tau_bg,i0(j),psi1,i2(j),psi2);

    %WMS Harmonics
    R1f = sqrt(X1f.^2 + Y1f.^2);
    R1f_bg = sqrt(X1f_bg.^2 + Y1f_bg.^2);

    S_2f1f(j) = sqrt( ((X2f./R1f) - (X2f_bg./R1f_bg)).^2 +  ((Y2f./R1f) - (Y2f_bg./R1f_bg)).^2 );
    S_2f(j) = sqrt((X2f-X2f_bg).^2 + (Y2f-Y2f_bg).^2);
    S_1f(j) = sqrt((X1f-X1f_bg).^2 + (Y1f-Y1f_bg).^2);
    S_4f(j) = sqrt((X4f-X4f_bg).^2 + (Y4f-Y4f_bg).^2);
    S_4f2f(j) = S_4f(j)/S_2f(j);

end


%% Plots

figure; set(gcf,'Position', [294 176 1174 755])

hAxis(1) = subplot(4,1,1);
plot(repmat(It_sim,1,5))
box off
legend('$I_t(t)$','box','on')

hAxs(2) = subplot(4,1,2);
plot(repmat(absor_scan,1,5))
box off
legend('$\alpha(\nu)$','box','on')

hAxis(3) = subplot(4,1,3);
hold on
plot(repmat(S_2f1f,1,5))
plot(repmat(S_4f * 3,1,5))
box off
legend('WMS-2$f$/1$f$','WMS-4$f$ $\times 3$','box','on')

hAxis(4) = subplot(4,1,4);
plot(repmat(i0,1,5))
box off
legend('$i_0(t)$','box','on')

linkaxes(hAxis,'x')


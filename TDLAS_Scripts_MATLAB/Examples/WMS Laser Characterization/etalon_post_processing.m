% This script loads in the etalon signal and labelled peaks and determines
% the modulation depth and the relative-phase shifts if the phase of the
% 1st and 2nd order fits to the intensity signal are specified. 
%
% The etalon peaks here have already been found and labeled depending on
% whether the wavenumber is increasing or decreasing. Each peak is evenly
% spaced by the FSR of the etalon so the peaks only need to be "counted"
% adding 1 when the wavenumber increases, subtracting 1 when the wavenumber
% decreases, or adding 0 after a frequency reversal event. The data here 
% illistrates the general process.
%
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

clc; clear; close all;

load('etalon_dat.mat') % loads in variable "dat"
rgb = [0,0.5,1];

%From intensity post-processing
phi1 = 0.9715*pi;
phi2 = 0.8672*pi;

%etalon detector signal
t = dat.t;
etalon = dat.etalon;
f = dat.f; %modulation frequency

%labelled etalon peaks
tpeaks = dat.tpeaks;
ypeaks = dat.ypeaks;
peak_label = dat.peak_label;
fsr = dat.fsr;

rel_freq = peak_label * fsr;

a = (max(rel_freq) - min(rel_freq))/2;
phi = pi;
b = mean(rel_freq); 
X0 = [a,phi,b];

fit_param = lsqnonlin(@(X) cosinefit_lsq(tpeaks,rel_freq,X(1),X(2),X(3),f),X0,[0,0,-1],[0.5,2*pi,1]);

a = fit_param(1);
a = round(a,3);
phi = fit_param(2);
b = fit_param(3);

%Relative phase-shift

psi_1 = phi1 - phi;
psi_1 = psi_1 / pi;

psi_2 = phi2 - 2*phi;
psi_2 = psi_2 / pi;

psi_1 = mod(psi_1 + 10,2);
psi_2 = mod(psi_2 + 10,2);

%% Plots

baseline = dat.baseline;
baseline = baseline - mean(baseline);
baseline = baseline./max(baseline);
m = 10;

figure; set(gcf,'Position',[307 211 1093 694])
hold on
plot(t(1:m:end),baseline(1:m:end),'k','LineWidth',2)
plot(t,cos(2*pi*f*t+phi),'color',rgb,'LineWidth',2)
title('FM/IM Phase shift')
xlabel('time (s)')
box off
legend('Intensity','Wavenumber','box','on')


figure; set(gcf,'Position',[307 211 1093 694])
hAxis(1) = subplot(3,1,1);
hold on
plot(t,etalon,'color',rgb)
plot(tpeaks,ypeaks,'k*')
hold off
box off
ylabel('Amplitude [V]')
legend('Etalon','Peaks','box','on')

hAxis(2) = subplot(3,1,2);
plot(tpeaks,peak_label,'k*')
grid on; grid minor;
box off
ylabel('Label')
legend('Peaks','box','on')
ylim([-10,10])

hAxis(3) = subplot(3,1,3);
hold on;
plot(tpeaks,rel_freq,'k*');
t_fit = linspace(tpeaks(1),tpeaks(end),2000);
plot(t_fit,a*cos(2*pi*f*t_fit+phi)+b,'r-','LineWidth',2)
xlabel('time (s)')
ylabel('rel. freq. $[\mathrm{cm^{-1}}]$')
box off
legend('Peaks','Best Fit','box','on')
ylim([-0.2,0.15])

linkaxes(hAxis,'x');

clc;

disp(['modulation_depth = ', num2str(a),' cm^-1'])
disp(['psi1 = ', num2str(psi_1),'*pi'])
disp(['psi2 = ', num2str(psi_2),'*pi'])

%%

function r = cosinefit_lsq(t,dat,a,phi,b,f)
r = dat - a*cos(2*pi*f*t + phi)- b;
end

% This script loads in the sample data and determines the intensity
% characterization parameters for WMS. First a cosine wave at 1f is fit to
% the baseline and then a cosine wave at 2f is fit to the residual. The DC
% normalized intensity modulation amplitudes are printed at the end.
% When the nonlinear solver fails to converge, generally chaning the inital
% phase guess solves the problem. The phase of the 1st and 2nd order fits
% are inputs into the etalon_post_processing.m script where the actual
% relative phase shifts will be determined.
% 
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu


clc;
clear;
close all;
rgb = [0,0.5,1];

%% Load Data

load('etalon_dat.mat') % loads in variable "dat"
t = dat.t;
baseline = dat.baseline;
etalon = dat.etalon;
f = dat.f;

figure
hold on;
plot(t,baseline,'k')
plot(t,etalon,'color',rgb)
legend('Baseline','Etalon','box','on')
hold off
box off
xlabel('time [s]')
ylabel('Amplitude [V]')

%% Intensity Characterization

%First order terms 

i0 = (max(baseline) - min(baseline))/2;
phi = 1.4*pi;
b = mean(baseline); 
X0 = [i0,phi,b];

fit_param = lsqnonlin(@(X) cosinefit_lsq(t,baseline,X(1),X(2),X(3),f),X0,[0.8*i0,0,0.8*b],[1.2*i0,2*pi,1.2*b]);

i0 = fit_param(1);
phi = fit_param(2);
b = fit_param(3);

figure; plot(t,baseline);
hold on;
plot(t,i0*cos(2*pi*f*t+phi)+b,'r--','LineWidth',2)
xlabel('time (s)')
ylabel('Amplitude [V]')
legend('Baseline','1st-order fit','box','on')

%Second order terms 

%residual
res = baseline - i0*cos(2*pi*f*t+phi) - b;

i2 = (max(res) - min(res))/2;
phi2 = 0.2*pi; 
b2 = 0;
X0 = [i2,phi2];
fit_param = lsqnonlin(@(X) cosinefit_lsq(t,res,X(1),X(2),b2,2*f),X0,[0.25*i2,0],[1.5*i2,2*pi]);

i2 = fit_param(1);
phi2 = fit_param(2);

figure; plot(t,res);
hold on;
plot(t,i2*cos(2*pi*2*f*t+phi2),'r-','LineWidth',2)
xlabel('time (s)')
ylabel('Amplitude [V]')
legend('Residual','2nd-order fit','box','on')
box off


%% Plots 

Idc = b;
scan_fit =  Idc + i0*cos(2*pi*f*t+phi) + i2*cos(2*pi*(2*f)*t+phi2);

figure;
hAxis(1) = subplot(2,1,1);
plot(t,baseline)
hold on
plot(t,scan_fit,'--')
legend('Baseline','Fit','box','on')
ylabel('Amplitude [V]')
box off

hAxis(2) = subplot(2,1,2);
plot(t,(scan_fit-baseline)*100./baseline)
ylabel('$\%$')
legend('Residual','box','on')
xlabel('time (s)')
box off

clc

%The actual relative phase-shifts are determined after post-processig the
%etalon signal

format short g
A = [i0/Idc, phi/pi, i2/Idc, phi2/pi, Idc];
T = array2table(A,...
    'VariableNames',{'i0','phi1/pi','i2','phi2/pi','Idc'})

%%

function r = cosinefit_lsq(t,dat,a,phi,b,f)
r = dat - a*cos(2*pi*f*t + phi)- b;
end

% This script loads in the aligned background and absorption intensity
% signals and extracts the background subtracted WMS harmonics. A power
% spectral density plot highlighting the hamonics and side-bands is also
% computed.
% 
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu


clc; clear; 
close all;

k = 10^3;
M = 10^6;

load('raw_data.mat')
fscan = raw.fscan;
fm1 = raw.fm1;
fm2 = raw.fm2;

load('aligned_data.mat')
t = dat(:,1);
sig = dat(:,2);
bg_sig = dat(:,3);

fs = 1/(t(2)-t(1));
N = length(sig);

%Lock-in Amplifier Details
filterOrder = 5;
fc = 2.2*fscan; %filter cutt-off frequency


%% Power Spectral Density of Raw Signal

figure
hold on
plot(t * 10^6,bg_sig)
plot(t * 10^6,sig)
box off
xlabel('time [$\mu$s]')
ylabel('Signal [V]')
legend('Background','Absorption Sig.','box','on')

% return

%******************** color PSD spectrum ***********************
%
PSD_signal = abs(fft(sig)).^2;
f_PSD = [0:N-1]/N*fs/10^6;

ind = f_PSD > 1 & f_PSD < 150;
PSD_signal = PSD_signal(ind);
f_PSD = f_PSD(ind);

%spectrum 1391.7 nm
ind = f_PSD > (fm1 - 2.2*fscan)/M & f_PSD < (fm1 + 2.2*fscan)/M;
PSD_1f = PSD_signal(ind);
freq_1f = f_PSD(ind);

ind = f_PSD > (2*fm1 - 2.2*fscan)/M & f_PSD < (2*fm1 + 2.2*fscan)/M;
PSD_2f = PSD_signal(ind);
freq_2f = f_PSD(ind);

%spectrum 1469.3 nm
ind = f_PSD > (fm2 - 2.2*fscan)/M & f_PSD < (fm2 + 2.2*fscan)/M;
PSD_1f_2 = PSD_signal(ind);
freq_1f_2 = f_PSD(ind);

ind = f_PSD > (2*fm2 - 2.2*fscan)/M & f_PSD < (2*fm2 + 2.2*fscan)/M;
PSD_2f_2 = PSD_signal(ind);
freq_2f_2 = f_PSD(ind);

%Plot
figure; 
semilogy(f_PSD,PSD_signal,'k');
hold on
semilogy(freq_1f,PSD_1f,'b');
semilogy(freq_2f,PSD_2f,'b');
semilogy(freq_1f_2,PSD_1f_2,'r');
semilogy(freq_2f_2,PSD_2f_2,'r');
xlabel('Frequency (MHz)')
ylabel('PSD ($\mathrm{V^2/Hz}$)')
xlim([15,150])
box off
xlim([30,93])
legend('','$\nu_o$ = 1391.7nm', '','$\nu_o$ = 1469.3nm','box','on')
%}
%****************************************************************


%%                          Extract WMS Harmonics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Laser 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,X1f,Y1f]  = WMS_harmonic(sig,t,fs,fc,fm1,1,filterOrder);
[~,X1f_bg,Y1f_bg]  = WMS_harmonic(bg_sig,t,fs,fc,fm1,1,filterOrder);

[~,X2f,Y2f]  = WMS_harmonic(sig,t,fs,fc,fm1,2,filterOrder);
[~,X4f,Y4f]  = WMS_harmonic(sig,t,fs,fc,fm1,4,filterOrder);

[~,X2f_bg,Y2f_bg]  = WMS_harmonic(bg_sig,t,fs,fc,fm1,2,filterOrder);
[~,X4f_bg,Y4f_bg]  = WMS_harmonic(bg_sig,t,fs,fc,fm1,4,filterOrder);

R1f = sqrt(X1f.^2 + Y1f.^2);
R1f_bg = sqrt(X1f_bg.^2 + Y1f_bg.^2);
WMS_2f1f_L1 = sqrt( ((X2f./R1f) - (X2f_bg./R1f_bg)).^2 +  ((Y2f./R1f) - (Y2f_bg./R1f_bg)).^2 );

WMS_2f_L1 = sqrt( (X2f - X2f_bg).^2 + (Y2f - Y2f_bg).^2 );
WMS_4f_L1 = sqrt( (X4f - X4f_bg).^2 + (Y4f - Y4f_bg).^2 );
WMS_4f2f_L1 = WMS_4f_L1./WMS_2f_L1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Laser 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The 4f harmonic (176 MHz) is outside the detectors bandwidth (150 MHz)

[~,X1f,Y1f]  = WMS_harmonic(sig,t,fs,fc,fm2,1,filterOrder);
[~,X1f_bg,Y1f_bg]  = WMS_harmonic(bg_sig,t,fs,fc,fm2,1,filterOrder);

[~,X2f,Y2f]  = WMS_harmonic(sig,t,fs,fc,fm2,2,filterOrder);
% [~,X4f,Y4f]  = WMS_harmonic(sig,t,fs,fc,fm2,4,filterOrder);

[~,X2f_bg,Y2f_bg]  = WMS_harmonic(bg_sig,t,fs,fc,fm2,2,filterOrder);
% [~,X4f_bg,Y4f_bg]  = WMS_harmonic(bg_sig,t,fs,fc,fm2,4,filterOrder);

R1f = sqrt(X1f.^2 + Y1f.^2);
R1f_bg = sqrt(X1f_bg.^2 + Y1f_bg.^2);
WMS_2f1f_L2 = sqrt( ((X2f./R1f) - (X2f_bg./R1f_bg)).^2 +  ((Y2f./R1f) - (Y2f_bg./R1f_bg)).^2 );

WMS_2f_L2 = sqrt( (X2f - X2f_bg).^2 + (Y2f - Y2f_bg).^2 );
% WMS_4f_L2 = sqrt( (X4f - X4f_bg).^2 + (Y4f - Y4f_bg).^2 );
% WMS_4f2f_L2 = WMS_4f_L2./WMS_2f_L2;

%% %%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

m = 20; %downsample data
Nsubplot = 2;

figure;
hAxis(1) = subplot(Nsubplot,1,1);
plot(t(1:m:end) * 10^6 ,WMS_2f1f_L1(1:m:end));
box off
ylabel('1391.7 nm')

hAxis(2) = subplot(Nsubplot,1,2);
plot(t(1:m:end) * 10^6 ,WMS_2f1f_L2(1:m:end),'r')
box off
ylabel('1469.3 nm')
linkaxes(hAxis,'x')
xlabel('time [$\mu$s]')


figure;
hAxis(1) = subplot(Nsubplot,1,1);
plot(t(1:m:end) * 10^6 ,WMS_2f_L1(1:m:end));
box off
ylabel('WMS-2$f$')

hAxis(2) = subplot(Nsubplot,1,2);
plot(t(1:m:end) * 10^6 ,WMS_4f_L1(1:m:end),'r')
ylim([0,0.005])
box off
ylabel('WMS-4$f$')
xlabel('time [$\mu$s]')
linkaxes(hAxis,'x')


% This script implements an algorithm based on translating the background
% intensity signal and computes the root-mean-squared-deviation. The 
% signals are aligned at the minimum RMSD. The script assumes that both
% modulation frequencies are divisible by fscan such every scan is periodic 
% and therefore only one period of the scanning frequency needs to be 
% searched to align the signals. 
% 
% This script loads in the sample data "raw_data.mat" and saves the
% "aligned_data.mat" file which is used by "plot_WMS_Harmonics_LPF.m" to
% extract the background subtracted WMS harmonics.
% 
% The script crops the background singal here to 26 scanning periods but
% this number can be chosen to be any length greater than 2 scanning
% periods upto the total length of the absorption signal. 
% 
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu


clc; clear; close all;

%Load Sample Data
load('raw_data.mat')
t = raw.t;
sig = raw.sig;
bg_sig = raw.bgsig;
fscan = raw.fscan;
fs = raw.fs;

Nsamples = round(fs/fscan); % samples per cycle
bg_sig = bg_sig(1:100*Nsamples); %crop the background signal

%make copies of the raw signals for saving
sig_copy = sig;
bg_sig_copy = bg_sig;

% plot initial unaligned signals
% figure;
% plot(sig(1:Nsamples))
% hold on
% plot(bg_sig(1:Nsamples))
% hold off
% box off
% title('Raw')
% ylabel('Signal [V]')

%Normalize data
sig = sig/max(sig(1:Nsamples));
bg_sig = bg_sig/max(bg_sig(1:Nsamples));

%plot initial unaligned signals that are now normalized
figure;
plot(sig(1:Nsamples))
hold on
plot(bg_sig(1:Nsamples))
hold off
box off
title('Normalized')


%% Find minimum residual

% This section of the code assumes that the modulation frequencies are both
% divisible by fscan, that is mod(fm1,fscan) = 0 and mod(fm2,fscan) = 0
% such that each scan is periodic. This allows us to align the signals by
% translating the background signal by "Nsamples" and locating the minium
% residual
%
% since the background signal is translated by Nsamples and the minium
% residual may occur at the end, Nsamples of the background signal are not
% useable for saving

Nbg = length(bg_sig) - Nsamples; %effective length of background signal
Nsig = length(sig);

P = floor(Nsig/Nbg); %Number of times we need to align the bg signal
R = Nsig - P*Nbg;    %Number of samples in remainder

if R == 0 %if the remainder is zero
    N = P;
else
    N = P + 1;
end

sig_bg  = zeros(Nsig,1);

if R ~= 0
    disp(['R = ',num2str(R)])
end

for i = 1 : N

    if i < N
        iL = Nbg*(i-1)+1;
        iR = i*Nbg;
    else
        iL = Nbg*(i-1)+1;
        iR = Nsig;

        if R ~= 0 
            Nbg = R;
            disp('Length of last section is R')
        end

    end

    RMSD = zeros(1,Nsamples);
    for j = 1:Nsamples
        SSR = sum((sig(iL:iR) - bg_sig([1:Nbg]+j)).^2);
        RMSD(j) = sqrt(1/Nbg * SSR);
    end

    [~,Nshift] = min(RMSD);
    ratio =  max(sig_copy(iL:iR)) / max(bg_sig_copy((1:Nbg)+Nshift));
    sig_bg(iL:iR) = bg_sig_copy((1:Nbg)+Nshift) * ratio; % scale the intensity

end

%% Aligned signals

figure;
hold on
plot(RMSD)
plot(Nshift,RMSD(Nshift),'ro','LineWidth',2)
xlabel('Sample Offset')
ylabel('RMSD')

dat = [t,sig_copy,sig_bg];

Nsig = round(10*fs/fscan);
figure;
plot(t(1:Nsig),sig_copy(1:Nsig)./max(sig_copy(1:Nsig)))
hold on
plot(t(1:Nsig),sig_bg(1:Nsig)./max(sig_bg(1:Nsig)))
hold off
legend('sig','bg sig')
xlabel('samples')

save('aligned_data.mat','dat')

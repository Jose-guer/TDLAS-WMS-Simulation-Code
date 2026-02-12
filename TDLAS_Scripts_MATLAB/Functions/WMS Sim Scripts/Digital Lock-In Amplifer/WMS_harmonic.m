function [WMS_nf,X_nf,Y_nf] = WMS_harmonic(S,t,fs,fc,fm,harmonic,filterOrder)
% This function implements a digital lock-in amplifer and returns the WMS
% harmonic signals extracted from a signal S
%
% Inputs
% S = detector signal
% t = time vector - seconds
% fs = sampling frequency - Hz
% fc = LPF cut off frequency (fpass) - Hz
% fm = modulation frequency - Hz
% harmonic = fn = n*fm, modulation frequency harmonic, n = 1, 2, 3, ...
% filterOrder = butterworth filter order, 5th order generally works well
%
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

n = harmonic;
w = 2*pi*fm;

%X-component
X = S.*cos(n*w*t); 
X_nf = filter_butter(X,fs,fc,filterOrder);

%Y-component
Y = S.*sin(n*w*t); 
Y_nf = filter_butter(Y,fs,fc,filterOrder);

%nf-harmonic magnitude
WMS_nf = sqrt(X_nf.^2 + Y_nf.^2);

end
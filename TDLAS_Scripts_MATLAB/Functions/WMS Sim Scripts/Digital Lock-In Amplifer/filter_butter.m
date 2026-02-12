function filteredSignal = filter_butter(inputSignal,fs,fc,N)
% This function filters an inputSignal using an Nth-order lowpass
% butterworth filter.
%
% Inputs
% fs = sampling frequency, Hz
% fc = cut off frequency, Hz
% N = filter order
%
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu


%butter is part of the signal processing toolbox and filter is built in
[b, a] = butter(N, fc/(fs/2)); %create a butterworth filter

filteredSignal = filter(b, a, inputSignal);

end

function Hk = wms_Hk(t,tau,k)
% This function returns the k-th Fourier coefficient
%
% Inputs
% t - time vector
% tau - transmittance = exp(-absorption)
% k - harmonic number
%
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

T = t(end)-t(1); %modulation period
w = 2*pi/T;

if k == 0
    Hk = 1/T*trapz(t,tau);
else
    Hk = 2/T*trapz(t,tau.*cos(k*w*t));
end

end
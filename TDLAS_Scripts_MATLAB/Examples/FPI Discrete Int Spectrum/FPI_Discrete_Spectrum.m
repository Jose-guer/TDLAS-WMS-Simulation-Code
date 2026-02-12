% This script simulates the discrete intensity spectrum of modulated light 
% following the theory outlined in :
% 
% Mathews, Garrett, and Christopher Goldenstein. 
% "Near-GHz scanned-wavelength-modulation spectroscopy for MHz thermometry 
% and H2O measurements in aluminized fireballs of energetic materials." 
% Applied Physics B 126.11 (2020): 189.
% 
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

close all; clear; clc;

%1469 nm Laser
lambda = 1469.28;
fscan = 500*10^3;
fmod = 45.5*10^6;

%Laser characteriation
a = 0.0795;
aHz = 2.385*10^9;
i0 = 0.175;
i2 = 8.394*10^-3;
psi = 5.88;
psi2 = 4.22;

c = 3*10^10; %cm/s
vc = 10^7/lambda; %cm^-1

M = i0/2;
beta = aHz/fmod;

K = 500;
n = 1;

for k = -K:K
    
    v(n) = vc + k*fmod/c;
    
    Jk = besselj(k,beta);
    Jkn = besselj(k-1,beta);
    Jkp = besselj(k+1,beta);
    
    E(n) = Jk + M/(2*1i)*exp(1i*psi)*Jkn - M/(2*1i)*exp(-1i*psi)*Jkp;
    n = n+1;
    
end

F = c*v/10^9;
F = F - mean(F);

intensity = abs(E).^2;
Im = max(intensity);
figure; stem(F,intensity/Im,'k','LineWidth',1,'Marker','none')
box off
xlabel('Relative frequency [GHz]')
ylabel('Discrete spectrum')
xlim([-3,3])



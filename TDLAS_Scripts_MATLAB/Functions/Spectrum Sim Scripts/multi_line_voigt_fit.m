function Resid =  multi_line_voigt_fit(vo,A,wL,T,wavenum,abs_meas)
% This function reconstructs the absorption spectrum using the input
% parameters and returns the residual 
%
% vo : vector of line centers to fit a Voigt profile
% A : Areas of voigt profiles for lines in vo
% wL : collisional broadening (HWHM) of lines
% T : temperature
% wavenum : wavenumber vector for measured absorption spectrum
% abs_meas : measured absorption spectrum
%
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

MW_rad = 18.01528; %H2O
wG  = (3.58115e-7)*vo.*sqrt(T./MW_rad); %Doppler broadening HWHM

VP = zeros(length(wavenum),length(vo)); %Voigt profile
for j = 1:length(vo)
    VP(:,j) = Voigt(wavenum, vo(j), wG(j), wL(j));
end

abs = A.*VP;

if length(vo) == 1
    abs = abs';
else
    abs = sum(abs'); %absorption spectrum
end

Resid = (abs_meas - abs);

end
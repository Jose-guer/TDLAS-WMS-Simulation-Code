function [X, absorption] = HITRAN_Spectrum_Plot(T,P,Xrad,L,dwn,wavenumLow,wavenumHigh,species)
% HITRAN_Spectrum_Plot(T [K], P [atm], Xrad, L [cm], dwn [cm^-1], wavenumLow [cm^-1], wavenumHigh [cm^-1], species)
%
% Example function call
% [X, absorption] = HITRAN_Spectrum_Plot(1200, 1, 0.15 , 10 , 0.001, 7185.2, 7186, 'H2O');
%
% This function simulates the absorbance spectrum and returns a wavenumber
% vector and the absorbance spectrum.
%
% Inputs
% T - Temperature [K]
% P - Pressure  [atm]
% Xrad - Radiating species mole fraction
% L - Path length [cm]
% dwn - wavenumber spacing [cm^-1]
% wavenumLow - Lower wavenumber bound [cm^-1]
% wavenumHigh - Upper wavenumber bound [cm^-1]
% species - Molecule name
%
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu


if strcmp('H2O',species) == 1
    qfilename = 'q_h2o_2020.txt';
    htd_filename = 'h2o_hitran2020.mat';
    MW_rad = 18.01528; %g/mole

elseif strcmp('CO',species) == 1
    qfilename = 'q_co_2020.txt';
    htd_filename = 'co_hitran2020.mat';
    MW_rad = 27.994915; %g/mole

elseif strcmp('CO2',species) == 1
    qfilename = 'q_co2_2020.txt';
    htd_filename = 'co2_hitran2020.mat';
    MW_rad = 43.989830; %g/mole

elseif strcmp('13CO2',species) == 1
    qfilename = 'q_13co2_2020.txt';
    htd_filename = '13co2_hitran2020.mat';
    MW_rad = 44.993185; %g/mole

elseif strcmp('CH4',species) == 1
    qfilename = 'q_ch4_2020.txt';
    htd_filename = 'ch4_hitran2020.mat';
    MW_rad = 16.031300; %g/mole

elseif strcmp('13CH4',species) == 1
    qfilename = 'q_13ch4_2020.txt';
    htd_filename = '13ch4_hitran2020.mat';
    MW_rad = 17.034655; %g/mole

elseif strcmp('NO',species) == 1
    qfilename = 'q_no_2020.txt';
    htd_filename = 'no_hitran2020.mat';
    MW_rad = 29.9979895; %g/mole

elseif strcmp('NO2',species) == 1
    qfilename = 'q_no2_2020.txt';
    htd_filename = 'no2_hitran2020.mat';
    MW_rad = 45.992904; %g/mole
    
elseif strcmp('H2O_HiTemp',species) == 1
    qfilename = 'q_h2o_2020.txt';
    htd_filename = 'H2O_2000_2015_HITEMP2010.mat';
    MW_rad = 18.01528; %g/mole

elseif strcmp('CO_HiTemp',species) == 1
    qfilename = 'q_co_2020.txt';
    htd_filename = 'CO_2000_2015_HITEMP2019.mat';
    MW_rad = 27.994915; %g/mole

elseif strcmp('CO2_HiTemp',species) == 1
    qfilename = 'q_co2_2020.txt';
    htd_filename = 'CO2_2380_2400_HITEMP2010.mat';
    MW_rad = 43.989830; %g/mole
end

load(htd_filename);
Part_Func = qfilename;
vo = htd.transitionWavenumber; %range of wavenumbers for which data is avaliable

wnL = wavenumLow - 0.5;
wnH = wavenumHigh + 0.5;

jj = vo >= wnL & vo <= wnH;

%Data from file that falls within the specified range
vo_vac = htd.transitionWavenumber(jj); %%vaccuum wavenumber [cm^-1]
STo = htd.lineIntensity(jj); % spectral line intensity cm^-1/(molecule-cm^-2)
Es = htd.lowerStateEnergy(jj); %lower state energy [cm^-1]
gamma_self = htd.selfBroadenedWidth(jj); %self broadening halfwidth [cm^-1/atm]
gamma_air = htd.airBroadenedWidth(jj); %air broadening halfwidth [cm-1/atm]
n_air = htd.temperatureDependence(jj); %Temperature dependence exponent
d_air = htd.pressureShift(jj); %air pressure shift [cm^-1/atm]


n_self = 0.75 * ones(size(gamma_self));

if strcmp('H2O',species) == 1
    %use measured values for these two transitions
    if vo_vac(1) > 7000
        if 7185.596 > vo_vac(1) && 7185.596 < vo_vac(end)
            [~,ind] = sort(abs (vo_vac - 7185.596));
            n_self(ind(1:2)) = 0.53;
        end
    elseif vo_vac(1) < 7000
        if 6806.03 > vo_vac(1) && 6806.03 < vo_vac(end)
            [~,ind] = sort(abs (vo_vac - 6806.03));
            n_self(ind(1:2)) = 0.86;
        end
    end
end

%% computes the Voigt profile and Linestrength for a given transition

%%% Physical constants
h = 6.6260680e-34;  % Planck constant, J-s
kB = 1.3806503e-23; % Boltzmann constant, J/K
c = 2.99792458e10; % Speed of light, cm/s
Tref = 296;  % reference temperature [K]

warning('off','all')

%Partition Function
dat = load(Part_Func);
pn = polyfit(dat(:,1),dat(:,2),6);
Q = @(T) polyval(pn,T);

%%

Xair = 1 - Xrad;

%Line Stength, cm^-1/(molecule-cm^-2)

% ST = STo .* (Q(Tref)./Q(T)).*exp(-h*c*Es/kB.*(1./T - 1/Tref));

ST = STo .* (Q(Tref)./Q(T)).*exp(-h*c*Es/kB.*(1./T - 1/Tref)).*...
    (1 - exp(-h*c*vo_vac./(kB*T)))./(1 - exp(-h*c*vo_vac/(kB*Tref)));

%pressure-shifted line center
vo_vac = vo_vac + P*Xair*d_air*(Tref./T)^0.96;

%Doppler Broadening HWHM
wG  = (3.58115e-7)*vo_vac.*sqrt(T./MW_rad);

%Collisional broadening HWHM
wL  = P.*(Xair.*gamma_air.*(Tref./T).^n_air + Xrad.*gamma_self.*(Tref./T).^n_self);

%wavenumber range
X = wnL:dwn:wnH;
SVP = zeros(size(X));

for j = 1:length(vo_vac)
    [VP] = Voigt(X, vo_vac(j), wG(j), wL(j)); % Voigt profile
    SVP = SVP + ST(j)*VP'; % Linestrength * Voigt Profile
end

absorption = SVP*((P*Xrad)/(T*kB))*L*(.101325);  % P/(KB*T) = number density

ind = X >= min(wavenumLow) - 0.02 & X <= max(wavenumHigh) + 0.02;

X = X(ind);
absorption = absorption(ind);


end
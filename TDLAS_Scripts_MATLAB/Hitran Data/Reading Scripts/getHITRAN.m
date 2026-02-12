function [vo_vac,Es,XS,So,gamma_air,gamma_self,n_air,d_air] = getHITRAN(T,Xrad,wavenumLow,wavenumHigh,threshold,species,wn_units,S_units)
%Inputs
% T - Temperature [K]
% Xrad - Radiating Species Mole Fraction
% wnL - Lower Bound Wavenumber [cm^-1]
% wnH - Upper Bound Wavenumber [cm^-1]
% threshold - Linestrength minimum value [cm^-2/atm]
% species - String of species name
% units - 1 for micromenter, 0 for wavenumbers [cm^-1]
%
%Outputs
% vo_vac - vaccuum wavenumber [cm^-1] or wavelength [um]
% Es - Lower state engery cm^-1
% XS - Linestrength at the species temperature scaled by Xrad
% So - Reference linestrength (T = 296 K) [cm^-2/atm]
% gamma_air/self - Air and self broadening parameters
% n_air - temperature dependence exponent
% d_air - air pressure shift parameter
%
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

%%% Physical constants
h = 6.6260680e-34;  % Planck constant, J-s
k = 1.3806503e-23; % Boltzmann constant, J/K
c = 2.99792458e10; % Speed of light, cm/s
Tref = 296; %[K]

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
vo = htd.transitionWavenumber; %Range of wavenumber for which data is avaliable

%wavenumber window
ind = vo >= wavenumLow & vo <= wavenumHigh;

%Data from file that falls within desired window
vo_vac = htd.transitionWavenumber(ind); %vacuum wavenumbers between the two bounds [cm^-1]... wavescan of data
So = htd.lineIntensity(ind); % spectral line intensity cm^-1/(molecule-cm^-2) (Also linestrength)
gamma_self = htd.selfBroadenedWidth(ind); %self broadened halfwidth [cm^-1/atm]
gamma_air = htd.airBroadenedWidth(ind); %air broadened halfwidth [cm-1/atm]
n_air = htd.temperatureDependence(ind); %Temperature dependence coefficient
Es = htd.lowerStateEnergy(ind); %lower state energy [cm^-1]
d_air = htd.pressureShift(ind); %air pressure shift [cm^-1/atm]

dat = load(qfilename);
pn = polyfit(dat(:,1),dat(:,2),6);
Q = @(T) polyval(pn,T);


if T == 296
    S = So;
else
    S = So .* (Q(Tref)./Q(T)).*exp(-h*c*Es/k.*(1./T - 1/Tref)).*...
        (1 - exp(-h*c*vo_vac./(k*T)))./(1 - exp(-h*c*vo_vac/(k*Tref))); %HITRAN Def.
end

if S_units == 1 %cm^-2/atm
    S = S*(7.34*10^21)./T;
    So = So*7.34*10^21./Tref;
end

if wn_units == 0
    vo_vac = 10^4./vo_vac; %units = micrometers
end

XS = Xrad*S;

ind = XS > threshold; %threshold

vo_vac = vo_vac(ind); %units = cm^-1
Es = Es(ind);
XS = XS(ind);
So = So(ind);
gamma_air = gamma_air(ind);
gamma_self = gamma_self(ind);
n_air = n_air(ind);
d_air = d_air(ind);

end


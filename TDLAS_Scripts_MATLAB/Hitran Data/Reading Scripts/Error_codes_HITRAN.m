%Hitran parameter uncertainty codes

% Code   Relative Uncertainty range
% 0      Not Reported
% 1      Default or constant
% 2      Average or estimate
% 3      >= 20%
% 4      >= 10% and < 20%
% 5      >= 5% and < 10%
% 6      >= 2% and < 5%
% 7      >= 1% and < 2%
% 8      < 1%

clc; 
clear;
close all;

%% wavenumber range

wnL = 7185.59;
wnH = 7185.6;
species = 'H2O';


%%

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

end

load(htd_filename);

vo = htd.transitionWavenumber; %Range of wavenumber for which data is avaliable
ST0 = htd.lineIntensity;
E = htd.lowerStateEnergy;
err = htd.errorCodes;           %this is an array of characters

jj = vo >= wnL & vo <= wnH;
vo_vac = vo(jj);
ST0 = ST0(jj);
E = E(jj);
err = err(jj,:); %error codes for each vo is a 1x6 char

Tab = zeros(length(vo_vac), 7);

for j = 1:length(vo_vac) %convert each digit in the character string to a number
    err_block = [str2double(err(j,2)), str2double(err(j,3)), str2double(err(j,4)), str2double(err(j,5))];

    Tab(j,:) = [vo_vac(j), ST0(j), E(j) err_block];

end

format long g
Err_codes = array2table(Tab,...
    'VariableNames',{'vo','S0','E','err-S(To)','err-g-air','err-g-self','err-n-air'})


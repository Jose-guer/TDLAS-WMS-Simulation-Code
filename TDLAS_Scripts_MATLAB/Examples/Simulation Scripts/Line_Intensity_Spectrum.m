% This script simulates the line intensity spectrum at a specified
% temperature and species mole fraction
% 
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

clear; clc;
close all;

T = 1800; %K
Xrad = 1;
threshold = 10^-4;
wavenumLow = 10^4/6.0; %wavenumber low
wavenumHigh = 10^4/1.0; %wavenumber high
wn_units = 1;
S_units = 1;


%% H2O 

species = 'H2O';

figure; set(gcf,'Position',[302 288 1194 498])
[vo,~,XS] =  getHITRAN(T,Xrad,wavenumLow,wavenumHigh,threshold,species,wn_units,S_units);
stem(vo,XS,'b','LineWidth',1.0,'Marker','none');
box off
set(gca,'yscal','log')
legend(species,'Location','Northwest','box','on')
ylabel('$X_i \cdot S(T) ~ \mathrm{[cm^{-2}/atm]}$')

if wn_units == 1
    xlim([wavenumLow,wavenumHigh]);
    xlabel('wavenumber $\mathrm{[cm^{-1}]}$')

elseif wn_units == 0
    xlim([1,6])
    xlabel('wavelength $(\mu$m)')
end


%% CO2


figure; set(gcf,'Position',[302 288 1194 498])
hold on;

species = 'CO2';
[vo,~,XS] =  getHITRAN(T,Xrad,wavenumLow,wavenumHigh,threshold,species,wn_units,S_units);
stem(vo,XS,'b','LineWidth',1.0,'Marker','none');


species = '13CO2';
[vo,~,XS] =  getHITRAN(T,Xrad,wavenumLow,wavenumHigh,threshold,species,wn_units,S_units);
stem(vo,XS,'r','LineWidth',1.0,'Marker','none');

box off
set(gca,'yscal','log')
legend('CO2','13CO2','Location','Northwest','box','on')
ylabel('$X_i \cdot S(T) ~ \mathrm{[cm^{-2}/atm]}$')

if wn_units == 1
    xlim([wavenumLow,wavenumHigh]);
    xlabel('wavenumber $\mathrm{[cm^{-1}]}$')

elseif wn_units == 0
    xlim([1,6])
    xlabel('wavelength $(\mu$m)')
end


% This script plots the normalized linestrength function S(T)/S(T0) over a
% range of temperatures for specified lower-state energies (Eg)
% 
% written by Jose Guerrero, University of Michigan - Aerospace Department 
% joseguer@umich.edu

clc; clear;
close all;

%% Simulation details

%Molecule
qfilename = 'q_h2o_2020.txt';

%Lower state energy values
Eg = [1045,3291]; 

%Temperature range
Trange = 300:10:3000;

%Plot configurations
ymax = 0;      %scale by peak amplitude
logplot = 1;   %use semilogy
S_P_units = 1; % 0 = cm^-1/(molecule - cm^-2), 1 = cm^-2/atm

%%

%Constants
h = 6.6260680e-34;  % Planck constant, J-s
k = 1.3806503e-23; % Boltzmann constant, J/K
c = 2.99792458e10; % Speed of light, cm/s
Tref = 296;

%Partition Function
Part_Func = qfilename;
dat = load(Part_Func);
pn = polyfit(dat(:,1),dat(:,2),6);
Q = @(T) polyval(pn,T);

LineStrengthRatio = zeros(length(Eg),length(Trange));

legstr = {};

for i = 1:length(Eg)

    if S_P_units == 1
        %cm^-2/atm
        LineStrengthRatio(i,:) =  (Q(Tref)./Q(Trange)).*Tref./Trange.*exp(-h*c*Eg(i)/k.*(1./Trange - 1/Tref));
    else
        %cm^-1/molec-cm^-2
        LineStrengthRatio(i,:) =  (Q(Tref)./Q(Trange)).*exp(-h*c*Eg(i)/k.*(1./Trange - 1/Tref));
    end

    legStr(i) = {['$E''''$ = ',num2str(Eg(i)),' $\mathrm{cm^{-1}}$']};
end

%% Plots

Smax = max(LineStrengthRatio');

figure; set(gcf,'position',[377 319 690 567]);

if ymax == 1
    if logplot == 1
        semilogy(Trange,LineStrengthRatio./Smax');
    else
        plot(Trange,LineStrengthRatio./Smax');
    end
else
    if logplot == 1
        semilogy(Trange,LineStrengthRatio);
    else
        plot(Trange,LineStrengthRatio);
    end
end
hold on

xline(1000,'k--','LineWidth',1)
xline(2500,'k--','LineWidth',1)
hold off
xlabel('T [K]')
legend(legStr,...
    'box','on','Location','southwest')
grid on; grid minor;
box off
xlim([min(Trange),max(Trange)])

if S_P_units == 1
    ylabel('$S(T)/S(T_0)$')
else
    ylabel('$S^*(T)/S^*(T_0)$')
end





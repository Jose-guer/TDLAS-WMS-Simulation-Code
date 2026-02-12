% This script implements a multi-line Voigt profile fitting routine to infer
% gas properties from noisy synthetic absorption spectra. The algorithm
% implements a nonlinear fitting routine where A, wL, and vo are varried.
% The fitting routine is initiated by guessing an inital temperature. The
% fitting routine will loop until the inferred temperature from the best-fit
% parameters is within a specified tolerance of the current iterations
% guess.
% 
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu

clc; clear;
close all;

load('data_noise.mat')
% load('data_noise_free.mat')

x1 = dat.x1;
y1 = dat.absor1;
x2 = dat.x2;
y2 = dat.absor2;

P = dat.P;
T = dat.T;
X = dat.X;
L = dat.L;

figure;
subplot(2,1,1)
plot(x1,y1)
box off
ylabel('$\alpha(\nu)$')

subplot(2,1,2)
plot(x2,y2,'r')
box off
xlabel('wavenumber [$\mathrm{cm^{-1}}$]')
ylabel('$\alpha(\nu)$')


%% %%% Infer gas properties

%using the sum of the doublets linestrength values in the Hitran database
Es1 = 1045.058; %cm^-1
STo1 = 0.01962; %cm^-2/atm

STo2 = 6.55e-07; %cm^-1
Es2 = 3291.15; %cm^-2/atm

Tsim = 1000:100:2500; %K

qfilename = 'q_h2o_2020.txt';
h = 6.6260680e-34;  % Planck constant, J-s
kB = 1.3806503e-23; % Boltzmann constant, J/K
c = 2.99792458e10; % Speed of light, cm/s
Tref = 296; %K
MW_rad = 18.01528;

%Partition Function
Part_Func = qfilename;
dat = load(Part_Func);
pn = polyfit(dat(:,1),dat(:,2),6);
Q = @(T) polyval(pn,T);

%pressure normalized definition of the linestrength function
ST1 = STo1.* (Q(Tref)./Q(Tsim)).*(Tref./Tsim).*exp(-h*c*Es1/kB.*(1./Tsim - 1/Tref));
ST2 = STo2.* (Q(Tref)./Q(Tsim)).*(Tref./Tsim).*exp(-h*c*Es2/kB.*(1./Tsim - 1/Tref));
R = ST1./ST2;

%%

%initial guesses 

Tguess = 1000;

vo1 = [7185.4,7185.596];
A1 = [0.0015,0.015];
wG1  = (3.58115e-7)*vo1.*sqrt(Tguess./MW_rad);
wL1 = 2*wG1;

vo2 = [6805.81,6806.033,6806.12];
A2 = [0.001,0.0047,0.0008];
wG2  = (3.58115e-7)*vo2.*sqrt(Tguess./MW_rad);
wL2 = 2*wG2;

errT = 100;
options =  optimset('display','off','Algorithm','levenberg-marquardt');

while errT > 0.1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    X0 = [vo1,A1,wL1];
    lb = [vo1 - 0.005,[A1,wL1]*0.5];
    ub = [vo1 + 0.005,[A1,wL1]*2];

    fit_param1 = lsqnonlin(@(X) multi_line_voigt_fit([X(1),X(2)],[X(3),X(4)],...
        [X(5),X(6)],Tguess,x1,y1),X0,lb,ub,options);

    vo1 = [fit_param1(1),fit_param1(2)];
    A1 = [fit_param1(3),fit_param1(4)];
    Area1 = A1(2);
    wL1 = [fit_param1(5),fit_param1(6)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    X0 = [vo2,A2,wL2];
    lb = [vo2 - 0.005,[A2,wL2]*0.5];
    ub = [vo2 + 0.005,[A2,wL2]*2];

    fit_param2 = lsqnonlin(@(X) multi_line_voigt_fit([X(1),X(2),X(3)],...
        [X(4),X(5),X(6)],[X(7),X(8),X(9)],Tguess,x2,y2),X0,lb,ub,options);

    vo2 = [fit_param2(1),fit_param2(2),fit_param2(3)];
    A2 = [fit_param2(4),fit_param2(5),fit_param2(6)];
    Area2 = A2(2);
    wL2 = [fit_param2(7),fit_param2(8),fit_param2(9)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Tmeas = interp1(R,Tsim,Area1/Area2);
    errT = abs(Tmeas-Tguess)/Tguess*100;

    wG1  = (3.58115e-7)*vo1.*sqrt(Tmeas./MW_rad);
    wG2  = (3.58115e-7)*vo2.*sqrt(Tmeas./MW_rad);

    %update guessed temperature
    Tguess = Tmeas;

end

%Species partial pressure
S1 = interp1(Tsim,ST1,Tmeas);
XH2O = Area1/(S1*P*L);

%%

%Reconstruct absorption spectra from best-fit parameters
%7185.596
VP = zeros(length(x1),length(vo1));
for j = 1:length(vo1)
    VP(:,j) = Voigt(x1, vo1(j), wG1(j), wL1(j));
end

abs_fit1 = A1.*VP;
abs_fit1 = sum(abs_fit1');

%6806.03
VP = zeros(length(x2),length(vo2));
for j = 1:length(vo2)
    VP(:,j) = Voigt(x2, vo2(j), wG2(j), wL2(j));
end

abs_fit2 = A2.*VP;
abs_fit2 = sum(abs_fit2');

% Plot results

figure
subplot(2,1,1)
plot(x1,y1,'k*','MarkerSize',1)
hold on
plot(x1,abs_fit1,'b','linewidth',2)
hold off
ylabel('$\alpha(\nu)$')
box off

subplot(2,1,2)
plot(x2,y2,'k*','MarkerSize',1)
hold on
plot(x2,abs_fit2,'r','linewidth',2)
hold off
ylabel('$\alpha(\nu)$')
box off
xlabel('wavenumber [$\mathrm{cm^{-1}}$]')


errT = (T-Tmeas)/T * 100;
errX = (X-XH2O)/X * 100;

A = [Tmeas,errT,XH2O,errX];

Tab = array2table(A,...
    'VariableNames',{'T [K]','err T [%]','X','err X [%]'})


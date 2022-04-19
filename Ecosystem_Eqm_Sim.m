clear
clc

% This code considers three generic outcomes. The first is where there are 
% priority effects. Annual and perennial have differences in litter production. 

%% Generic Parameters
% Annual fitness parameters
R0A = 48;

gA = 1;
sA = 0.5;
eA = 1;
yA = R0A*(1-sA*(1-gA))/(gA*eA);

% Perennial fitness parameters
R0P = 23;

gP = 0.5;
sP = 0.5;
eP = 1;
f = 0.5;

T_P = 4; % <-- T_P = pS/(1-pP), the expected lifespan of a perennial plant
         % Restriction is that T_P > p_S
pS = 0.1;
pP = 1 - pS/T_P;
yP = R0P*(1 - sP*(1-gP))/(gP*eP*(f + pS/(1-pP)));

% Competition parameters
alphaA = 0.1;   alphaPprime = 0.1;
% Here we fix alphaPprime, which is a composite parameter that acts like
% alphaA, and then choose values of alphaP and gamma that give alphaPprime. 
alphaP = 0.1; gamma = (alphaPprime/alphaP - 1)*T_P;

% Litter sensitivity parameters
betaA = 0.13; betaP = 0.05;

% Ecosystem parameters
d = linspace(0.001,1,500); 
bT = 0;

%% bA >> bP (Broad region of stable coexistence)

% Litter production parameters
bA = 0.5; bPprime = 0.1; 
% Here we choos bPprime, which is a composite parameter that acts just like
% bA. Once we have chosen bPprime, we choose values of bP and delta that
% give bPprime (once already choosing T_P above).
bP = 0.1; delta = bPprime/bP - 1/T_P;

% Calculate the value of litter where species coexist.
LeqmCoex = (R0A - R0P)/(betaA*R0P - betaP*R0A);

% Calculate the critical points of decomposition where invasion of
% different species occurs. 
dcritP = 1/LeqmCoex*(bT + bA/alphaA*(R0A/(1 + betaA*LeqmCoex) - 1));
dcritA = 1/LeqmCoex*(bT + bP*(delta + 1/T_P)/(alphaP*(1+gamma/T_P))*(R0P/(1+betaP*LeqmCoex) - 1));

% Given that betaA > betaP, 
%   the perennial can invade when d < dcritP
%   the annual can invade when d > dcritA.
% Therefore, the species coexist over some range of d if dcritP > dcritA
% This is the case we consider here.

% Plot the 3 outcomes.
% 1. d < dcritA < dcritP (perennial only)
dvec1 = linspace(0,dcritA,100);
for i = 1:length(dvec1)
    Leqm_NoA(i) = LitterEq(R0P, alphaPprime, betaP, bPprime/dvec1(i), bT/dvec1(i));
end
figure
subplot(1,3,3)
plot(dvec1, Leqm_NoA, 'Color', 'black', 'LineWidth', 3);
hold on

% 2. dcritA < d < dcritP (coexistence)
dvec2 = linspace(dcritA,dcritP,100);
plot(dvec2, LeqmCoex*ones(size(dvec2)), 'Color', 'black', 'LineWidth', 3);

% 3. d > dcritP > dcritA (annual only)
dvec3 = linspace(dcritP, 1, 100);
for i = 1:length(dvec1)
    Leqm_NoP(i) = LitterEq(R0A, alphaA, betaA, bA/dvec3(i), bT/dvec3(i));
end
plot(dvec3, Leqm_NoP, 'Color', 'black', 'LineWidth', 3);
hold off

xlabel('Decomposition Fraction, {\itd}')
ylabel('Equilbirium Litter Density');
xline(dcritA, '--');
xline(dcritP, '--');
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
xlim([0,1]); ylim([0,150]);


%% bA == bP (Coexistence region has measure zero)

% For illustrative purposes, make it such that the invasion boundaries are
% both at d = 0.5.

% Calculate the value of b that gives dcrit = 0.5.
dcrit = 0.5;
b = (dcrit*LeqmCoex - bT)*alphaA/(R0A/(1 + betaA*LeqmCoex) - 1);

% Litter production parameters
bA = b; bPprime = b; 
% Here we choos bPprime, which is a composite parameter that acts just like
% bA. Once we have chosen bPprime, we choose values of bP and delta that
% give bPprime (once already choosing T_P above).
bP = b; delta = bPprime/bP - 1/T_P;

% Calculate the critical points of decomposition where invasion of
% different species occurs. 
dcritP = 1/LeqmCoex*(bT + bA/alphaA*(R0A/(1 + betaA*LeqmCoex) - 1));
dcritA = 1/LeqmCoex*(bT + bP*(delta + 1/T_P)/(alphaP*(1+gamma/T_P))*(R0P/(1+betaP*LeqmCoex) - 1));

% Given that betaA > betaP, 
%   the perennial can invade when d < dcritP
%   the annual can invade when d > dcritA.
% Therefore, the species coexist over some range of d if dcritP > dcritA
% In this case, dcritP == dcritA. Therefore, we only have two outcomes

% Plot the 2 outcomes.
% 1. d < dcritA = dcritP (perennial only)
dvec1 = linspace(0,dcritA,100);
for i = 1:length(dvec1)
    Leqm_NoA(i) = LitterEq(R0P, alphaPprime, betaP, bPprime/dvec1(i), bT/dvec1(i));
end

subplot(1,3,2)
plot(dvec1, Leqm_NoA, 'Color', 'black', 'LineWidth', 3);
hold on

% 2. d > dcritP = dcritA (annual only)
dvec2 = linspace(dcritP, 1, 100);
for i = 1:length(dvec1)
    Leqm_NoP(i) = LitterEq(R0A, alphaA, betaA, bA/dvec2(i), bT/dvec2(i));
end
plot(dvec2, Leqm_NoP, 'Color', 'black', 'LineWidth', 3);
hold off

xlabel('Decomposition Fraction, {\itd}')
ylabel('Equilbirium Litter Density');
xline(dcritA, '--');
xline(dcritP, '--');
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
xlim([0,1]); ylim([0,150]);

%% bA << bP (Broad region of priority effects)

% Litter production parameters
bA = 0.1; bPprime = 0.5; 
% Here we choos bPprime, which is a composite parameter that acts just like
% bA. Once we have chosen bPprime, we choose values of bP and delta that
% give bPprime (once already choosing T_P above).
bP = 0.5; delta = bPprime/bP - 1/T_P;

% Calculate the critical points of decomposition where invasion of
% different species occurs. 
dcritP = 1/LeqmCoex*(bT + bA/alphaA*(R0A/(1 + betaA*LeqmCoex) - 1));
dcritA = 1/LeqmCoex*(bT + bP*(delta + 1/T_P)/(alphaP*(1+gamma/T_P))*(R0P/(1+betaP*LeqmCoex) - 1));

% Given that betaA > betaP, 
%   the perennial can invade when d < dcritP
%   the annual can invade when d > dcritA.
% Therefore, the species coexist over some range of d if dcritP > dcritA
% In this case, dcritA > dcritP. In the region between dcritA and dcritP,
% neither species can invade the other. This means that there is a priority
% effect. 

% Plot the 3 outcomes.
% 1. d < dcritA < dcritP (perennial only)
dvec1 = linspace(0,dcritP,100);
for i = 1:length(dvec1)
    Leqm_NoA(i) = LitterEq(R0P, alphaPprime, betaP, bPprime/dvec1(i), bT/dvec1(i));
end
subplot(1,3,1)
plot(dvec1, Leqm_NoA, 'Color', 'black', 'LineWidth', 3);
hold on

% 2. dcritP < d < dcritA (priority effect)
dvec2 = linspace(dcritP,dcritA,100);
for i = 1:length(dvec2)
    Leqm_NoP(i) = LitterEq(R0A, alphaA, betaA, bA/dvec2(i), bT/dvec2(i));
    Leqm_NoA(i) = LitterEq(R0P, alphaPprime, betaP, bPprime/dvec2(i), bT/dvec2(i));
end
plot(dvec2, Leqm_NoP, 'Color', 2/3*ones(1,3), 'LineWidth', 3);
plot(dvec2, Leqm_NoA, 'Color', 2/3*ones(1,3), 'LineWidth', 3);

% 3. d > dcritP > dcritA (annual only)
dvec3 = linspace(dcritA, 1, 100);
for i = 1:length(dvec1)
    Leqm_NoP(i) = LitterEq(R0A, alphaA, betaA, bA/dvec3(i), bT/dvec3(i));
end
plot(dvec3, Leqm_NoP, 'Color', 'black', 'LineWidth', 3);
hold off

xlabel('Decomposition Fraction, {\itd}')
ylabel('Equilbirium Litter Density');
xline(dcritA, '--');
xline(dcritP, '--');
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
xlim([0,1]); ylim([0,150]);

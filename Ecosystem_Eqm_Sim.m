clear
clc

% This script plots 

%% Consider Strong and Weak Tradeoffs

% Strong Tradeoffs

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

% Litter production parameters
bA = 0.5; bPprime = 0.1; bP = 0.1; 
delta = (T_P*bPprime/bP - 1)/(T_P+1);

% Competition parameters
alphaA = 0.1;   alphaPprime = 0.1;
alphaP = 0.1; gamma = (alphaPprime/alphaP - 1)*T_P;

% Litter sensitivity parameters
betaA = 0.13; betaP = 0.05;

% Ecosystem parameters
d = linspace(0.001,1,500); 
bT = 0;


gen = 100000;


eqm = zeros(4,length(d));

s = [sA, sP, pS, pP];   y = [yA, yP, f];    g = [gA, gP];
e = [eA, eP];           alpha = [alphaA, alphaP, gamma];
beta = [betaA, betaP];

for i = 1:length(d)
    decay = [bA, bP, d(i), bT, delta];
    parameters = {s, y, g, e, decay, alpha, beta};
    
%    if i == 1
        init_cond = ones(4,1);
%    else
%        init_cond = sys(:,end);
%    end
    
    sys = APL_Sim_Tree(gen, init_cond, parameters);
    
    eqm(:,i) = sys(:,end);
end

dcrit(1) = d(find(log(eqm(1,:)) > -2,1));
dcrit(2) = d(find(log(ones(1,2)*eqm(3:4,:)) < -2, 1));

%% Figures 

figure
subplot(1,2,2)
plot(d, eqm(2,:), 'Color', 'black', 'LineWidth', 3);
xlabel('Decomposition Fraction, {\itd}')
ylabel('Equilbirium Litter Density');
xline(dcrit, '--');
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';

for i = 1:length(d)
    Leqm_NoA(i) = LitterEq(R0P, alphaPprime, betaP, bPprime/d(i), bT/d(i));
    Leqm_NoP(i) = LitterEq(R0A, alphaA, betaA, bA/d(i), bT/d(i));
end

hold on
plot(d, Leqm_NoA)
plot(d, Leqm_NoP)
LeqmCoex = -(R0A - R0P)/(betaP*R0A - betaA*R0P);
yline(LeqmCoex)
hold off


%% Weak Tradeoffs
bPprime = 0.1;

% Litter production parameters
bA = 0.1; bP = 0.1; 
delta = (T_P*bPprime/bP - 1)/(T_P+1);

% Competition parameters
alphaA = 0.1;   alphaPprime = 0.1;
alphaP = 0.1; gamma = (alphaPprime/alphaP - 1)*T_P;

eqm = zeros(4,length(d));

s = [sA, sP, pS, pP];   y = [yA, yP, f];    g = [gA, gP];
e = [eA, eP];           alpha = [alphaA, alphaP, gamma];
beta = [betaA, betaP];

parfor i = 1:length(d)
    decay = [bA, bP, d(i), bT, delta];
    parameters = {s, y, g, e, decay, alpha, beta};
    
%    if i == 1
        init_cond = ones(4,1);
%    else
%        init_cond = sys(:,end);
%    end
    
    sys = APL_Sim_Tree(gen, init_cond, parameters);
    
    eqm(:,i) = sys(:,end);
end

%% Figures 


dcrit(1) = d(find(log(eqm(1,:)) > -2,1));
dcrit(2) = d(find(log(ones(1,2)*eqm(3:4,:)) < -2, 1));


subplot(1,2,1)
plot(d, eqm(2,:), 'Color', 'black', 'LineWidth', 3);
xlabel('Decomposition Fraction, {\itd}')
ylabel('Equilbirium Litter Density');
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
xline(dcrit, '--')

for i = 1:length(d)
    Leqm_NoA(i) = LitterEq(R0P, alphaPprime, betaP, bPprime/d(i), bT/d(i));
    Leqm_NoP(i) = LitterEq(R0A, alphaA, betaA, bA/d(i), bT/d(i));
end

hold on
plot(d, Leqm_NoA)
plot(d, Leqm_NoP)
LeqmCoex = -(R0A - R0P)/(betaP*R0A - betaA*R0P);
yline(LeqmCoex)
hold off
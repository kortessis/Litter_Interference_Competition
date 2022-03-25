clear
clc

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
bA = 0.5; bPprime = 0.1; 
% Here we choos bPprime, which is a composite parameter that acts just like
% bA. Once we have chosen bPprime, we choose values of bP and delta that
% give bPprime (once already choosing T_P above).
bP = 0.1; delta = bPprime/bP - 1/T_P;

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


gen = 100000;

eqm = zeros(4,length(d));

LeqmCoex = (R0A - R0P)/(betaA*R0P - betaP*R0A);

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

dcritP = 1/LeqmCoex*(bT + bA/alphaA*(R0A/(1 + betaA*LeqmCoex) - 1));
dcritA = 1/LeqmCoex*(bT + bP*(delta + 1/T_P)/(alphaP*(1+gamma/T_P))*(R0P/(1+betaP*LeqmCoex) - 1));

%% Figures 

figure
subplot(1,2,2)
plot(d, eqm(2,:), 'Color', 'black', 'LineWidth', 3);
xlabel('Decomposition Fraction, {\itd}')
ylabel('Equilbirium Litter Density');
xline(dcritA, '--');
xline(dcritP, '--');
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';

 % Uncomment these lines if you want to see the resident equlibrium litter
 % for each species. They track the simulation results exactly except where
 % the two species coexist.
 
% for i = 1:length(d)
%     Leqm_NoA(i) = LitterEq(R0P, alphaPprime, betaP, bPprime/d(i), bT/d(i));
%     Leqm_NoP(i) = LitterEq(R0A, alphaA, betaA, bA/d(i), bT/d(i));
% end
% 
% hold on
% plot(d, Leqm_NoA)
% plot(d, Leqm_NoP)
% yline(LeqmCoex)
% hold off


%% Weak Tradeoffs
bPprime = 0.1;

% Litter production parameters
bA = 0.1; bP = 0.1; 
delta = bPprime/bP - 1/T_P;

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

dcritP = 1/LeqmCoex*(bT + bA/alphaA*(R0A/(1 + betaA*LeqmCoex) - 1));
dcritA = 1/LeqmCoex*(bT + bP*(delta + 1/T_P)/(alphaP*(1+gamma/T_P))*(R0P/(1+betaP*LeqmCoex) - 1));


subplot(1,2,1)
plot(d, eqm(2,:), 'Color', 'black', 'LineWidth', 3);
xlabel('Decomposition Fraction, {\itd}')
ylabel('Equilbirium Litter Density');
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
xline(dcritA, '--')
xline(dcritP, '--')

 % Uncomment these lines if you want to see the resident equlibrium litter
 % for each species. They track the simulation results exactly.
 
% for i = 1:length(d)
%     Leqm_NoA(i) = LitterEq(R0P, alphaPprime, betaP, bPprime/d(i), bT/d(i));
%     Leqm_NoP(i) = LitterEq(R0A, alphaA, betaA, bA/d(i), bT/d(i));
% end

% hold on
% plot(d, Leqm_NoA)
% plot(d, Leqm_NoP)
% yline(LeqmCoex)
% hold off
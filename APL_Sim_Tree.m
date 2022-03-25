function sys = APL_Sim_Tree(gen, init_cond, parameters)

% Function for running the dynamics of the model with one annual, one
% perennial, and litter, including the effect of litter produced by trees. 

s = parameters{1}; % Values of survival for both species. This is a 4-d 
                   % vector with entries [annual seed survival, perennial
                   % seed survival, first year survival, and perennial adult survival.         
                   %        [sA, sP, pS, pP]
sA = s(1); sP = s(2); pS = s(3); pP = s(4);
y = parameters{2}; % Values of reproduction for both species. 
                   %        [yA, yP, f]
yA = y(1); yP = y(2); f = y(3);
g = parameters{3}; % Germination when no litter is present. 
                   %        [gA, gP]
gA = g(1); gP = g(2);
e = parameters{4}; % Establishment succes when no litter is present
                   %        [eA, eP]
decay = parameters{5}; % Litter production and decay rates. 
                       %    [bA, bP, d, bT, delta];
bA = decay(1); bP = decay(2); d = decay(3); bT = decay(4); delta = decay(5);
alpha = parameters{6}; % Sensitivities to competition. 
                       %    [alphaA, alphaP, gamma]
alphaA = alpha(1); alphaP = alpha(2); gamma = alpha(3);
beta = parameters{7};  % Sensitivities to litter. 
                       %    [betaA, betaP]

% Create vectors to hold state variables
NA = zeros(1,gen);
L = zeros(1,gen);
NP = zeros(2,gen);

% Set initial conditions
NA(1) = init_cond(1);
L(1) = init_cond(2);
NP(:,1) = init_cond(3:4);

% Variables for germination and competition
E = zeros(2,gen);
C = zeros(1,gen);

% Germination and reproduction functions, accounting for both litter and
% competition.
E(:,1) = e./(1+beta*L(1));
C(1) = 1 + alphaA*E(1,1)*gA*NA(1)...
         + alphaP*gamma*gP*E(2,1)*NP(1,1)...
         + alphaP*NP(2,1);

% Loop for dynamics
for t = 2:gen
    
    % Annual Dynamics
    NA(t) = NA(t-1)*(sA*(1 - gA) + gA*E(1,t-1)*yA/C(t-1));
    
    % Litter Dynamics
    L(t) = bA*NA(t-1)*E(1,t-1)*gA + ...
        bP*(gP*E(2,t-1)*NP(1,t-1)*(pS*delta + 1 - pS) + NP(2,t-1)*(pP*delta + 1 - pP))...
        + (1-d)*L(t-1) + bT;
    
    % Perennial Dyanmics
    M = [sP*(1-gP) + E(2,t-1)*gP*yP*f/C(t-1), ...
         yP/C(t-1);...
         gP*E(2,t-1)*pS,...
         pP];
     
    NP(:,t) = M*NP(:,t-1);
    
    % Updating litter and reproduction
    E(:,t) = e./(1+beta*L(t));
    C(:,t) = 1 + alphaA*E(1,t)*gA*NA(t)...
               + alphaP*gamma*gP*E(2,t)*NP(1,t)...
               + alphaP*NP(2,t);

    
end

% Full dyanmics of four species system in a matrix of dimensions 4 x gen.
% Rows are in order
    % 1. Annual seeds
    % 2. Litter
    % 3. Perennial seeds
    % 4. Perennial adults
sys = [NA;L;NP];
end
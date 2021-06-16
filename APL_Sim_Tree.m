function sys = APL_Sim_Tree(gen, init_cond, parameters)

% Function for running the dynamics of the model with one annual, one
% perennial, and litter, including the effect of litter produced by trees. 

s = parameters{1}; % Values of survival for both species. This is a 4-d 
                   % vector with entries [annual seed survival, perennial
                   % seed survival, first year survival, and perennial adult survival.         
                   %        [sA, sP, p1, p2]
y = parameters{2}; % Values of reproduction for both species. 
                   %        [yA, yP, f]
g = parameters{3}; % Germination when no litter is present. 
                   %        [gA, gP]
e = parameters{4}; % Establishment succes when no litter is present
                   %        [eA, eP]
decay = parameters{5}; % Litter production and decay rates. 
                       %    [bA, bP, d, BT];
alpha = parameters{6}; % Sensitivities to competition. 
                       %    [alphaA, alphaP, gamma]
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
C(1) = 1 + alpha(1)*E(1,1)*g(1)*NA(1)...
         + alpha(2)*alpha(3)*g(2)*E(2,1)*NP(1,1)...
         + alpha(2)*NP(2,1);

% Loop for dynamics
for t = 2:gen
    
    % Annual Dynamics
    NA(t) = NA(t-1)*(s(1)*(1 - g(1)) + g(1)*E(1,t-1)*y(1)/C(t-1));
    
    % Litter Dynamics
    L(t) = decay(1)*NA(t-1)*E(1,t-1)*g(1) + NP(2,t-1)*decay(2) + (1-decay(3))*L(t-1) + decay(4);
    
    % Perennial Dyanmics
    M = [s(2)*(1-g(2)) + E(2,t-1)*g(2)*y(2)*y(3)/C(t-1), ...
         y(2)/C(t-1);...
         g(2)*E(2,t-1)*s(3),...
         s(4)];
     
    NP(:,t) = M*NP(:,t-1);
    
    % Updating litter and reproduction
    E(:,t) = e./(1+beta*L(t));
    C(:,t) = 1 + alpha(1)*E(1,t)*g(1)*NA(t)...
               + alpha(2)*alpha(3)*g(2)*E(2,t)*NP(1,t)...
               + alpha(2)*NP(2,t);

    
end

% Full dyanmics of four species system in a matrix of dimensions 4 x gen.
% Rows are in order
    % 1. Annual seeds
    % 2. Litter
    % 3. Perennial seeds
    % 4. Perennial adults
sys = [NA;L;NP];
end
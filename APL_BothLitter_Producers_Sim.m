function sys = APL_BothLitter_Producers_Sim(gen, init_cond, parameters)

s = parameters{1}; % Values of survival for both species. This is a 3-d 
                   % vector with entries [annual seed survival, perennial
                   % seed survival, and establishment probability of a 
                   % germinated perennial seed, and perennial adult 
                   % survival.         
                   %        [sA, sP, s]
y = parameters{2}; % Values of reproduction for both species. 
                   %        [yA, y1 = yP*f, yP]
g = parameters{3}; % Germination absent litter. 
                   %        [gA, gP]
e = parameters{4}; % Establishment succes in the absence of litter
                   %        [eA, eP]
decay = parameters{5}; % Litter decay rates. 
                       %    [bA, bP, d];
alpha = parameters{6}; % Sensitivities to competition. 
                       %    [alphaA, gamma*alphaP = alphaF, alphaP]
beta = parameters{7};  % Sensitivities to litter. 
                       %    [betaA, betaP]

% Vectors to hold densities
NA = zeros(1,gen);
L = NA;
N = zeros(2,gen);

% Initial conditions in the order of [NA(0), L(0), NS(0), NP(0)]
NA(1) = init_cond(1);
L(1) = init_cond(2);
N(:,1) = init_cond(3:4);

% Variables for germination and competition
E = zeros(2,gen);
C = zeros(1,gen);

% Germination and reproduction functions, accounting for both litter and
% competition.
E(:,1) = e./(1+beta*L(1));
C(1) = 1 + alpha(1)*g(1)*E(1,1)*NA(1) + alpha(2)*g(2)*E(2,1)*N(1,1) + alpha(3)*N(2,1);

% Loop for dynamics
for t = 2:gen
    
    % Annual Dynamics
    NA(t) = NA(t-1)*(s(1)*(1 - g(1)) + g(1)*E(1,t-1)*y(1)/C(t-1));
    % Litter Dynamics
    L(t) = decay(1)*NA(t-1)*g(1)*E(1,t-1) + N(1,t-1)*decay(2) + (1-decay(3))*L(t-1);
    
    % Perennial Dyanmics
    M = [s(2)*(1-g(2)) + g(2)*E(2,t-1)*y(2)/C(t-1), y(3)/C(t-1);...
         g(2)*E(2,t-1),                                 s(3)];
    N(:,t) = M*N(:,t-1);
    
    % Updating litter and reproduction
    E(:,t) = e./(1+beta*L(t));
    C(:,t) = 1 + alpha(1)*g(1)*E(1,t)*NA(t) + alpha(2)*g(2)*E(2,t)*N(1,t) + alpha(3)*N(2,t);

    
end

% Full dyanmics of four species system in a matrix of dimensions 4 x gen.
% Rows are in order
    % 1. Annual seeds
    % 2. Litter
    % 3. Perennial seeds
    % 4. Perennial adults
sys = [NA;L;N];
end
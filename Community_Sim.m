clear
clc

%% Plotting Parameters
x = viridis(4);     colors = x(2:3,:);
axFontSize = 30;

%% Simulation Parameters

n = 10;        % Number of perennial species in the community
gen = 3000;    % Total time steps of the simulation
tint = 1000;   % Initial time steps before annual is introduced
betaA = 1/1.3; % Sensitivity of annual to litter
betaP = 1/15;  % Sensitivity of perennials to litter
beta = [betaA; betaP*ones(n,1)]; 

rho = exp(-4); % Niche overlap among the perennials

% Mean fitness of perennials and a measure of variability across species.
lambda_mean = 5; lambda_var = 1;

% Annual Fitness
lambda_A = lambda_mean*30;

% Perennial Fitness
lambda_P = lambda_mean*exp(linspace(-lambda_var,lambda_var,n));

% Distribution of litter production
b = 1*[1, 0.01*linspace(0,1,n)];
d = linspace(0.001,1,32);

% Fitness Parameters for the Annual
eA = 1; gA = 1; sA = 0.5;
yA = lambda_A*(1 - sA*(1-gA))/(gA*eA);

% Fitness Parameters for the Perennial
eP = 0.5; gP = 0.5; sP = 0.5;
f = 0.1;  p1 = 1; p2 = 0.2;
gamma = 0.1;
yP = lambda_P*(1 - sP*(1-gP))/(gP*eP*(f + p1/(1 - p2)));

% Competitive effects;
alpha_A = 0.2;
alpha_P = 0.2;

% Species-specific competitive effects of perennials to allow for
% coexistence.
A_PP = alpha_P*(rho*ones(n) + (1-rho)*eye(n));

% Row and column vectors of the competition coefficient matrix
A_AP = alpha_A*ones(1,n); A_PA = alpha_P*ones(n,1);

% The entire competition matrisx.
A = [alpha_A, A_AP; A_PA, A_PP];

% Cacluate all possible fitness ratios for all species and only count those
% greater than 1. The distribution of ratios is symmetric on the log scale
% and so this is sufficient information to describe the entire distribution
% of fitness ratios for all species pairs.
fitness_inequal = reshape(lambda_P'*(1./lambda_P),1, []);
fitness_inequal = fitness_inequal(fitness_inequal > 1);


% Make figure S2 in the supplement
figure(1)

%subplot(1,2,1)
% Distribution of species fitnesses
scatter([lambda_A, lambda_P], b, 200, 'filled', 'CData', [zeros(1,3);viridis(n)])
xlabel('\it\lambda'); ylabel('\itb');
axis([1,150,0,1.05])
ax = gca; ax.FontSize = 30; ax.XTick = [1,5,10,50,100];
ax.XTickLabel = [1,5,10,50,100];
set(gca, 'xscale', 'log'); box on;
title({'Community Fitness', 'and Litter Production'})
tx = text(lambda_A, b(1), 'Annual  ');
tx.HorizontalAlignment = 'right';  tx.VerticalAlignment = 'middle';
tx.FontSize = 30;

tx = text(0,1,'(a)', 'Units', 'Normalized'); tx.FontSize = 40;
tx.HorizontalAlignment = 'right'; tx.VerticalAlignment = 'bottom';


axes('Position',[.175 0.35 .175 .25])
scatter(lambda_P, b(2:end), 100,'filled', 'CData', viridis(n));
set(gca, 'xscale', 'log')
title('Perennial Community'); box on;

% Now make the figure of lambda_ratio vs rho.
subplot(1,2,2)
rhovec = linspace(0,1,100);
semilogy(rhovec, rhovec, 'Color', 'black', 'LineWidth', 3);
hold on
semilogy(rhovec, 1./rhovec, 'Color', 'black', 'LineWidth', 3);
hold on
yline(max(lambda_P)/min(lambda_P), ':'); 
yline(min(lambda_P)/max(lambda_P), ':'); xline(rho, '--'); 
xlabel('\rho'); ylabel('{\it\lambda_{P_j}}/{\it\lambda_{P_k}}');
ax = gca; ax.FontSize = 30;
axis([0,1,1/20,20])
ax.YTickLabels = 10.^(-1:1);
tx = text(0,1,'(b)', 'Units', 'Normalized'); tx.FontSize = 40;
tx.HorizontalAlignment = 'right'; tx.VerticalAlignment = 'bottom';


%% Run simulation for the perennials only (i.e., L = 0 NA = 0)

% Set vectors to hold simulation output
NP = zeros(n, gen); NS = zeros(n, gen);
C = zeros(n, gen);
    
% Set the initial conditions for the perennial species.
NP(:,1) = rand(n,1); NS(:,1) = rand(n,1);
C(:,1) = 1 + (NP(:,1)' + gamma*NS(:,1)')*alpha_P*((1-rho)*eye(n) + rho*ones(n));

% Calculate growth rates and new densities every time step
for t = 2:gen
    
    for i = 1:n % Caclulate for each perennial species i
        N = [NP(i,:); NS(i,:)];
        M = [sP*(1-gP) + gP*eP*f*yP(i)/C(i,t-1),...
            yP(i)/C(i,t-1);...
            gP*eP*p1,...
            p2];
        N(:,t) = M*N(:,t-1);
        NP(i,t) = N(1,t);   NS(i,t) = N(2,t);
    end
    
    % Recaculate the intensity of competition in each year
    C(:,t) = 1 + (NP(:,t)' + NS(:,t)'*gamma)*alpha_P*((1-rho)*eye(n) + rho*ones(n));
    
    
end

% Plot figure S3 in the supplement
figure(2)
p = semilogy(NP'); set(p, {'Color'}, num2cell(viridis(n),2));
xlabel('Time'); ylabel('Density');

%% Run Simulations for the full community
% Run for different values of decomposition rate, d
figure(3)

tl = tiledlayout('flow', 'TileSpacing', 'none');

for k = 1:length(d)
    
    % Set vecotrs to hold output.
    NP = zeros(n, gen);     NS = zeros(n, gen);
    C = zeros(n+1, gen);    D = zeros(n+1, gen); % D = 1 + beta*L
    L = zeros(1, gen);      NA = zeros(1, gen);
    Plant_Vec = zeros(n+1,gen); % Total plants in a year
    Comp_Vec = zeros(n+1,gen);
    
    % Set initial conditions
    NP(:,1) = rand(n,1);    NS(:,1) = rand(n,1);
    L(1) = 0;               NA(1) = 0;
    D(:,1) = 1 + beta*L(1);
    Plant_Vec(:,1) = [eA/D(1,1)*gA*NA(1); NP(:,1)];
    Comp_Vec(:,1) = [eA/D(1,1)*gA*NA(1); NP(:,1) + gamma*NS(:,1)];
    C(:,1) = 1 + Comp_Vec(:,1)'*A;
    
    BT = 0;
    
    % Calculate growth rates and population densities each year t
    for t = 2:gen            
        
        % For litter
        L(t) = L(t-1)*(1-d(k)) + b*Plant_Vec(:,t-1) + BT;
        % For the annual
        NA(t) = NA(t-1)*(sA*(1-gA) + gA*eA/D(1,t-1)*yA/C(1,t-1));
        
        % For each perennial species i
        for i = 1:n
            N = [NS(i,:); NP(i,:)];
            M = [sP*(1-gP) + gP*eP*f*yP(i)/(D(i+1,t-1)*C(i+1,t-1)),...
                yP(i)/C(i+1,t-1);...
                p1*gP*eP/D(i+1,t-1),...
                p2];
            N(:,t) = M*N(:,t-1);
            NP(i,t) = N(2,t);   NS(i,t) = N(1,t);
        end
        
        % Add the annual at time step tint
        if t == tint
            NA(t) = 0.1;
        end
        
        % Calculate intensity of litter (D), intensity of competition (C)
        % and the number of seedlings each year.
        D(:,t) = 1 + beta*L(t);
        Plant_Vec(:,t) = [eA/D(1,t)*gA*NA(t); NP(:,t)];
        Comp_Vec(:,t) = [Plant_Vec(1,t); NP(:,t) + gamma*NS(:,t)];
        C(:,t) = 1 + Comp_Vec(:,t)'*A;
        
        
    end
    
    % Determine end community of plants
    Community_end = [NA(gen); NP(:,end)+NS(:,end)];
    % Density of perennials (natives in this scenario)
    Native_Density(k) = sum(Community_end(2:end));
    % Partition out the adults from seeds
    Native_Adults(k) = sum(NP(:,end));
    % Annual seed density
    NA_Density(k) = Community_end(1);
    % Perennial richness
    Native_Diversity(k) = sum(Community_end(2:end) > 0.001);
    % Ask if the invader is present
    Invasive_Diversity(k) = Community_end(1) > 0.001;
    % Take final biomass
    L_Biomass(k) = L(end);
    
    % Plot figure S4 in the supplement
    if any(k == [1:2:length(d)])
        nexttile
        p = semilogy(NP'); set(p, {'Color'}, num2cell(viridis(n),2));
        set(p, {'LineWidth'}, {2});
        hold on
        semilogy(NA, 'Color', 'red', 'LineWidth', 2); 
        semilogy(L, 'Color', 'black', 'LineWidth', 2);
        hold off
        axis([0, 2000, 10^(-3), max(max([Community_end; L(end)]))])
        ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 15;
        title(['d = ',num2str(d(k))])
    end
    
end
xlabel(tl, 'Time', 'FontSize', 50);
ylabel(tl, 'Density', 'FontSize', 50);


% Calcualte the number of adult plants at the final equilibrium. NA is seed
% density, so we have to multiply this by the germination and establishment
% fractions to get growing plant density.
NA_Plants = NA_Density*gA*eA./(1 + betaA*L_Biomass);

% Create figure 7 in the main text.

thresh1 = find(Invasive_Diversity>0,1);
thresh_d1 = mean(d(thresh1-1:thresh1));
thresh2 = find(Native_Diversity==0,1);
thresh_d2 = mean(d(thresh2-1:thresh2));
figure(4)
subplot(2,1,1)
fill([-1,-1,thresh_d1*ones(1,2)], [0,15,15,0], 0.9*ones(1,3),...
    'HandleVisibility', 'off');
hold on
fill([thresh_d2*ones(1,2),2,2], [0,15,15,0], 0.9*ones(1,3),...
    'HandleVisibility', 'off');
b = bar(d, [Native_Diversity;Invasive_Diversity]', 'stacked');
b(1).FaceColor = 'black'; b(1).EdgeColor = zeros(1,3);
b(2).FaceColor = 0.5*ones(1,3); b(2).EdgeColor = zeros(1,3);
hold off
ylabel('Species Richness');
title('Diversity');
set(gca,'FontSize', axFontSize); axis([-0.05,1.05,0,11]);
txt = text(0,1,'a', 'Units', 'Normalized');
txt.FontSize = 50; txt.FontName = 'Helvetica';
txt.HorizontalAlignment = 'right'; txt.VerticalAlignment = 'bottom';

subplot(2,1,2)
Ave_Native_Adults = Native_Adults./Native_Diversity;
Ave_Native_Adults(Native_Diversity == 0) = 0;
fill([-1,-1,thresh_d1*ones(1,2)], [0,1500,1500,0], 0.9*ones(1,3),...
    'HandleVisibility', 'off');
hold on
fill([thresh_d2*ones(1,2),2,2], [0,1500,1500,0], 0.9*ones(1,3),...
    'HandleVisibility', 'off');
b = plot(d, [Ave_Native_Adults; NA_Plants; L_Biomass], '-o');
hold off
b(1).Color = 0.6*ones(1,3); 
b(2).Color = 'black';
b(3).Color = zeros(1,3);
ylabel('Plant Density'); title('Ecosystem State');
set(gca,'FontSize', axFontSize); axis([-0.05,1.05,0,65]);
xlabel('Decomposition Fraction, {\itd}'); 
txt = text(0,1,'b', 'Units', 'Normalized');
txt.FontSize = 50; txt.FontName = 'Helvetica';
txt.HorizontalAlignment = 'right'; txt.VerticalAlignment = 'bottom';


lg = legend('Perennials', 'Annual', 'Litter');
lg.Location = 'southeast'; lg.FontSize = 25;
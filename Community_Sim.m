clear
clc

%% Plotting Parameters
x = viridis(4);     colors = x(2:3,:);
axFontSize = 30;

%% Simulation Parameters

n = 10;        % Number of perennial species in the community
gen = 3000;    % Total time steps of the simulation
tint = 500;   % Initial time steps before annual is introduced
panel_num = 4; % Square-root of the number of panels in the figure
betaA = 1/1.3; % Sensitivity of annual to litter
betaP = 1/15;  % Sensitivity of perennials to litter
beta = [betaA; betaP*ones(n,1)]; 

rho = exp(-2.5); % Niche overlap among the perennials

% Mean fitness of perennials and a measure of variability across species.
lambda_mean = 5; lambda_var = 1;

% Annual Fitness
lambda_A = lambda_mean*25;

% Perennial Fitness
lambda_P = lambda_mean*exp(linspace(-lambda_var,lambda_var,n));

% Distribution of litter production
b = [1, 0.01*linspace(0,1,n)];
d = linspace(0.001,1,panel_num^2);

% Fitness Parameters for the Annual
eA = 1; gA = 1; sA = 0.5;
yA = lambda_A*(1 - sA*(1-gA))/(gA*eA);

% Fitness Parameters for the Perennial
eP = 1; gP = 1; sP = 0.5;
f = 0.1;  s = 0.8;
yP = lambda_P*(1 - sP*(1-gP))/(gP*eP*(f + 1/(1 - s)));

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

subplot(1,2,1)
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
C(:,1) = 1 + NP(:,1)'*alpha_P*((1-rho)*eye(n) + rho*ones(n));

% Calculate growth rates and new densities every time step
for t = 2:gen
    
    for i = 1:n % Caclulate for each perennial species i
        N = [NP(i,:); NS(i,:)];
        M = [sP*(1-gP) + gP*eP*f*yP(i)/C(i,t-1),...
            yP(i)/C(i,t-1);...
            gP*eP,...
            s];
        N(:,t) = M*N(:,t-1);
        NP(i,t) = N(1,t);   NS(i,t) = N(2,t);
    end
    
    % Recaculate the intensity of competition in each year
    C(:,t) = 1 + NP(:,t)'*alpha_P*((1-rho)*eye(n) + rho*ones(n));
    
    
end

% Plot figure S3 in the supplement
figure(2)
p = semilogy(NP' + NS'); set(p, {'Color'}, num2cell(viridis(n),2));
xlabel('Time'); ylabel('Density');

%% Run Simulations for the full community

% Run for different values of decomposition rate, d
for k = 1:length(d)
    
    % Set vecotrs to hold output.
    NP = zeros(n, gen); NS = zeros(n, gen);
    C = zeros(n+1, gen);   D = zeros(n+1, gen); % D = 1 + beta*L
    L = zeros(1, gen);  NA = zeros(1, gen);
    Seedling_Vec = zeros(n+1,gen); % Total seedlings in a year
    
    % Set initial conditions
    NP(:,1) = rand(n,1); NS(:,1) = rand(n,1);
    L(1) = 0;           NA(1) = 0;
    D(:,1) = 1 + beta*L(1);
    Seedling_Vec(:,1) = [eA/D(1,1)*gA*NA(1); NP(:,1)];
    C(:,1) = 1 + Seedling_Vec(:,1)'*A;
    
    % Calculate growth rates and population densities each year t
    for t = 2:gen            
        
        % For litter
        L(t) = L(t-1)*(1-d(k)) + b*Seedling_Vec(:,t-1);
        % For the annual
        NA(t) = NA(t-1)*(sA*(1-gA) + gA*eA/D(1,t-1)*yA/C(1,t-1));
        
        % For each perennial species i
        for i = 1:n
            N = [NP(i,:); NS(i,:)];
            M = [sP*(1-gP) + gP*eP*f*yP(i)/(D(i+1,t-1)*C(i+1,t-1)),...
                yP(i)/C(i+1,t-1);...
                gP*eP/D(i+1,t-1),...
                s];
            N(:,t) = M*N(:,t-1);
            NP(i,t) = N(1,t);   NS(i,t) = N(2,t);
        end
        
        % Add the annual at time step tint
        if t == tint
            NA(t) = 10;
        end
        
        % Calculate intensity of litter (D), intensity of competition (C)
        % and the number of seedlings each year.
        D(:,t) = 1 + beta*L(t);
        Seedling_Vec(:,t) = [eA/D(1,t)*gA*NA(t); NP(:,t)];
        C(:,t) = 1 + Seedling_Vec(:,t)'*A;
        
        
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
    figure(3)
    subplot(panel_num,panel_num,k)
    p = semilogy(NP'); set(p, {'Color'}, num2cell(viridis(n),2));
    hold on
    semilogy(NA, 'Color', 'red'); semilogy(L, 'Color', 'black')
    hold off
    axis([0, 2000, 10^(-3), max(max([Community_end; L(end)]))])
    title(['d = ',num2str(d(k))])
    if k == 1
        tx = text(0,1, 'Density', 'Units', 'Normalized'); 
        tx.FontSize = 50; tx.Rotation = 90;
    end
    if k == length(d)
        tx = text(0,1, 'Time', 'Units', 'Normalized');
        tx.FontSize = 50;
    end
end

% Calculate the number of adult plants at the final equilibrium. NA is seed
% density, so we have to multiply this by the germination and establishment
% fractions to get growing plant density.
NA_Plants = NA_Density*gA*eA./(1 + betaA*L_Biomass);

% Create figure 7 in the main text.
figure(4)
subplot(2,1,1)
b = bar(d, [Native_Diversity;Invasive_Diversity]', 'stacked');
b(1).FaceColor = colors(2,:); b(1).EdgeColor = zeros(1,3);
b(2).FaceColor = colors(1,:); b(2).EdgeColor = zeros(1,3);
ylabel('Species Richness');
title('Community Diversity');
set(gca,'FontSize', axFontSize); axis([-0.05,1.05,0,8]);
txt = text(0,1,'(a)', 'Units', 'Normalized');
txt.FontSize = 40; txt.FontName = 'Times New Roman';
txt.HorizontalAlignment = 'right'; txt.VerticalAlignment = 'bottom';

subplot(2,1,2)
b = bar(d, [Native_Adults; 30*NA_Plants; L_Biomass], 'stacked');
b(1).FaceColor = colors(2,:); b(1).EdgeColor = zeros(1,3);
b(2).FaceColor = colors(1,:); b(2).EdgeColor = zeros(1,3);
b(3).FaceColor = zeros(1,3);     b(3).EdgeColor = zeros(1,3);
ylabel('Plant Density'); title('Ecosystem State');
set(gca,'FontSize', axFontSize); axis([-0.05,1.05,0,850]);
xlabel('Decomposition Fraction, {\itd}'); 
txt = text(0,1,'(b)', 'Units', 'Normalized');
txt.FontSize = 40; txt.FontName = 'Times New Roman';
txt.HorizontalAlignment = 'right'; txt.VerticalAlignment = 'bottom';


lg = legend('Perennials', 'Annual', 'Litter');
lg.Location = 'southeast'; lg.FontSize = 30;
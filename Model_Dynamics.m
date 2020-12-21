clear
clc

% This script plots figure 5 in the main text
% Script requires:
%   "viridis.m" <- a color palette function
%   "APL_BothLitter_Producers_Sim.m" <- function to simulate the
%   annual-perennial dynamics model


%% Parameters

% The coexistence conditions require
%   1. lambdaA and lambdaP,   2. betaA and betaP, and   3. L_P^* and L_A^*

% To find the values of L* for each species as resident, we need
%   1. alphaA and alphaP' = alphaP*(1 + gamma*(1-s))
%   2. The ratios bP/d and bA/d
%   3. betaA and betaP
%   4. lambdaA and lambdaP

% In sum total, we need 4 parameters for each species
%   1. lambda;     2. alpha;
%   3. beta;       4. b/d;

% This is a vector for creating figures in lambdaP vs lambdaA space
lambda = linspace(1,10,100);
% Parameters shared across all simulations
alpha = 1;     bA_d = 50; bP_d = 0;


%% Create figures 5a and 5b

% Consider the two cases where betaA > betaP and betaA < betaP
betaAvec = [1.5,1];      betaPvec = [1,1.5];

% Labels for outcomes in different regions of parameter space
txt_labls = {'Priority Effect','Coexistence',...
    'Annual Excludes Perennial','Perennial Excludes Annual'};


for j = 1:2
    % Set the specific beta values for each case
    betaA = betaAvec(j); betaP = betaPvec(j);
    
    %% Analytical Calculation for equilibrium litter
    % Calculate the equilibrium L*
    LeqA = -(1/betaA + bA_d*(1/alpha))/2 + sqrt(...
        (1/betaA + bA_d*(1/alpha))^2/4 + bA_d*(lambda-1)/(betaA*alpha));
    
    LeqP = -(1/betaP + bP_d*(1/alpha))/2 + sqrt(...
        (1/betaP + bP_d*(1/alpha))^2/4 + bP_d*(lambda-1)/(betaP*alpha));
    
    
    % Calculate the invasion boundaries for each species
    lambdaAinv = lambda.*(1+betaA*LeqP)./(1+betaP*LeqP);
    lambdaPinv = lambda.*(1+betaP*LeqA)./(1+betaA*LeqA);
    
    % Make figures
    labl = {'(a)', '(b)'};
    % Set colors
    figure(1)
    subplot(3,2,j)
    plot(lambdaAinv, lambda, '-', 'Color', 'black', 'LineWidth', 3);
    hold on
    plot(lambda, lambdaPinv, '--', 'Color', 'black', 'LineWidth', 3);
    hold off
    axis([1,max(lambda),1,max(lambda)]);
    ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
    ax.XTick = 1:2:max(lambda); ax.YTick = 1:2:max(lambda);
    xlabel('\it\lambda_A'); ylabel('\it\lambda_P');
    pan_labl = text(0,1,labl(j), 'Units', 'normalized');
    set(pan_labl, 'FontSize', 30); set(pan_labl, 'FontName', 'Times New Roman');
    set(pan_labl, 'VerticalAlignment', 'bottom');
    set(pan_labl, 'HorizontalAlignment', 'right');
    
    % Add in points corresponding to the parameters used in the example
    % dynamics panels (c-f)
    hold on
    if j == 1
        kvec = 2:4;
        scatter(7*ones(1,3), [8,4,5.5], 'filled');
        title('\beta_{\itA} > \beta_{\itP}')
    else
        kvec = [1,3:4];
        scatter(7, 8, 'filled');
        title('\beta_{\itA} < \beta_{\itP}')
    end
    hold off
    
    for k = kvec
        tx = text(0.5,0.5, txt_labls(k),'Units','normalized');
        tx.FontSize = 20; tx.FontName = 'Times New Roman';
    end
    
end

%% Now simulation dynamics in these specific cases

%Plotting Parameters
x = viridis(4);
colors = x(2:3,:);

% General Parameters
gen = [200,700]; d = 0.01;

% Cases
% a. Stable Coexistence
% b. Perennial Excludes Annual
% c. Annual Excludes Perennial
% d. Priority Effects
cases = 4; indx = [3,4,5,6]; %indx here is a plotting parameter
% Figure labels
title_str = {'{\itN_A} excludes {\itN_P}','{\itN_P} excludes {\itN_A}', ...
    'Stable Coexistence', 'Priority Effect'};
subplot_str = {'(c)', '(d)', '(e)', '(f)'};

% Setting the specific parameters to use in the 4 cases
lambdaA = 7*ones(1,4); lambdaP = [4, 8, 5.5, 8];
betaA = [1.5*ones(1,3),1]; betaP = [1*ones(1,3), 1.5];

% Annual Parameters
sA = 0.9;   gA = 0.09;   eA = 1;
bA = 0.5;   alphaA = 1;

% Perennial Parameters
sP = 0.2;   s = 0.2;    f = 0.1;    gP = 0.8;   eP = 0.8;
bP = 0;     alphaP = 1; gamma = 0.1;

% Calculate the stable stage distribution for the perennial
ssd = (1-s)/(gP*eP);

% Initial Conditions
% Each row represents the different panels in teh figure
% row 1 = panel c;          % row 3 = panel e;
% row 2 = panel d;          % row 4 = panel f;
% Each column is a different model component
% Col 1 = NA;   Col 2 = L;  Col 3 = NS;     Col 4 = NP;
init_cond1 = [0,0,ssd/(1 + ssd), 1/(1+ssd);...
              1,0,0,0;...
              0,0,ssd/(1 + ssd), 1/(1+ssd);...
              0,0,ssd/(1 + ssd), 1/(1+ssd)]; 

% Loop over the four cases
for i = 1:cases
    %% Collect parameters into a single cell frame

    % Calculate seed yield values
    yA = lambdaA(i)*(1 - sA*(1-gA))/(gA*eA);
    yP = lambdaP(i)*(1 - sP*(1-gP))/(gP*eP*(f + 1/(1-s)));
    
    % Collect all the parameters into vectors to pass into simulation
    % function
    S = [sA, sP, s];    y = [yA, yP*f, yP];     g = [gA, gP];   e = [eA, eP];
    decay = [bA, bP, d];   alpha = [alphaA, gamma*alphaP, alphaP];
    beta = [betaA(i), betaP(i)];
    
    % Collect all parameters together
    parameters = {S, y, g, e, decay, alpha, beta};
    
    % Delete any previous calculated dynamics
    clear NA L NP
    
    % Run the simulation before the other species is introduced
    sys = APL_BothLitter_Producers_Sim(gen(1), init_cond1(i,:), parameters);
    
    % Specify output back into NA, L, and NP
    NA = sys(1,:); L = sys(2,:);    NP = sum(sys(3:4,:));

    % Specify which species is introduced second in each case
    % If it is case b (i.e., i == 2), the perennial is introduced second.
    % Otherwise the annual is introduced first. Introduction density is
    % exp(-2).
    if any(i == 2)
        NPinit = exp(-1)/2*[1, gP*eP/(1-s)/(1+betaP(i)*L(end))];
        init2 = [sys(1:2,end)', NPinit];
    else
        NAinit = exp(-2);       Linit = 0;
        init2 = [NAinit, Linit, sys(3:4,end)'];
    end
    
    % Run the simulation once the second species is introduced
    sys = APL_BothLitter_Producers_Sim(gen(2)-gen(1), init2, parameters);
    
    % Compile the dynamics since introduction of the second species and the
    % dynamics prior to the introduction of the second species into a
    % single vector of dynamics for each species.
    newNA = sys(1,2:end); newNP = sum(sys(3:4,2:end)); newL = sys(2,2:end);
    NA = [NA, newNA]; L = [L,newL]; NP = [NP, newNP];
        
    % If it is case d (priority effects) run another set of simulations
    % with a second introduction with NA = exp(-0.5) at the second
    % introduction.
    if i == 4
        sys = APL_BothLitter_Producers_Sim(400, [exp(-0.5), sys(2:4,400)'], parameters);
        NA = [NA(1:600), sys(1,2:end)]; L = [L(1:600), sys(2,2:end)];
        NP = [NP(1:600), sum(sys(3:4,2:end))];
    end
    
    % Plot the dynamics as different panels
    figure(1)
    subplot(3,2,indx(i))

    % Plotting specifics for case where perennial excludes annual
    if i == 2
        p(1) = plot(NA); hold on; 
        p(2) = plot((gen(1)+1):gen(2), NP(gen(1):gen(2)-1));
        l = plot(L, '--', 'Color', 'black');
        hold off;
    end
    % Plotting specifics for case where there are priority effects
    if i == 4
        p(1) = plot(gen(1)+1:length(NA), NA(gen(1)+1:length(NA)));
        hold on;
        p(2) = plot(NP);
        l = plot((gen(1)+1):length(L), L(gen(1)+1:length(L)), '--', 'Color', 'black');
        xline(600, ':', 'HandleVisibility', 'off');
        hold off;
    end
    
    % Plotting specifics for case a and c annual excludes perennial and
    % when species stably coexist
    if any(i == [1,3])
        p(1) = plot((gen(1)+1):gen(2), NA(gen(1):gen(2)-1)); hold on; 
        p(2) = plot(NP);
        l = plot((gen(1)+1):gen(2), L(gen(1):gen(2)-1), '--', 'Color', 'black');
        hold off;
    end
    
    % Extra cosmetic features of the figure panels
    set(p, {'Color'}, num2cell(colors,2)); set(p, {'LineWidth'}, {3});
    l.LineWidth = 3;        xline(gen(1),':', 'HandleVisibility', 'off');
    hold off
    xlabel('Year {\itt}'); ylabel('Density'); title(title_str(i));
    ax = gca; ax.FontSize = 20; ax.FontName = 'Times New Roman';
    tx = text(-0.16,1.05, subplot_str(i), 'Units', 'normalized');
    tx.FontName = 'Times New Roman'; tx.FontSize = 30;
    
end

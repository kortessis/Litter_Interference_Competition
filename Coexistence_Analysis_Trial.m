clear
clc

% This script plots the coexistence regions in R0A v R0P space, which is
% figure 3 in the main text. 

% This script requires 
%   1. LitterEq.m
%   2. viridis.m

%% Parameters

% The coexistence conditions require 
%   1. R0A and R0P
%   2. betaA and betaP
%   3. L*

% To plot the coexistence conditions, we need the parameters to determine
% equilibrium litter without a species there.

%   L* itself requires 5 parameters
%       1. beta; 2. alpha; 3. The ratio b/d; 4. R0; 
%       and 5. The ratio b_T/d.

R0 = linspace(1,100,500);
alphaA = 1;  
alphaP = 1; % This is the effective alpha for the perennial. It is 
            % alpha_P*[1 + gamma*(1-p_P)/p_S.

% Now define the ecosystem parameters            
d = 0.9; bT = 0;
% Now define the sensitivity to 
betaA = 1;      betaP = 0;


%% Create a figure showing ratios of bs.

% Made the strength of the differences betweeen bs of species.
b_ratio_vec = [1.5, 5, 10, 100];
% Fix the geometric mean of the b's (the average on the log scale, which is
% appropriate for ratios). 
geo_mean_b = 0.025;

% Define the colors for making the figure
cols = viridis(length(b_ratio_vec)+1);
cols = cols(1:(end-1),:);

% Start the plot of the figure.
figure(1)
subplot(1,2,2)
plot(0,0, 'Color', 'none', 'HandleVisibility', 'off')
axis([1,max(R0),1,max(R0)]);
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
ax.XTick = 0:10:max(R0); ax.YTick = 0:10:max(R0);
xlabel('{\itR}_{0{\itA}}'); ylabel('{\itR}_{0{\itP}}');

% Loop over the different ratios of bA/bP
for i = 1:length(b_ratio_vec)
    % Define bA and bP for each loop
    bA = geo_mean_b*b_ratio_vec(i); bP = geo_mean_b/b_ratio_vec(i);

a;slkdfja ;sdkfjas;dlf 
    
    %% Calculate the equilibrium L* for each species across the range of R0
    
    LAeq = LitterEq(R0, alphaA, betaA, bA/d, bT/d);
    LPeq = LitterEq(R0, alphaP, betaP, bP/d, bT/d);
    
    
    %% Find the invasion boundaries for each value of R0
    
    Ainv = R0.*(1 + betaA*LPeq)./(1 + betaP*LPeq);
    Pinv = R0.*(1 + betaP*LAeq)./(1 + betaA*LAeq);
    
    hold on
    % Plot the invasion boundary for the annual
    plot(Ainv, R0, '-', 'Color', cols(length(b_ratio_vec)+1-i,:),...
        'LineWidth', 2, 'HandleVisibility', 'off');
    % Plot the invasion boundary for the perennial
    plot(R0, Pinv, '-', 'Color', cols(length(b_ratio_vec)+1-i,:), 'LineWidth', 2);
    hold off
end

% Make a legend
l = legend(strsplit(num2str(b_ratio_vec.^2)), 'Location', 'northwest');
l.Box = 'off';
title(l, '{\itb_A/b_P}')

%% Create a figure showing effect of ratios of betas.

% Fix the values of b for each species
bA = 1; bP = 0.01;

% Define the different ratios of betaA/betaP and fix the geometric mean of
% the two (which fixes the mean of lnbetaA and lnbetaP, which is
% appropriate for ratios).
geo_mean_beta = 2; beta_ratio_vec = [1/5, 1/1.2, 1.2, 5];

% Plotting parameters
ltype = {'--', '--', '-', '-'};
lcol = [0,0.5,0.5,0]'*ones(1,3);

% Start the plot
subplot(1,2,1)
plot(0,0,'Color', 'none', 'HandleVisibility', 'off')
axis([1,max(R0),1,max(R0)]);
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
ax.XTick = 0:10:max(R0); ax.YTick = 0:10:max(R0);
xlabel('{\itR}_{0{\itA}}'); ylabel('{\itR}_{0{\itP}}');

% Loop over the different betaA/betaP ratios
for i = 1:length(beta_ratio_vec)
    
    % For each ratio, calculate the specific values of betas.
    betaA = geo_mean_beta*beta_ratio_vec(i); 
    betaP = geo_mean_beta/beta_ratio_vec(i);

    
    %% Calculate the equilibrium L*
    
    LAeq = LitterEq(R0, alphaA, betaA, bA/d, bT/d);
    LPeq = LitterEq(R0, alphaP, betaP, bP/d, bT/d);
        
    %% Find the invasion boundaries
    
    Ainv = R0.*(1 + betaA*LPeq)./(1 + betaP*LPeq);
    Pinv = R0.*(1 + betaP*LAeq)./(1 + betaA*LAeq);
    
    hold on
    % Plot the invasion boundary for the annual
    plot(Ainv, R0, ltype{i}, 'Color', lcol(i,:),...
        'LineWidth', 1, 'HandleVisibility', 'off');
    % Plot the invasion boundary for the perennial
    plot(R0, Pinv, ltype{i}, 'Color', lcol(i,:), 'LineWidth', 1);
    hold off
end

    

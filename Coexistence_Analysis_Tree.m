clear
clc

% This script plots the coexistence regions in lambdaA v lambdaP space


%% Parameters

% The coexistence conditions require 
%   1. lambdaA and lambdaP
%   2. betaA and betaP
%   3. L*

% To plot the coexistence conditions, we need the parameters to determine
% equilibrium litter without a species there.

%   L* itself requires 5 parameters
%       1. beta; 2. alpha; 3. The ratio b/d; 4. lambda; 
%       and 5. The ratio BT/d.

lambda = linspace(1,100,500);
alphaA = 1;  alphaP = 1;

d = 0.9; BT = 0;
betaA = 1;      betaP = 0;


%% Create a figure showing ratios of bs.
b_ratio_vec = [1.5, 5, 10, 100];
geo_mean_b = 0.025;

cols = viridis(length(b_ratio_vec)+1);
cols = cols(1:(end-1),:);

figure(1)
subplot(1,2,2)
plot(0,0, 'Color', 'none', 'HandleVisibility', 'off')
axis([1,max(lambda),1,max(lambda)]);
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
ax.XTick = 0:10:max(lambda); ax.YTick = 0:10:max(lambda);
xlabel('\it\lambda_A'); ylabel('\it\lambda_P');

for i = 1:length(b_ratio_vec)
    bA = geo_mean_b*b_ratio_vec(i); bP = geo_mean_b/b_ratio_vec(i);

    
    %% Calculate the equilibrium L*
    
    if betaA == 0
        LAeq = (bA/d)*(lambda - 1)/alphaA + BT/d;
    else
        if d*alphaA == 0
            LAeq = (lambda - 1)/betaA;
        else
            X1A = (1/alphaA)*bA/d + 1/betaA - BT/d;
            X2A = 1/(betaA)*(bA*(lambda-1)/(alphaA*d) + BT/d);
            LAeq = -0.5*X1A + sqrt(0.25*X1A.^2 + X2A);
        end
    end
    
    if betaP == 0
        LPeq = (bP/d)*(lambda - 1)/alphaP + BT/d;
    else
        if d*alphaP == 0
            LPeq = (lambda - 1)/betaP;
        else
            X1P = (1/alphaP)*bP/d + 1/betaP - BT/d;
            X2P = 1/(betaP)*(bP*(lambda-1)/(alphaP*d) + BT/d);
            LPeq = -0.5*X1P + sqrt(0.25*X1P.^2 + X2P);
        end
    end
    
    
    
    %% Find the invasion boundaries
    
    Ainv = lambda.*(1 + betaA*LPeq)./(1 + betaP*LPeq);
    Pinv = lambda.*(1 + betaP*LAeq)./(1 + betaA*LAeq);
    
    hold on
    % Plot the invasion boundary for the annual
    plot(Ainv, lambda, '-', 'Color', cols(length(b_ratio_vec)+1-i,:),...
        'LineWidth', 2, 'HandleVisibility', 'off');
    % Plot the invasion boundary for the perennial
    plot(lambda, Pinv, '-', 'Color', cols(length(b_ratio_vec)+1-i,:), 'LineWidth', 2);
    hold off
end

l = legend(strsplit(num2str(b_ratio_vec.^2)), 'Location', 'northwest');
l.Box = 'off';
title(l, '{\itb_A/b_P}')

%% Create a figure showing effect of ratios of betas.
bA = 1; bP = 0.01;
geo_mean_beta = 2; beta_ratio_vec = [1/5, 1/1.2, 1.2, 5];

ltype = {'--', '--', '-', '-'};
lcol = [0,0.5,0.5,0]'*ones(1,3);

subplot(1,2,1)
plot(0,0,'Color', 'none', 'HandleVisibility', 'off')
axis([1,max(lambda),1,max(lambda)]);
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
ax.XTick = 0:10:max(lambda); ax.YTick = 0:10:max(lambda);
xlabel('\it\lambda_A'); ylabel('\it\lambda_P');

for i = 1:4
    betaA = geo_mean_beta*beta_ratio_vec(i); 
    betaP = geo_mean_beta/beta_ratio_vec(i);

    
    %% Calculate the equilibrium L*
    
    if betaA == 0
        LAeq = (bA/d)*(lambda - 1)/alphaA + BT/d;
    else
        if d*alphaA == 0
            LAeq = (lambda - 1)/betaA;
        else
            X1A = (1/alphaA)*bA/d + 1/betaA - BT/d;
            X2A = 1/(betaA)*(bA*(lambda-1)/(alphaA*d) + BT/d);
            LAeq = -0.5*X1A + sqrt(0.25*X1A.^2 + X2A);
        end
    end
    
    if betaP == 0
        LPeq = (bP/d)*(lambda - 1)/alphaP + BT/d;
    else
        if d*alphaP == 0
            LPeq = (lambda - 1)/betaP;
        else
            X1P = (1/alphaP)*bP/d + 1/betaP - BT/d;
            X2P = 1/(betaP)*(bP*(lambda-1)/(alphaP*d) + BT/d);
            LPeq = -0.5*X1P + sqrt(0.25*X1P.^2 + X2P);
        end
    end
        
    %% Find the invasion boundaries
    
    Ainv = lambda.*(1 + betaA*LPeq)./(1 + betaP*LPeq);
    Pinv = lambda.*(1 + betaP*LAeq)./(1 + betaA*LAeq);
    
    hold on
    % Plot the invasion boundary for the annual
    plot(Ainv, lambda, ltype{i}, 'Color', lcol(i,:),...
        'LineWidth', 1, 'HandleVisibility', 'off');
    % Plot the invasion boundary for the perennial
    plot(lambda, Pinv, ltype{i}, 'Color', lcol(i,:), 'LineWidth', 1);
    hold off
end

    

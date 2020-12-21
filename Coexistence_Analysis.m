clear
clc

% This script plots the coexistence regions in lambdaA v lambdaP space
% The function requires the function "viridis.m", which is just a color
% palette.

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

% The state space here is lambdaP versus lambdaA, so we specify a range of
% values of lambda to look at.
lambda = linspace(1,10,100);
alpha = 1;  % Set the competition coefficient here that applies to both sps.
            % Changing this value changes L* for both species, according to
            % the expressions given in the supplement. 


%% Subplot of the the figure showing the effect of the trade-off strength
% The trade-off here is characterized by the ratio of betas, betaA/betaP
    % When betaA/betaP > 1, stable coexistence is possible
    % When betaA/betaP < 1, contingent exclusion is possible
    % When betaA/betaP = 1, neutral coexistence is possible

% The following lines evenly compare ratios of beta where betaA > betaP and
% where betaP > betaA. It does this while keeping the average sensitivity
% to litter across species constant, i.e. beta_ bar = (betaA + betaP)/2 is
% constant.
to_st = [1/5,1/1.25,1.25,5]; beta_bar = 5;
bA_d = 10; bP_d = 0.1;
% Make new color plotting parameters
colors = viridis(length(to_st) + 2);
clear p

% Plot figure 6a in the main text
figure(1)
subplot(1,2,1)
plot(lambda, lambda, 'Color', 'none', 'HandleVisibility', 'off')
hold on

for i = 1:length(to_st)
    % Find the values of beta given a specific tradeoff ratio (to_st)
    betaA = 2*beta_bar/(1+to_st(i)); betaP = betaA/to_st(i);
    
    % Calculate the equilibrium L value
    LeqA = -(1/betaA + bA_d*(1/alpha))/2 + sqrt(...
    (1/betaA + bA_d*(1/alpha))^2/4 + bA_d*(lambda-1)/(betaA*alpha));
    LeqP = -(1/betaP + bP_d*(1/alpha))/2 + sqrt(...
    (1/betaP + bP_d*(1/alpha))^2/4 + bP_d*(lambda-1)/(betaP*alpha));
    
    % Calculate the invasion boundary (these are just the two inequalities
    % in the coexistence condition)
    lambdaPinv = lambda.*(1 + betaP*LeqA)./(1 + betaA*LeqA);
    lambdaAinv = lambda.*(1 + betaA*LeqP)./(1 + betaP*LeqP);
    
    % Specify that the invasion boundaries as different line types to 
    % differentiate stable coexistence from priority effects
    if to_st(i) < 1
        linetype = ':';
    else
        linetype = '--';
    end
    
    % Plot the invasion boundaries for the perennial and the annual
    pP(i) = plot(lambda, lambdaPinv, 'LineStyle', linetype);
    pA(i) = plot(lambdaAinv, lambda, 'LineStyle', linetype, 'HandleVisibility', 'off');

end

% A bunch of plotting parameters to make the figure look nicer
hold off
pan_labl(1) = text(0,1,'(a)', 'Units', 'normalized');

set(pP, {'Color'}, num2cell(colors(2:end-1,:),2));
set(pP, {'LineWidth'}, {3})
set(pA, {'Color'}, num2cell(colors(2:end-1,:),2));
set(pA, {'LineWidth'}, {3})
xlabel('\it\lambda_A'); ylabel('\it\lambda_P')
axis([1,max(lambda), 1, max(lambda)])
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
ax.XTick = 1:2:max(lambda); ax.YTick = 1:2:max(lambda);

lg = legend(strsplit(num2str(to_st))); lg.Title.String = '$\beta_A/\beta_P$';
lg.Interpreter = 'Latex';









%% Subplot of the figure showing the effect of the b/d

% Making figure 6b in the main text
clear pA pP
d = 0.5;    bbar = d; % Set d and the geometric mean b across species
bratio = [200,30,2];  % This is the ratio b_A/b_d
betaA = 3; betaP = 1; % Specify litter sensitivities
colors = viridis(length(bratio)+1);


subplot(1,2,2)
% Create axes
plot(lambda, lambda, 'Color', 'none', 'HandleVisibility', 'off')
hold on

% Plot the asymptotic invasion boundaries as b/d -> \ifty
pend1 = plot(lambda, 1 + (lambda - 1)*betaP/betaA, 'Color', 'black');
pend2 = plot(lambda, lambda, 'Color', 'black', 'HandleVisibility', 'off');

% Here we calculate and plot the invasion boundaries for the values of
% bratio specified above
for i = 1:length(bratio)
    
    % Fix the values of bA and bP for each bratio value
    bA_d = exp(log(bbar) + 0.5*log(bratio(i))); 
    bP_d = exp(log(bbar) - 0.5*log(bratio(i)));
    
    % Calculate equilibrium litter values
    LeqA = -(1/betaA + bA_d*(1/alpha))/2 + sqrt(...
        (1/betaA + bA_d*(1/alpha))^2/4 + bA_d*(lambda-1)/(betaA*alpha));
    
    LeqP = -(1/betaP + bP_d*(1/alpha))/2 + sqrt(...
        (1/betaP + bP_d*(1/alpha))^2/4 + bP_d*(lambda-1)/(betaP*alpha));
    
    % Calculate invasion boundaries
    lambdaPinv = lambda.*(1 + betaP*LeqA)./(1 + betaA*LeqA);
    lambdaAinv = lambda.*(1 + betaA*LeqP)./(1 + betaP*LeqP);
    
    % Plot the boundaries
    pP(i) = plot(lambda, lambdaPinv, 'LineStyle', ':');
    pA(i) = plot(lambdaAinv, lambda, 'LineStyle', ':', 'HandleVisibility', 'off');
end

% A bunch of plotting parameters to make the figures look nice.
set([pend1,pend2], {'LineWidth'}, {4});
hold off
set(pA, {'Color'}, num2cell(colors(1:(end-1),:),2));
set(pA, {'LineWidth'}, {3})
set(pP, {'Color'}, num2cell(colors(1:(end-1),:),2));
set(pP, {'LineWidth'}, {3})
xlabel('\it\lambda_A'); ylabel('\it\lambda_P')
axis([1,max(lambda), 1, max(lambda)])
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
ax.XTick = 1:2:max(lambda); ax.YTick = 1:2:max(lambda);

pan_labl(2) = text(0, 1,'(b)', 'Units', 'normalized'); 
set(pan_labl, {'FontSize'}, {30}); set(pan_labl, {'FontName'}, {'Times New Roman'});
set(pan_labl, {'VerticalAlignment'}, {'bottom'});
set(pan_labl, {'HorizontalAlignment'}, {'right'});

lgd = legend(['\infty', strsplit(num2str(bratio))]);
lgd.Title.String = '{\itb_A}/{\itb_P}'; lgd.FontSize = 25;
lgd.FontName = 'Times New Roman'; lgd.Box = 'off';
lgd.Location = 'NorthWest';




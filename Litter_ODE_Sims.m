% This script generates figure 4 of the main text
% The script requires the following functions
%   "viridis.m" <- a color palette function
%   "Litter_Resources_ODE.m" <- function specifying the ODE for the ode solver 

clear
clc

%% Environment Parameters
S = 5;
c = 2;
r = 1;

tsteps = 1000;

%% Species response parameters
Lstar = [0.5,2];
Rstar = [0.5,2];

%% Species effect Parameters
m = 1*[0.1,0.1];
meanv = log(1);
to_value = 2; % Specify the tradeoff value for fig 4b. 
              % This tradeoff value is the ratio of v values, v1/v2. When
              % to_value = 1, there is no tradeoff. Since we assume that
              % L*2 > L*1, species 2 is the better litter competitor and so
              % v1/v2 > 1 implies species 1 effects litter at lot, thereby
              % limiting its own growth rate. This is a stable coexistence
              % configuration. Setting v1/v2 < 1 is more likely to yield
              % priority effects.
a = [3, 3];

%% Solving for extra parameters
w = m./(a.*Rstar);
beta = (S./Rstar - 1)./Lstar;

%% Find the analytical Equilbiria
Req = (beta(2) - beta(1))/(beta(2)/Rstar(1) - beta(1)/Rstar(2));

Leq = (Rstar(2) - Rstar(1))/(beta(1)*Rstar(1) - beta(2)*Rstar(2));

%% Now we plot the size of the coexistence region and the environment vector
% Start off by creating figure 4a of the main text.

% Tradeoffs in vm
to = 10.^(linspace(-2,2,200)); 
% This creates a vector of evenly spaced tradeoff values on the log scale,
% which is appropriate for a ratio.

% Here we just calculate the middle term of the coexistence condition and
% call it the environment vector
env_vec = Leq/(S/Req - 1);

% Specify the values of v for each species given each tradeoff value
% specified above and determine the correspond effect vector for each
% species.
for i = 1:length(to)
    v = exp(meanv + log(to(i))/2*[1,-1]);
    B_vec(i,:) = v.*m.*(1 + beta*Leq)./a;
end

% Plot the effect vectors
figure(1)
subplot(1,2,1)
loglog(to, B_vec(:,1), 'Color', 'black', 'LineWidth', 3)
hold on
loglog(to, B_vec(:,2), '-.', 'Color', 'black', 'LineWidth', 3)
hold off

% Add in labels for the effect vectors and parameter space outcomes
xloc = [0.15, 0.5, 0.5, 0.8];
yloc = [0.5, 0.2, 0.8, 0.5];
txt_String = {['PriorityEffect'], '{\itB}_1 Excludes {\itB}_2',...
    '{\itB}_2 Excludes {\itB}_1',   ['StableCoexistence']};
for i = 1:length(xloc)
    tx = text(xloc(i), yloc(i), txt_String(i), 'Units', 'Normalized');
    tx.FontSize = 20;       tx.FontName = 'Helvetica';
    tx.HorizontalAlignment = 'center';
end

pos_indx = 175;
txB2 = text(to(pos_indx), B_vec(pos_indx,1), '{\itB}_2 Effect Ratio');
txB1 = text(to(pos_indx), B_vec(pos_indx,2), '{\itB}_1 Effect Ratio');
set([txB1,txB2], {'FontSize'}, {20});
set([txB1,txB2], {'FontName'}, {'Helvetica'});
txB1.Rotation = 26;
txB2.Rotation = -26;

% Add axis labels
xlabel('\nu_1/\nu_2'); ylabel('{\itc}/{\itr}'); title('Competitive OUtcomes');
axis([min(to), 0.6*max(to), 0.1*env_vec, max(B_vec,[],'all')]);
ax = gca; ax.FontSize = 30; ax.FontName = 'Helvetica';
ax.YMinorTick = 'off';

% Here we convert the y-axis to be on the scale of c/r. This is justified
% by the fact that the environment vector is constant w.r.t. v1/v2, the
% x-axis. First, we specify the values we want to present and then create
% corrsponding ticks on the y-axis and give them the appropriate labels.
cr_y_axis = 10.^[-1,0,1];
ax.XTick = 10.^(-2:2); ax.YTick = cr_y_axis*env_vec;
ax.XTickLabel = 10.^(-2:2); ax.YTickLabel = cr_y_axis;

SubLabl = text(0,1,'(a)', 'Units', 'Normalized');
SubLabl.FontSize = 40; SubLabl.FontName = 'Helvetica';
SubLabl.HorizontalAlignment = 'right'; SubLabl.VerticalAlignment = 'bottom'; 


%% Gradient Figure
% This creates figure 4b
tsteps = 10000; % Specify that simulations run for this length of time.

% Caclualte values of v and w.
v = exp(meanv + log(to_value)/2*[1,-1]);
w = m./(a.*Rstar);

% Specify gradient of decomposition rate
c = exp(linspace(log(0.01), log(10), 100));
r = 1;

% Loop over decomposition rates
for i = 1:length(c)
    % Set initial conditions [B1(0), B2(0), R(0), L(0)];
    % Plants have initial density 0.01, R(0) = S, and litter is absent
    init_cond = [0.01,0.01,S,0];
        
    % Solve the ODE in Litter_Resources_ODE function
    [T,X] = ode45(@(t,X) Litter_Resources_ODE(t,X,w,a,beta,m,v,S,r,c(i)),...
        [0 tsteps], init_cond);
    i
    % Take the final value to be the equilibrium
    EQ(i,:) = X(end,:);
end

% Plot the output
colors = viridis(4); 
subplot(1,2,2)
p = semilogx(c, EQ);
p(1).Color = colors(2,:); p(2).Color = colors(3,:); 
set(p, {'LineWidth'}, {3}); set(p(3:4), {'Color'}, {'black'});
p(4).LineStyle = '--';
xlabel('Decomposition Rate, {\itc}'); ylabel('Equilibrium Densities');
title('Ecosystem State')

ax = gca; ax.FontSize = 30; ax.FontName = 'Helvetica';
axis([min(c),max(c),0,S])
ax.XTick = 10.^[-3:2]; ax.XTickLabels = 10.^[-3:2];

SubLabl = text(0,1,'(b)', 'Units', 'Normalized');
SubLabl.FontSize = 40; SubLabl.FontName = 'Helvetica';
SubLabl.HorizontalAlignment = 'right'; SubLabl.VerticalAlignment = 'bottom'; 


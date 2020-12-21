% Script to make figure 2 in the main text.
% Script requires the functions
%   "viridis.m" <- a color palette function
%   "Litter_ODE.m" <- a function specifying the ODE for the ode solver

% Specify colors
dummy_color = viridis(5);
colors = [dummy_color(4, :); dummy_color(3,:); dummy_color(2,:)];
 
%% Create Figure 2, panel a
% Negative Litter effects

% Species parameters
beta = [0.8, 2, 5]; 
m = 0.2*ones(1,3);
r = [0.35, 0.5, 0.7];
waS = r + m;

% Calculate L* for each species and create labels for plotting
Lstar = (waS./m - 1)./beta;
Ltext = {'$L^*_1$', '$L^*_2$', '$L^*_3$'}; 

% Create a vector of litter values for plotting
L = 0:0.01:1.2*max(Lstar);

% Make figure 2a
figure(1)
subplot(2,2,1)

% Plot each species growth rates
p = plot(L, waS'*ones(size(L))./(1 + beta'*L) - m', 'LineWidth', 3);
hold on
% Specify the Lstars on thh figure
scatter(Lstar, zeros(size(beta)), 150, 'filled', 'CData', colors);
yline(0, ':', 'Color', 'black', 'LineWidth', 2);
% Add lables
for i = 1:length(beta)
    plot(Lstar(i)*ones(1,2), [-10*m(i),0], '--', 'Color', colors(i,:));
    txt = text(Lstar(i), -0.3, Ltext(i));
    txt.Interpreter = 'Latex'; txt.FontSize = 30;
    txt.Color = colors(i,:);
end
hold off
% Extra plotting parameters
xlabel('Litter Density');
ylabel('Per-unit biomass growth rate');
title('Plant Growth Rates')
ax = gca; ax.XTick = 0; ax.YTick = 0;
ax.FontSize = 25; ax.FontName = 'Helvetica';
set(p, {'Color'}, num2cell(colors, 2));
axis([0, max(L), -0.4, 1.05*max(r)])
labl_text = text(0,1,'(a)', 'Units', 'Normalized');
labl_text.FontName = 'Helvetica';
labl_text.FontSize = 30; labl_text.HorizontalAlignment = 'right';
labl_text.VerticalAlignment = 'bottom'; 

%% Create Figure 2, panel b
% Dynamics from model 

% Specify the values of decomposition rate and litter conversion
c = 0.04; v = 0.1*[1:3];

% Solve the ODE with initial conditions
%   L(0) = 0 and B_i(0) = 0.25 for all species
% Solve the ODE to time t=1000;
[T,X] = ode45(@(t,X) Litter_ODE(t, X, waS, beta, v, m, c),...
    [0 1000], 0.25*[ones(1,length(beta)),0]);

% Plot the outputs of species densities and litter densities
subplot(2,2,2)
Bplots = plot(T, X(:,1:3)); hold on; Lplot = plot(T, X(:,end));
yline(max(Lstar),':', 'Color', 0.25*ones(1,3)); hold off;
set(Bplots, {'Color'}, num2cell(colors, 2)); set(Bplots, {'LineWidth'}, {3});
Lplot.Color = 'black'; Lplot.LineWidth = 3; Lplot.LineStyle = '--';
xlabel('Time'); ylabel('Densities'); title('Dynamics');
axis([0,150,0,max(X,[],'all')])
ax = gca; ax.FontSize = 25; ax.FontName = 'Helvetica';

labl_text = text(0,1,'(b)', 'Units', 'Normalized');
labl_text.FontName = 'Helvetica';
labl_text.FontSize = 30; labl_text.HorizontalAlignment = 'right';
labl_text.VerticalAlignment = 'bottom'; 

% Here we are adding in labels on the plots
x = [76,76,40,76]; y = [0.6, 0.125, 0.025, 0.85];
txt_vec = {'$B_1$', '$B_2$', '$B_3$', '$L$'};
for i = 1:4
    ln_lab = text(x(i), y(i), txt_vec{i});
    ln_lab.FontSize = 26; ln_lab.Interpreter = 'Latex';
    if any(i == 1:3)
        ln_lab.Color = colors(i,:);
    else
        ln_lab.Color = 'black';
    end
end

%% Create figure 2, panel c
% Effects of decomposition on equilibirum

subplot(2,2,3)
% Set parameters
cvec = 0.01:0.01:1; vm = 0.5; Lstar = 1;
% Plot equilbirium outcomes
Bplots = plot(cvec, Lstar*cvec/vm); 
hold on; Lplot = yline(Lstar); hold off;
Bplots.Color = colors(1,:); set(Bplots, 'LineWidth', 3);
Lplot.Color = 'black'; Lplot.LineWidth = 3; Lplot.LineStyle = '--';
xlabel('Decomposition rate, {\itc}'); ylabel('Equilibrium Densities'); 
title('Equilibrium Dynamics');
ax = gca; ax.FontSize = 25; ax.FontName = 'Helvetica';

% Add lables
labl_text = text(0,1,'(c)', 'Units', 'Normalized');
labl_text.FontName = 'Helvetica';
labl_text.FontSize = 30; labl_text.HorizontalAlignment = 'right';
labl_text.VerticalAlignment = 'bottom'; 


y = Lstar*[1, cvec(end)/vm]; text_vec = {'$\hat{L}$', '$\hat{B}$'};

for i = 1:2
    ln_lab = text(cvec(end), y(i), text_vec{i});
    set(ln_lab, 'FontSize', 26); set(ln_lab, 'Interpreter', 'latex');
    ln_lab.HorizontalAlignment = 'left';
    ln_lab.VerticalAlignment = 'middle';
end

%% Create figure 2, panel d
% Effects of litter conversion on equilibrium
subplot(2,2,4)
% Set parameters
vmvec = 0.25:0.01:2; c = 0.5; Lstar = 1;
% Plot equilibrium outcomes
Bplots = plot(vmvec, Lstar*c./vmvec); 
hold on; Lplot = yline(Lstar); hold off;
Bplots.Color = colors(1,:); set(Bplots, 'LineWidth', 3);
Lplot.Color = 'black'; Lplot.LineWidth = 3; Lplot.LineStyle = '--';
xlabel('Litter production rate, {\it\nu m}'); ylabel('Equilibrium Densities'); 
title('Equilibrium Dynamics');
ax = gca; ax.FontSize = 25; ax.FontName = 'Helvetica';

% Add labels
labl_text = text(0,1,'(d)', 'Units', 'Normalized');
labl_text.FontName = 'Helvetica';
labl_text.FontSize = 30; labl_text.HorizontalAlignment = 'right';
labl_text.VerticalAlignment = 'bottom'; 

y = Lstar*[1, c/vmvec(end)]; text_vec = {'$\hat{L}$', '$\hat{B}$'};

for i = 1:2
    ln_lab = text(vmvec(end), y(i), text_vec{i});
    set(ln_lab, 'FontSize', 26); set(ln_lab, 'Interpreter', 'latex');
    ln_lab.HorizontalAlignment = 'left';
    ln_lab.VerticalAlignment = 'middle';
end
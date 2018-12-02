%% Wave Pinning
%% entire phase space, Size Limit WP, x = 6.429 ... 6.857, y 1.035, 1.056
clear all;
si = [linspace(6, 8, 15), linspace(8.1, 21, 35)];
conc_tot = [linspace(0.95, 1.03, 14), linspace(1.035, 1.08, 22), linspace(1.095, 1.2, 14)]; % concentrations
Da = 0.1;
workspace_name = 'wp_10000.mat';
%% only CPSS for L = 2fold
clear all;
si = 2*linspace(6.429, 6.857, 5);
conc_tot = linspace(1.035, 1.056, 5);
Da = 0.4;
workspace_name = 'wp_10000_2fold.mat';
%% only CPSS for L = 4fold
clear all;
si = 4*linspace(6.429, 6.857, 5);
conc_tot = linspace(1.035, 1.056, 5);
Da = 1.6;
workspace_name = 'wp_10000_4fold.mat';
%% only CPSS for L = 0.5fold
clear all;
si = 0.5*linspace(6.429, 6.857, 5);
conc_tot = linspace(1.035, 1.056, 5);
Da = 0.025;
workspace_name = 'wp_10000_0_5fold.mat';
%%
tic
delta = 1/9;%0.005;
gamma = 1/9;
Dc = 100000;
K = 1*ones(1, length(conc_tot));%*conc_tot;
k0 = 0.067/9;

u1_all = cell(1, length(si));
u2_all = cell(1, length(si));

for j = 1:length(si)
    u1 = zeros(length(conc_tot), 200);
    u2 = zeros(length(conc_tot), 200);
    parfor i = 1:length(conc_tot)
        [u, v] = wave_pin(Da, Dc, si(j), delta, gamma, K(i), conc_tot(i), k0);
        u1(i, :) = u(end, :);
        u2(i, :) = v(end, :);
    end
    u1_all{j} = u1;
    u2_all{j} = u2;
end
toc

logical = cell2mat(cellfun(@(x) (max(x') - min(x'))./max(x')<0.05, u1_all, 'uni', 0)');

% slope (a.u./microns):
% max(abs(diff(u(end, :))/(max(u(end, :) - min(u(end, :))))))*200/30
% 0.1519
[~, top]=min(fliplr(logical),[],2);
size_logical = size(logical);
top = size_logical(2) - top + 1;
[~, bot]=min(logical,[],2);
save(workspace_name);
%% Otsuji entire phase space
clear all;
si = [linspace(5, 5.3, 10), linspace(5.4, 21, 40)]; % cell sizes
conc_tot = logspace(-1, 1.7, 50); % concentrations
workspace_name = 'ot_10000s.mat';
D1 = 0.1;
%% Otsuji 2 fold
clear all;
si = 2*linspace(5.1, 5.2, 5); % cell sizes
conc_tot = linspace(1.436, 2.385, 5); % concentrations
workspace_name = 'ot_10000s_2fold.mat';
D1 = 0.4;
%% Otsuji 4 fold
clear all;
si = 4*linspace(5.1, 5.2, 5); % cell sizes
conc_tot = linspace(1.436, 2.385, 5); % concentrations
workspace_name = 'ot_10000s_4fold.mat';
D1 = 1.6;
%% Otsuji 0.5 fold
clear all;
si = 0.5*linspace(5.1, 5.2, 5); % cell sizes
conc_tot = linspace(1.436, 2.385, 5); % concentrations
workspace_name = 'ot_10000s_0_5fold.mat';
D1 = 0.025;
%%
a1 = 1;
a2 = 0.7*ones(1, length(conc_tot));%*conc_tot; % a2 needs to be scaled by concentration for right antagonism
s = 1;
D2 = 100000;
tic
for j = 1:length(si)
    u1 = zeros(length(conc_tot), 200);
    u2 = zeros(length(conc_tot), 200);
    parfor i = 1:length(conc_tot)
        [u, v] = otsuji(si(j), D1, D2, a1, a2(i), s, conc_tot(i));
        u1(i, :) = u(end, :);
        u2(i, :) = v(end, :);
    end
    u1_all{j} = u1;
    u2_all{j} = u2;
end
toc

logical = cell2mat(cellfun(@(x) (max(x') - min(x'))./max(x')<0.05, u1_all, 'uni', 0)');

[~, top]=min(fliplr(logical),[],2);
size_logical = size(logical);
top = size_logical(2) - top + 1;
[~, bot]=min(logical,[],2);
save(workspace_name);
%% Top: get smoothed outline by averaging all sizes with same concentration
u = unique(top);
u = u(u<length(top));
sizes_smooth_top = zeros(1, length(u));
conc_smooth_top = zeros(1, length(u));

for i = 1:length(u)
    sizes_smooth_top(i) = mean(si(top==u(i)));
    conc_smooth_top(i) = mean([conc_tot(u(i)),conc_tot(u(i)+1)]);
end
%% Bottom: same as above
u = unique(bot);
u = u(u>1);
sizes_smooth_bot = zeros(1, length(u));
conc_smooth_bot = zeros(1, length(u));

for i = 1:length(u)
    sizes_smooth_bot(i) = mean(si(bot==u(i)));
    conc_smooth_bot(i) = mean([conc_tot(u(i)),conc_tot(u(i)-1)]);
end
%% Otsuji files for Python
csvwrite('otsuji_matlab_sizeDosageTop_10000.csv', [sizes_smooth_top', conc_smooth_top']);
csvwrite('otsuji_matlab_sizeDosageBot_10000.csv', [sizes_smooth_bot', conc_smooth_bot']);
%% WP files for Python
csvwrite('WP_matlab_sizeDosageTop.csv', [sizes_smooth_top', conc_smooth_top']);
csvwrite('WP_matlab_sizeDosageBot.csv', [sizes_smooth_bot', conc_smooth_bot']);
                                                
%% Goryachev Nate entire phase space
clear all;
si = linspace(5, 20, 20); % cell sizes, DO NOT USE LOG SPACING, OTHERWISE
                            % boundary calculation below (average) doesn't
                            % work anymore!
conc_tot = linspace(0.1, 50, 20); % concentrations
D1 = 0.1;
workspace_name = 'goryachev_10000.mat';
%% Goryachev L 2 fold
clear all;
si = 2*linspace(8.9, 11.3, 5); % cell sizes, DO NOT USE LOG SPACING, OTHERWISE
                            % boundary calculation below (average) doesn't
                            % work anymore!
conc_tot = linspace(0.1, 8, 8); % concentrations
D1 = 0.4;
workspace_name = 'goryachev_10000_2fold.mat';
%% Goryachev L 4 fold
clear all;
si = 4*linspace(8.9, 11.3, 5); % cell sizes, DO NOT USE LOG SPACING, OTHERWISE
                            % boundary calculation below (average) doesn't
                            % work anymore!
conc_tot = linspace(0.1, 8, 8); % concentrations
D1 = 1.6;
workspace_name = 'goryachev_10000_4fold.mat';
%% Goryachev L 0.5 fold
clear all;
si = 0.5*linspace(8.9, 11.3, 5); % cell sizes, DO NOT USE LOG SPACING, OTHERWISE
                            % boundary calculation below (average) doesn't
                            % work anymore!
conc_tot = linspace(0.1, 8, 8); % concentrations
D1 = 0.025;
workspace_name = 'goryachev_10000_0_5fold.mat';
%%
a1 = 0.67/101;
a2 = 0.33/101;
a3 = 0.01;
D2 = 10000000;

u1_all = cell(1, 1);
u2_all = cell(1, 1);

tic
for j = 1:length(si)
    u1 = zeros(length(conc_tot), 200);
    u2 = zeros(length(conc_tot), 200);
    parfor i = 1:length(conc_tot)
        [x, u, v, mass] = goryachev(si(j), D1, D2, a1, a2, a3, conc_tot(i));
        u1(i, :) = u(end, :);
        u2(i, :) = v(end, :);
    end
    u1_all{j} = u1;
    u2_all{j} = u2;
end
toc
logical = cell2mat(cellfun(@(x, y) ((max(x') - min(x'))./...
                   (sum(x')+sum(y'))*200) < 0.01, u1_all, u2_all, 'uni', 0)');

[~, top]=min(fliplr(logical),[],2);
size_logical = size(logical);
top = size_logical(2) - top + 1;
[~, bot]=min(logical,[],2);
% csvwrite('goryachev_matlab_sizeDosageTopBot.csv', [conc_tot', si', top, bot]);
save(workspace_name);
%% Plot all simulations with one figure per size
for j = 1 :  length(u1_all)
    figure; hold on;
    for i = 1:n
    % plot(u*ones(size(u1(i, :))));
%     plot(u1_all{j}(i, :)/max(u1_all{j}(i, :))); % Normalized to max
    plot(u1_all{j}(i, :), '--b');
%     plot(u2_all{j}(i, :), '-r');
    % rho_tot = sum(u1(i, :))+sum(u2(i, :));
    end
end
%% Plot entire phase space
% n=30; n_s = 30;
sizes = repmat(si', 1, n);
concs = repmat(conc_tot', 1, n_s)';
% logical = cell2mat(cellfun(@(x) (max(x') - min(x'))./max(x')<0.05, u1_all, 'uni', 0)');

figure; hold on;
plot(sizes(logical==0), concs(logical==0), '.', 'MarkerSize', 15);
plot(sizes(logical==1), concs(logical==1), '.', 'MarkerSize', 15);
ax = gca;
ax.FontSize=18;
xlabel('system size (\mum)');
ylabel('concentration (a.u.)');
%% Calculate CPSS vs D
% Breakdown for Otsuji happens between 1.0 and 1.25, verify this for
% different Ds and sizes. This can be done using Lcpss =? A*sqrt(D), with 
% A = 1.125/sqrt(0.1) = 3.5576 and turns out to be true if cytoplasmic 
% diffusion is infinite. Otherwise breakdown deviates from this equation 
% for obvious reasons!
Ds = [0.9, 3.6, 14.4, 57.6, 230.4];
Ls = 3.5576*sqrt(Ds);

for i = 1:length(Ds)
     L1 = Ls(i)-0.1*Ls(i);
     L2 = Ls(i)+0.1*Ls(i);
     u_below = otsuji(L1, Ds(i), 10000000, 25, a2, s);
     u_above = otsuji(L2, Ds(i), 10000000, 25, a2, s);
     close all;
     figure; hold on;
     plot(u_below(end, :));
     plot(u_above(end, :));
     pause();
end
f = fit(Ds(1:end)', Ls(1:end)', 'a*sqrt(x)');
figure; hold on;
plot(Ds, f(Ds), 'o');
plot(Ds, Ls);
csvwrite('fit_D.csv', f(Ds));
csvwrite('size.csv', Ls');
csvwrite('Ds.csv', Ds');

%% Get graphs for different diffusion rates
Ds = [0.1, 0.9, 3.6];
figure; hold on;
for i = 1:length(Ds)
     u_below = otsuji(10, Ds(i), 10000000, 25, a2, s);
     plot(u_below(end, :));
end
%%
% u_python = zeros(200, 3);
figure; hold on;
for i = 1:length(Ds)
     u = otsuji(10, Ds(i), 10000000, 25, a2, s);
     u_python(:, i) = u(end, :);
     plot(u_python(:, i));
end
csvwrite('/Users/hubatsch/Desktop/PARmodelling/MatlabDet/ActivatorInhibitor/ints_otsuji_off.csv', u_python);
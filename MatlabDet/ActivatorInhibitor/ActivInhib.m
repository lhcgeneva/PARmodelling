%% Otsuji
n_s = 25; % Precision of sampling of cell sizes
si = linspace(5, 20, n_s); % cell sizes
n = 25; % Precision of sampling of concentrations
conc_tot = logspace(-0.3, 2, n); % concentrations
a1 = 0.5;
a2 = 0.7;
s = 1;
D1 = 0.1;
D2 = 1000000;

for j = 1:length(si)
    u1 = zeros(n, 200);
    u2 = zeros(n, 200);
    parfor i = 1:n
    [u, v] = otsuji(si(j), D1, D2, a1, a2, s, conc_tot(i));
    u1(i, :) = u(end, :);
    u2(i, :) = v(end, :);
    end
    u1_all{j} = u1;
    u2_all{j} = u2;
end

%% Plot all simulations with one figure per size
figure; hold on;
for j = 1 : length(u1_all)
    figure; hold on;
    for i = 1:n
    % plot(u*ones(size(u1(i, :))));
%     plot(u1_all{j}(i, :)/max(u1_all{j}(i, :))); % Normalized to max
    plot(u1_all{j}(i, :));
    % plot(u2(i, :));
    % rho_tot = sum(u1(i, :))+sum(u2(i, :));
    end
end
%% Plot entire phase space
sizes = repmat(si', 1, n);
concs = repmat(conc_tot', 1, n_s)';
logical = cell2mat(cellfun(@(x) (max(x') - min(x')) < 0.01, u1_all, 'uni', 0)');

figure; hold on;
plot(sizes(logical==0), concs(logical==0), '.', 'MarkerSize', 15);
plot(sizes(logical==1), concs(logical==1), '.', 'MarkerSize', 15);
ax = gca;
ax.FontSize=18;
xlabel('system size (\mum)');
ylabel('concentration (a.u.)');
%% Goryachev Nate
n_s = 25; % Precision of sampling of cell sizes
si = linspace(10, 40, n_s); % cell sizes
n = 25; % Precision of sampling of concentrations
conc_tot = logspace(-0.3, 2, n); % concentrations
a1 = 0.1;
a2 = 0.1;
a3 = 0.005;
D1 = 0.1;
D2 = 10000;
u1_all = cell(1, 1);
u2_all = cell(1, 1);
% m = zeros(100001, n);

for j = 1:length(si)
    u1 = zeros(n, 200);
    u2 = zeros(n, 200);
    parfor i = 1:n
    [x, u, v, mass] = goryachev(si(j), D1, D2, a1, a2, a3, conc_tot(i));
    u1(i, :) = u(end, :);
    u2(i, :) = v(end, :);
%     m(:, i) = mass;
    end
    u1_all{j} = u1;
    u2_all{j} = u2;
end
%% Plot entire phase space
sizes = repmat(si', 1, n);
concs = repmat(conc_tot', 1, n_s)';
% logical = cell2mat(cellfun(@(x) (max(x') - min(x')) < 10, u1_all, 'uni', 0)');
logical = cell2mat(cellfun(@(x, y) ((max(x') - min(x'))./(sum(x')+sum(y'))*200) < 0.01, u1_all, u2_all, 'uni', 0)');

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
% % Activator Inhibitor model
% Da = linspace(0.1,0.6, 6);
% % Da = [0.001, Da];
% Dh = 15;
% mu_a = 0.001;
% mu_h = 0.001;
% rho = 0.06;
% rho_a = 0.06;
% rho_h = 0.0006;
% si = 40;
% 
% figure; hold on;
% for i = 1:10
% u1 = activ_inhib(si, Da(i), Dh, mu_a, mu_h, rho, rho_a, rho_h);
% plot(u1(end, :));
% end
%% Otsuji
% si = 1.295:0.004:1.39;
si = 1.0:0.25:5;
% si = 100;
u1 = zeros(length(si), 200);
u2 = zeros(length(si), 200);
a2 = 0.7;
s = 1;
for i=1:length(si)
[u, v] = otsuji(si(i), 0.1, 100, 25, a2, s);
u1(i, :) = u(end, :);
u2(i, :) = v(end, :);
end
%%
figure; hold on;
for i = 1:length(si)
plot(u*ones(size(u1(i, :))));
plot(u1(i, :));
plot(u2(i, :));
% rho_tot = sum(u1(i, :))+sum(u2(i, :));
end
%%
rho_tot = sum(u1(i, :))+sum(u2(i, :));
rho_tot = rho_tot/200;
v = rho_tot/(a2*s*rho_tot+1)^2
u = (rho_tot - v)
plot((u1(2:end, 1)-u)./(u1(1:end-1, 1)-u));
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
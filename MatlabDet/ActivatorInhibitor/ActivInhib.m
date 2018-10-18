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
%%

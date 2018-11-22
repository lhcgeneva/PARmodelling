function [u1, u2] = wave_pin(Da, Dc, si, delta, gamma, K, conc_tot)
% Mass conserved model for wave pinning

x = linspace(0, si, 200);
t = 0:1:100000;
ic = @(x) wp_ic(x, si, conc_tot);
pde = @(x, t, u, DuDx) wp_pde(x, t, u, DuDx, Da, Dc, delta, gamma, K);
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
sol = pdepe(0, pde, ic, @wp_bc, x, t, options);

u1 = sol(:,:,1);
u2 = sol(:,:,2);
% u1 = sol(end,:,1);
% u2 = sol(end,:,2);
% figure;
% plot(x, u1(end, :));
% axis([0, si, 0, 1.5]);

% PDE
function [c,f,s] = wp_pde(x,t,u,DuDx,Da, Dc, delta, gamma, K)

c = [1; 1]; 
f = [Da; Dc] .* DuDx; 
F = u(2)*(0.067 + gamma*u(1)^2/(K^2 + u(1)^2))-delta*u(1);
s = [F; -F]; 

% Initial conditions
function u0 = wp_ic(x, si, conc_tot)
if x < si/2
    u0 = conc_tot*[2*0.2683312; 2.0]; 
else
    u0 = conc_tot*[0; 2.0];
end
% u0 = sigmf(x, 

% Boundary conditions
function [pl,ql,pr,qr] = wp_bc(xl,ul,xr,ur,t)
pl = [0; 0]; 
ql = [1; 1]; 
pr = [0; 0]; 
qr = [1; 1]; 
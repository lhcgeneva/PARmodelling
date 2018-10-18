function u1 = activ_inhib(si, Da, Dh, mu_a, mu_h, rho, rho_a, rho_h)
% This hasn't been tested and might not work correctly (i.e. as intended
% here: http://www.scholarpedia.org/article/Gierer-Meinhardt_model
x = linspace(0, si, 500);
t = 0:1:10000;
ic = @(x) ai_ic(x, si);
pde = @(x, t, u, DuDx) ai_pde(x, t, u, DuDx, Da, Dh, mu_a, mu_h, rho,...
                              rho_a, rho_h);

sol = pdepe(0, pde, ic, @ai_bc, x, t);

u1 = sol(:,:,1);
u2 = sol(:,:,2);
% figure;
% plot(x, u1(end, :));
% axis([0, si, 0, 1.5]);

% PDE
function [c,f,s] = ai_pde(x, t, u, DuDx, Da, Dh, mu_a, mu_h, rho,...
                          rho_a, rho_h)
% Da = 0.1;
Db = 1000000;
% delta = 1;
% gamma = 1;

c = [1; 1]; 
f = [Da; Dh] .* DuDx; 
Ra = rho*u(1)^2/u(2) - mu_a*u(1)+rho_a;
Rh = rho*u(1)^2 - mu_h*u(2)+rho_h;
s = [Ra; Rh]; 

% Initial conditions
function u0 = ai_ic(x, si)
if x < si/2
    u0 = [2*0.2683312; 2.0]; 
else
    u0 = [0; 2.0];
end

% Boundary conditions
function [pl,ql,pr,qr] = ai_bc(xl,ul,xr,ur,t)
pl = [0; 0]; 
ql = [1; 1]; 
pr = [0; 0]; 
qr = [1; 1]; 
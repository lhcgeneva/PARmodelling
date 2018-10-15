function u1 = wave_pin(Da, si, delta)
m = 0;
x = linspace(0, si, 500);
t = 0:1:2000;
ic = @(x) wp_ic(x, si);
pde = @(x, t, u, DuDx) wp_pde(x, t, u, DuDx, Da, delta);
tic
sol = pdepe(m, pde, ic, @wp_bc, x, t);
toc 
u1 = sol(:,:,1);
u2 = sol(:,:,2);
% figure;
% plot(x, u1(end, :));
% axis([0, si, 0, 1.5]);

% PDE
function [c,f,s] = wp_pde(x,t,u,DuDx,Da, delta)
% Da = 0.1;
Db = 1000000;
% delta = 1;
% gamma = 1;
gamma = delta;
k0 = 0.067*delta;
c = [1; 1]; 
f = [Da; Db] .* DuDx; 
F = u(2)*(k0 + gamma*u(1)^2/(1^2 + u(1)^2))-delta*u(1);
s = [F; -F]; 

% Initial conditions
function u0 = wp_ic(x, si)
if x < si/2
    u0 = [2*0.2683312; 2.0]; 
else
    u0 = [0; 2.0];
end

% Boundary conditions
function [pl,ql,pr,qr] = wp_bc(xl,ul,xr,ur,t)
pl = [0; 0]; 
ql = [1; 1]; 
pr = [0; 0]; 
qr = [1; 1]; 
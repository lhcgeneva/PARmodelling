function [x, u1, u2, mass] = goryachev(si, D1, D2, a1, a2, a3, conc_tot)
% Mass conserved model from goryachev FEBS letters

x = linspace(0, si, 200);
t = 0:1:10000;
ic = @(x) ai_ic(x, si, conc_tot);
pde = @(x, t, u, DuDx) ai_pde(x, t, u, DuDx, D1, D2, a1, a2, a3);
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
sol = pdepe(0, pde, ic, @ai_bc, x, t,options);

u1 = sol(:,:,1);
u2 = sol(:,:,2);
mass = sum(u1, 2)+sum(u2, 2);



% PDE
function [c,f,s] = ai_pde(x, t, u, DuDx, D1, D2, a1, a2, a3)

c = [1; 1]; 
f = [D1; D2] .* DuDx; 
F = a1*(u(1)^2)*u(2)+a2*u(1)*u(2)-a3*u(1);
s = [F; -F]; 

% Initial conditions
function u0 = ai_ic(x, si, conc_tot)
if x < si/2
    u0 = conc_tot*[1; 1]; 
else
    u0 = conc_tot*[0; 1];
end

% Boundary conditions
function [pl,ql,pr,qr] = ai_bc(xl,ul,xr,ur,t)
pl = [0; 0]; 
ql = [1; 1]; 
pr = [0; 0]; 
qr = [1; 1]; 




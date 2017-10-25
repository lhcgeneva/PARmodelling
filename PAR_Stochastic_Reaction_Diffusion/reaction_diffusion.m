close all;
% Calculate surface to volume ratio for prolate spheroid
a = 15; % in mu
b = 27; % in mu
e = sqrt(1 - a^2 / b^2);
S = 2 * 3.1415 * a^2 * (1 + b/(a * e) * asin(e)); % in mu^2
V = 4/3 * 3.1415 * a^2 * b;                       % in mu^2
Psi = S/V;

T = 2000;          % Simulation time in s
t_plot = 10;     % Tolerance time till next plot in s
K = 50;        % number of compartments
L = 67.3;       % in mu
h = L/K;        % Compartment length in mu
Si = h^2;       % Compartment area in mu^2
Vi = Si/Psi;    % Cytoplasmic volume associated with each compartment based on Psi in mu^3
S_total = Si*K;
V_total = S_total/Psi;
t = 0;          % Start time
t_print = 10;   % Time interval for printing in s
print_im = 'off';% Sets whether graphs are printed every t_print 
k = 0;          % Event counter
t_progress = 0; % 
scale_factor = 1; % Rescales the number of particles in system
nM = 6.022 * 10^23 * 10^-15 * 10^-9; % particles per mu^3: 1mol*10^-15mu^3*nano
% Nates concentrations
P_molar = 20;%79; %nm
J_conv = P_molar * 0.60221413; %convert nM 
A_molar = P_molar * 1.56; % nm * density ratio (see Paper table S1)
A_tot = V_total * 0.6022 * A_molar * scale_factor; % Total A
P_tot = V_total * 0.6022 * P_molar * scale_factor; % Total B

% My concentrations
% A_cyto  = floor(145 * nM * scale_factor * V_total);
% P_cyto  = floor(41 * nM * scale_factor * V_total);
% A_membr = (177 * nM - 145 * nM) * V_total * scale_factor;
% P_membr = (79 * nM - 41 * nM) * V_total * scale_factor;
%%
% First species
alpha       = 1;
D_A         = 0.28;         % mu^2/s
d_A         = D_A/h^2;
k_off_A     = 0.0054; 
k_on_A      = 0.00858;
k_on_A_s    = k_on_A * Si/(V_total);         
k_ant_AP    = 0.19/(Si*J_conv)^alpha; 
Conc_A      = zeros(1, K);  % Initial conditions
% Nates concentrations
Conc_A(1,1:25) = floor(k_on_A_s * A_tot /(k_on_A_s*K+k_off_A));
A_cyto = A_tot - sum(Conc_A);
% My concentrations
% Conc_A(1:floor(K/2)) = floor(A_membr/floor(K/2));
% Conc_A(floor(K/2) + 1) = floor(mod(A_membr,floor(K/2)));

% Second species
beta        = 2;
D_P         = 0.150;         % mu^2/s
d_P         = D_P/h^2;
k_off_P     = 0.0073; 
k_on_P      = 0.0474;
k_on_P_s    = k_on_P * Si/(V_total);         
k_ant_PA    = 2.0/(Si*J_conv)^beta;
Conc_P      = zeros(1, K);  % Initial conditions 
% Nates concentrations
Conc_P(26:50) = floor(k_on_P_s * P_tot /(k_on_P_s*K+k_off_P));
P_cyto = P_tot - sum(Conc_P);
% My concentrations
% Conc_P(ceil(K/2):end) = floor(P_membr/floor(K/2));
% Conc_P(ceil(K/2) - 1) = floor(mod(P_membr,floor(K/2)));

%%
tic;    % Start timer
%%% Gillespie ftw
while t < T
    r_1 = rand(1);
    r_2 = rand(1);
    
   %%% Propensity functions
    alpha_A_diff = Conc_A * d_A;
   % alpha_A_diff = gpuArray(alpha_A_diff);
    alpha_P_diff = Conc_P * d_P;
   % alpha_P_diff = gpuArray(alpha_P_diff);
    alpha_A_off  = Conc_A * k_off_A;
    %alpha_A_off  = gpuArray(alpha_A_off);
    alpha_P_off  = Conc_P * k_off_P;
   % alpha_P_off  = gpuArray(alpha_P_off);
    alpha_A_on   = ones(1,K) * A_cyto * k_on_A_s; 
   % alpha_A_on   = gpuArray(alpha_A_on);
    alpha_P_on   = ones(1,K) * P_cyto * k_on_P_s; 
   % alpha_P_on   = gpuArray(alpha_P_on);
    alpha_A_ant  = Conc_A .* Conc_P.^alpha * k_ant_AP;
   % alpha_A_ant  = gpuArray(alpha_A_ant);
    alpha_P_ant  = Conc_P .* Conc_A.^beta * k_ant_PA;
   % alpha_P_ant  = gpuArray(alpha_P_ant);
    
    sum_A_diff_1    = sum(alpha_A_diff(1:end-1));
    sum_A_diff_2    = sum(alpha_A_diff(2:end));
    sum_P_diff_1    = sum(alpha_P_diff(1:end-1));
    sum_P_diff_2    = sum(alpha_P_diff(2:end));
    sum_A_off       = sum(alpha_A_off(1:end));
    sum_P_off       = sum(alpha_P_off(1:end));
    sum_A_on        = sum(alpha_A_on(1:end));
    sum_P_on        = sum(alpha_P_on(1:end));
    sum_A_ant       = sum(alpha_A_ant(1:end));
    sum_P_ant       = sum(alpha_P_ant(1:end));
    
   %%% Sum of propensity functions of all reactions
    alpha_0 =   sum_A_diff_1 + sum_A_diff_2 + sum_P_diff_1 + sum_P_diff_2 + ...
                sum_A_off + sum_P_off + ...
                sum_A_on + sum_P_on + ...
                sum_A_ant + sum_P_ant;
            
   %%% Cumulative sum of all propensities      
    alpha_cum_sum = cumsum([alpha_A_diff(1:end-1)/alpha_0, ...
                            alpha_A_diff(2:end)/alpha_0, ...
                            alpha_P_diff(1:end-1)/alpha_0, ...
                            alpha_P_diff(2:end)/alpha_0, ...
                            alpha_A_off(1:end)/alpha_0,...
                            alpha_P_off(1:end)/alpha_0,...
                            alpha_A_on/(alpha_0),...
                            alpha_P_on/(alpha_0),...
                            alpha_A_ant/alpha_0,...
                            alpha_P_ant/alpha_0]); 
                            
   %%% Next time step
    tau = 1/alpha_0*log(1/r_1);

   %%% Probability intervals (cumulative, increasing)
    cum_sum_A_diff_1    =                     sum_A_diff_1/alpha_0;
    cum_sum_A_diff_2    = cum_sum_A_diff_1  + sum_A_diff_2/alpha_0;
    cum_sum_P_diff_1    = cum_sum_A_diff_2  + sum_P_diff_1/alpha_0;
    cum_sum_P_diff_2    = cum_sum_P_diff_1  + sum_P_diff_2/alpha_0;
    cum_sum_A_off       = cum_sum_P_diff_2 	+ sum_A_off/alpha_0;
    cum_sum_P_off       = cum_sum_A_off     + sum_P_off/alpha_0;
    cum_sum_A_on        = cum_sum_P_off 	+ sum_A_on/alpha_0;
    cum_sum_P_on        = cum_sum_A_on      + sum_P_on/alpha_0;
    cum_sum_A_ant       = cum_sum_P_on      + sum_A_ant/alpha_0;
    cum_sum_P_ant       = cum_sum_A_ant     + sum_P_ant/alpha_0;
       
   %%% Choose next reaction
    %%% Diffusion
    if r_2 < cum_sum_A_diff_1
        j = find(alpha_cum_sum > r_2, 1, 'first'); 
        Conc_A(j) = Conc_A(j) - 1;
        Conc_A(j + 1) = Conc_A(j + 1) + 1;
    elseif r_2 >= cum_sum_A_diff_1 && r_2 < cum_sum_A_diff_2
        j = find(alpha_cum_sum > r_2, 1, 'first') - K + 2; 
        Conc_A(j) = Conc_A(j) - 1;
        Conc_A(j - 1) = Conc_A(j-1) + 1;
    elseif r_2 >= cum_sum_A_diff_2 && r_2 < cum_sum_P_diff_1
        j = find(alpha_cum_sum > r_2, 1, 'first') - 2*K + 2;
        Conc_P(j) = Conc_P(j) - 1;
        Conc_P(j + 1) = Conc_P(j + 1) + 1;
    elseif r_2 >= cum_sum_P_diff_1 && r_2 < cum_sum_P_diff_2
        j = find(alpha_cum_sum > r_2, 1, 'first') - 3*K + 4;
        Conc_P(j) = Conc_P(j) - 1;
        Conc_P(j - 1) = Conc_P(j - 1) + 1; 
     %%% Off rate
     elseif r_2 >= cum_sum_P_diff_2 && r_2 < cum_sum_A_off
         j = find(alpha_cum_sum > r_2, 1, 'first') - 4*K + 4;
        Conc_A(j) = Conc_A(j) - 1;
        A_cyto = A_cyto + 1;
     elseif r_2 >= cum_sum_A_off && r_2 < cum_sum_P_off
        j = find(alpha_cum_sum > r_2, 1, 'first') - 5*K + 4;
        Conc_P(j) = Conc_P(j) - 1;
        P_cyto = P_cyto + 1;
    %%% On rate
    elseif r_2 >= cum_sum_P_off && r_2 < cum_sum_A_on
        j = find(alpha_cum_sum > r_2, 1, 'first') - 6*K + 4;
        Conc_A(j) = Conc_A(j) + 1;
        A_cyto = A_cyto - 1;
    elseif r_2 >= cum_sum_A_on && r_2 < cum_sum_P_on
        j = find(alpha_cum_sum > r_2, 1, 'first') - 7*K + 4;
        Conc_P(j) = Conc_P(j) + 1;
        P_cyto = P_cyto - 1;
    %%% Antagonism
    elseif r_2 >= cum_sum_P_on && r_2 < cum_sum_A_ant
        j = find(alpha_cum_sum > r_2, 1, 'first') - 8*K + 4;
        Conc_A(j) = Conc_A(j) - 1;
        A_cyto = A_cyto + 1;
    elseif r_2 >= cum_sum_A_ant && r_2 < cum_sum_P_ant
        j = find(alpha_cum_sum > r_2, 1, 'first') - 9*K + 4;
        Conc_P(j) = Conc_P(j) - 1;
        P_cyto = P_cyto + 1;
    else
        disp('Problem here is.');
        break;
    end
    
    %%% Print every t_tolerance seconds
%     if t+tau >= t_print && strcmp(print_im, 'on')
%         f = figure('visible', 'off');
%         hold on;
%         plot(1:K, Conc_A, 'r');
%         plot(1:K, Conc_P, 'b');
%         bar(1:K, Conc_A, 'r', 'visible', 'off');
%         bar(1:K, Conc_P, 'b', 'visible', 'off');
%         hold off;
%         axis([0, 100, 0, 100]);
%         M(k) = getframe(gcf);
%         saveas(f, [num2str(k),'.png']);
%         t_print = t_print + 1;
%     end
%     
%     %%% Plot tolerance
%     if t > t_progress
%         disp(sum(Conc_A));
%         disp(sum(Conc_P));
%         disp(sum(Conc_A) + sum(A_cyto));
%         disp(sum(Conc_P) + sum(P_cyto));
%         close all; 
%         disp(t); 
%         hold on;
%         bar(1:K, Conc_A, 'r');
%         bar(1:K, Conc_P, 'b');
%         hold off;
%         drawnow;
%         t_progress = t_progress + t_plot;
%     end
     
    %%% Update time
    t = t + tau;
    k = k + 1;
end
        hold on;
        plot(1:K, Conc_A, 'r');
        plot(1:K, Conc_P, 'b');
        bar(1:K, Conc_A, 'r', 'visible', 'off');
        bar(1:K, Conc_P, 'b', 'visible', 'off');
        hold off;
time_elapsed = toc;
shg;
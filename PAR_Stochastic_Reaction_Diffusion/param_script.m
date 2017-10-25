numSimus  = 10; %Number of times each simulation is carried out

simparams.totalT    = 3000;
simparams.bins      = 50;
simparams.boundPos  = 25;

params.L        = 134.6;
params.psi      = 0.174;
params.J        = 79;
params.Bconc 	= 79;
params.DA       = 0.28;
params.DB       = 0.15;
params.alpha    = 1;
params.beta     = 2;
params.koffA    = 0.0054; % A > cytoA
params.konA     = 0.00858; % cytoA > A
params.koffB    = 0.0073; % B > cytoB
params.konB     = 0.0474; % cytoB > B
params.kcA      = 0.19; % A + alpha*B > Acyto + alpha*B
params.kcB      = 2; % B + beta*A > Bcyto + beta*A
params.p        = 0.0;

%% CHANGE LENGTH AND INITIAL DOMAIN SIZE
L_arr   = 20:20:300;
p_arr   = 0.05:0.05:0.25;

params_t = struct2table(params);
for i = 1 : length(L_arr)
    params_t.L = L_arr(i);
    for j = 1 : length(p_arr)
        params_t.p = p_arr(j);
        T((i-1)*(length(p_arr))+j, :) = params_t(1, :);
    end
end
%%
tic
if isempty(gcp('nocreate')); parpool(8); end;
L_arr = 30:10:300;
params_t = struct2table(params);

for i = 1:length(L_arr)
    params_t.L = L_arr(i);
    T(i, :) = params_t(1, :);
end

% parfor i = 1 : length(L_arr)
%     s(i) = PARsRDv1_0('Domain_stability', 3, simparams, T(i, :));
% end
s = cell(2,1);
if isempty(gcp('nocreate')); parpool(8); end;
for i = 1 : length(L_arr)
Sim_cell = cell(numSimus,length(p_arr));
    for j = 1 : length(p_arr)
        for k = 1 : numSimus
            Sim_cell{k,j} = PARsRDv1_0('NightlyRun3', 3, simparams, T((i-1)*(length(p_arr))+j, :));
        end
    end
    s(i) = {Sim_cell};
end
toc
%%
params_a(4) = params;
for i = 1:length(params_a)
    params_a(i) = params;
end
%% Wildtype
params_a(1).psi = 0.162;
params_a(1).L = 137.1;

params_a(2).psi = 0.278;
params_a(2).L = 89.0;

params_a(3).psi = 0.321;
params_a(3).L = 62.7;

params_a(4).psi = 0.497;
params_a(4).L = 41.5;
% 
% %% 1.5 folder bigger diffusion coefficients
% params_a(5).psi = 0.162;
% params_a(5).L = 137.1/2;
% params_a(5).DB       = 0.225;
% 
% params_a(6).psi = 0.278;
% params_a(6).L = 89.0/2;
% params_a(6).DB       = 0.225;
% 
% params_a(7).psi = 0.321;
% params_a(7).L = 62.7/2;
% params_a(7).DB       = 0.225;
% 
% params_a(8).psi = 0.497;
% params_a(8).L = 41.5/2;
% params_a(8).DB       = 0.225;

%%
% for i = 1:length(L_arr)
%     params_a(1, i).DB      = params.DB;
%     params_a(2, i).DB      = 2 * params.DB;
%     params_a(3, i).DB      = 0.5 * params.DB;
% 
%     params_a(4, i).koffB   = 2 * params.koffB;
%     params_a(5, i).koffB   = 0.5 * params.koffB;
% 
%     params_a(6, i).konB   = 2 * params.konB;
%     params_a(7, i).konB   = 0.5 * params.konB;
%     for j = 1 : 7
%         params_a(j, i).L = L_arr(i)/2;
%         params_a(j, i).psi = StoV(i);
%     end
% end
L_arr = 30:10:300;
for i = 1:length(L_arr)
    for j = 1 : 7
        params_a(j, i).L = L_arr(i)/2;
        params_a(j, i).psi = StoV(i);
    end
end
tic
% parpool(8);
parfor i = 1:8
    for j = 1 : 4
        s{i} = PARsRDv1_0(['todelete' num2str(j), '_', num2str(i)], 2, 'PER', simparams, params_a(j));
    end
end
% delete(gcp);
toc
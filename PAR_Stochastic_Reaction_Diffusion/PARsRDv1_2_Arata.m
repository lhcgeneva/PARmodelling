function Simulation = PARsRDv1_0(name, numDom, BoundaryC, simparams, params)

% Stochastic RD model based on the PAR RD system described in Goehring et
% al (2011).  
%
% Last updated 4 NOV 2014
%
% This version responds as expected to changes in bin and protein amounts.
% Accounts for normalisation factor J (i.e. reverse normalisation of RhoP =
% 1 used in the paper) and conversion from concentrations to N.

if nargin > 3
    totalT  = simparams.totalT;
%     Either variable or fix bin size
%     if keepBinSize
%         bins = params.L/1;%1 micron/bin
%         boundPos = round(bins/2);
%     else
%         bins = simparams.bins;
%         simparams.boundPos;
%     end
    bins    = params.bins;
    boundPos= round(bins/2);
    L       = params.L;
    psi     = params.psi;
    J       = params.J;
    Bconc   = params.Bconc;
    DA      = params.DA;
    DB      = params.DB;
    alpha   = params.alpha;
    beta    = params.beta;
    koffA   = params.koffA; % A > cytoA
    konA    = params.konA; % cytoA > A
    koffB   = params.koffB; % B > cytoB
    konB    = params.konB; % cytoB > B
    kcA     = params.kcA; % A + alpha*B > Acyto + alpha*B
    kcB     = params.kcB; % B + beta*A > Bcyto + beta*A
    p       = params.p; % Percentage of membrane taken by B in beginning
    ratioAB = params.ratioAB;
else
    totalT = 10; % Total simulation time
    bins = 100; % Discretization bins
    boundPos = 25; % Initial boundary position in bin #
    
    L   = 134.6/2; %P0: L = 134.6/2; % Length of system in um
    psi = 0.174; % S/V ration from Goehring et al (2011).
    J   = 79; 
    
    Bconc   = 79; % in nM
    DA      = 0.28;
    DB      = 0.15;
    alpha   = 1;
    beta    = 2;
    koffA   = 0.0054; % A > cytoA
    konA    = 0.00858; % cytoA > A
    koffB   = 0.0073; % B > cytoB
    konB    = 0.0474; % cytoB > B
    kcA     = 0.19; % A + alpha*B > Acyto + alpha*B
    kcB     = 2; % B + beta*A > Bcyto + beta*A
    p       = 0; % Percentage of membrane taken by B in beginning
    ratioAB = 1.56;
end

% 

J = J * 0.60221413; % Normalisation factor equal [PAR-2] in N / um^3
Bconc = Bconc * 0.60221413; %  nM to N / um^3
Aconc = Bconc * ratioAB; % Aconc calculated from PAR2/6 ratio

% Parameter adjustments and calculations

h = L/bins; % Bin length
w = 1; % Bin width
Si = h*w; % Area of bin
totalS = Si * bins; % Total surface area
totalV = totalS / psi; % Total volume (assuming 

dA = DA/h^2;
dB = DB/h^2;
% The below adjustment for the on rate is per compartment. Check out how on
% rate is summed for every compartment below, if in doubt about S/V
% conversion
konA = konA * Si / totalV; % Adjust for concentration to N conversion
konB = konB * Si / totalV; % Adjust for concentration to N conversion
kcA = kcA / J^alpha / Si^alpha; % Adjust for concentration to N conversion
kcB = kcB / J^beta / Si^beta; % Adjust for concentration to N conversion

% Initial conditions

% Plot interval definition and counter
k = 1;
plotInt = totalT/10;
plotInts = 0:plotInt:totalT;

% Initiate time variable
t = 0;
plotTime(1) = t;

% Calculate total amounts
A0 = totalV * Aconc; % Total A
B0 = totalV * Bconc; % Total B

% Initialize A and B distributions
A(:,1) = zeros(bins,1);
B(:,1) = zeros(bins,1);

%%% Initial conditions
if numDom == 2
    A(1:boundPos,1) = 1 * floor(konA * A0 / (konA*bins+koffA)); % Init distribution 
    B(boundPos+1:end,1) = 1 * floor(konB * B0 / (konB*bins+koffB));
elseif numDom == 4
    A(1:floor(boundPos/2),1) = 1 * floor(konA * A0 / (konA*bins+koffA)); % Init distribution 
    B(floor(boundPos/2)+1:2*floor(boundPos/2),1) = 1 * floor(konB * B0 / (konB*bins+koffB));
    A(2*floor(boundPos/2)+ 1:3*floor(boundPos/2),1) = 1 * floor(konA * A0 / (konA*bins+koffA)); % Init distribution 
    B(3*floor(boundPos/2)+1:end,1) = 1 * floor(konB * B0 / (konB*bins+koffB));
elseif numDom == 3
    A(floor(p * bins) + 1:floor((0.5+p)*bins), 1) = 1 * floor(konA * A0 / (konA*bins+koffA));
    B(1:floor(p * bins),1) = 1 * floor(konB * B0 / (konB*bins+koffB));
    B(floor((0.5+p)*bins)+1:end,1) = 1 * floor(konB * B0 / (konB*bins+koffB));
else disp('Number of domains not supported');
end
% Define initial values for species A, cytoA, B, cytoB
cytoA = A0 - sum(A(:,1));
cytoB = B0 - sum(B(:,1));

nowA = A(:,1);
nowB = B(:,1);

temp = zeros(8, bins);

% %%%% ONLY DIFFUSION, NO REACTIONS, COMMENT OUT FOR NORMAL RUNS %%%%%%%%%%%%%%%%%
%     koffA   = 0.00; % A > cytoA
%     konA    = 0.00; % cytoA > A
%     koffB   = 0.00; % B > cytoB
%     konB    = 0.00; % cytoB > B
%     kcA     = 0.00; % A + alpha*B > Acyto + alpha*B
%     kcB     = 0.00; % B + beta*A > Bcyto + beta*A
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa_t = 0;
bb_t = 0;
cc_t = 0;
dd_t = 0;
aa_c = 0;
bb_c = 0;
cc_c = 0;
dd_c = 0;
% Simulation

h = @(x) 0.5*erf((-x+50)/50)+1.3;

while t < totalT
    
    % Select 2 random numbers
    randA = rand;
    randB = rand;
    
    % Calculate exchange probabilities
    offA = nowA .* (koffA + kcA * nowB.^alpha);
    offB = nowB .* (koffB.*h(nowA) + kcB * nowA.^beta); 
    onA = ones(bins,1) * konA * cytoA;
    onB = ones(bins,1) * konB * cytoB;
    
    % Calculate diffusion probabilities
    difA = nowA .* dA;
    difB = nowB .* dB;
            
    % Cummulative propensities
    
    offAsum = cumsum(offA);
    offBsum = cumsum(offB);
    onAsum = cumsum(onA);
    onBsum = cumsum(onB);

    
    %%%PERIODIC OR REFLECTIVE BOUNDARY CONDITIONS%%%
    if strcmp(BoundaryC, 'PER')
        difAR = cumsum(difA(1:end));
        difAL = cumsum(difA(1:end));
        difBR = cumsum(difB(1:end));
        difBL = cumsum(difB(1:end));
    elseif strcmp(BoundaryC, 'REF')
        difAR = cumsum(difA(1:end-1));
        difAL = cumsum(difA(2:end));
        difBR = cumsum(difB(1:end-1));
        difBL = cumsum(difB(2:end));
    else
        disp('Boundary condition type not supported');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calculate probability variables
    a0 = offAsum(end) + offBsum(end) + onAsum(end)+ onBsum(end)...
        + difAR(end) + difAL(end) + difBR(end) + difBL(end);
    
    offAP = offAsum(end) / a0;
    offBP = offBsum(end) / a0;
    onAP = onAsum(end) / a0;
    onBP = onBsum(end) / a0;
    difARP = difAR(end) / a0;
    difALP = difAL(end) / a0;
    difBRP = difBR(end) / a0;
    difBLP = difBL(end) / a0;    
    
    % Find next timestep
    tau = 1 / a0 * log(1/randA);

%%%%%%%%%%%%%%%%     Choose process     %%%%%%%%%%%%%%%%%%%
    if randB < offAP
        probMatrix = offAsum/a0;
        j = find(probMatrix >= randB, 1, 'first');
        nowA(j) = nowA(j) - 1;
        cytoA = cytoA + 1;
       
    elseif randB < offAP + offBP
        probMatrix = offBsum/a0 + offAP;
        j = find(probMatrix >= randB, 1, 'first');
        nowB(j) = nowB(j) - 1;
        cytoB = cytoB + 1;
        
    elseif randB < offAP + offBP + onAP
        probMatrix = onAsum/a0 + offAP + offBP;
        j = find(probMatrix >= randB, 1, 'first');
        nowA(j) = nowA(j) + 1;
        cytoA = cytoA - 1;
        
    elseif randB < offAP + offBP + onAP + onBP
        probMatrix = onBsum/a0 + offAP + offBP + onAP;
        j = find(probMatrix >= randB, 1, 'first');
        nowB(j) = nowB(j) + 1;
        cytoB = cytoB - 1;
        
    %DIFFUSION
    elseif randB >= offAP + offBP + onAP + onBP
        if strcmp('REF', BoundaryC)
            if randB < offAP + offBP + onAP + onBP + difARP
                probMatrix = difAR/a0 + offAP + offBP + onAP + onBP;
                j = find(probMatrix >= randB, 1, 'first');
                nowA(j) = nowA(j) - 1;
                nowA(j+1) = nowA(j+1) + 1;
                aa_t = aa_t + 1;
            elseif randB < offAP + offBP + onAP + onBP + difARP + difALP
                probMatrix = difAL/a0 + offAP + offBP + onAP + onBP + difARP;
                j = find(probMatrix >= randB, 1, 'first');
                j = j+1;
                nowA(j) = nowA(j) - 1;
                nowA(j-1) = nowA(j-1) + 1;
                bb_t = bb_t + 1;
            elseif randB < offAP + offBP + onAP + onBP + difARP + difALP + difBRP
                probMatrix = difBR/a0 + offAP + offBP + onAP + onBP + difARP + difALP;
                j = find(probMatrix >= randB, 1, 'first');
                nowB(j) = nowB(j) - 1;
                nowB(j+1) = nowB(j+1) + 1;
                cc_t = cc_t + 1;
            elseif randB <= 1
                probMatrix = difBL/a0 + offAP + offBP + onAP + onBP + difARP + difALP + difBRP;
                j = find(probMatrix >= randB, 1, 'first');
                j = j+1;
                nowB(j) = nowB(j) - 1;
                nowB(j-1) = nowB(j-1) + 1; 
                dd_t = dd_t + 1;  
            end
        elseif strcmp('PER', BoundaryC)
            if randB < offAP + offBP + onAP + onBP + difARP
                probMatrix = difAR/a0 + offAP + offBP + onAP + onBP;
                j = find(probMatrix >= randB, 1, 'first');
                nowA(j) = nowA(j) - 1;
                try 
                    nowA(j+1) = nowA(j+1) + 1;
                    aa_t = aa_t + 1;
                catch
                    nowA(1) = nowA(1) + 1;
                    aa_c = aa_c + 1;
                end
            elseif randB < offAP + offBP + onAP + onBP + difARP + difALP
                probMatrix = difAL/a0 + offAP + offBP + onAP + onBP + difARP;
                j = find(probMatrix >= randB, 1, 'first');
                nowA(j) = nowA(j) - 1;
                try 
                    nowA(j-1) = nowA(j-1) + 1;
                    bb_t = bb_t + 1;
                catch 
                    nowA(end) = nowA(end) + 1;
                    bb_c = bb_c + 1;
                end
            elseif randB < offAP + offBP + onAP + onBP + difARP + difALP + difBRP
                probMatrix = difBR/a0 + offAP + offBP + onAP + onBP + difARP + difALP;
                j = find(probMatrix >= randB, 1, 'first');
                nowB(j) = nowB(j) - 1;
                try
                    nowB(j+1) = nowB(j+1) + 1;
                    cc_t = cc_t + 1;
                catch
                    nowB(1) = nowB(1) + 1;
                    cc_c = cc_c + 1;
                end
            elseif randB <= 1
                probMatrix = difBL/a0 + offAP + offBP + onAP + onBP + difARP + difALP + difBRP;
                j = find(probMatrix >= randB, 1, 'first');
                nowB(j) = nowB(j) - 1;
                try 
                    nowB(j-1) = nowB(j-1) + 1;
                    dd_t = dd_t + 1;
                catch 
                    nowB(end) = nowB(end) + 1;
                    dd_c = dd_c + 1;
                end
            end
        end
    end
    
%%%%%%%%%%%%%%   Update time, output select timepoints   %%%%%%%%%%%%%%%%%%
    t = t + tau; 
    
    if t > plotInts(k)
        k = k + 1;
        A(:,k) = nowA;
        B(:,k) = nowB;
        plotTime(k) = t;
%         disp(['simulation time: ', num2str(t)]);
    end

end
% aa_t + aa_c
% bb_t + bb_c
% cc_t + cc_c
% dd_t + dd_c
% aa_c
% bb_c
% cc_c
% dd_c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

xpos = linspace(0, 1, bins);
endA = A(:,end);
endB = B(:,end);

subplot(1,2,1), plot(xpos, endA,'r-'), hold on, plot(xpos', endB,'b-'), xlim([0 1]), axis square
xlabel('Position in cell (a.u.)');
ylabel('Concentration (a.u.)');
I = cat(3,A'/max(max(A)),B'/max(max(B)),zeros(size(A')));
subplot(1,2,2), imagesc(I), axis square;
xlabel('Position in cell (a.u.)');
ylabel('Time (in plotInterval * s)');

Simulation.name = name;
if nargin > 2
    Simulation.params = params;
    Simulation.A    = A;
    Simulation.B    = B;
    Simulation.plotTime = plotTime;
end
end
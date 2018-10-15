function [A, B, plotTime] = PARsRDv1_0_NateUnchanged

% Stochastic RD model based on the PAR RD system described in Goehring et
% al (2014).  
%
% Last updated 4 NOV 2014
%
% This version responds as expected to changes in bin and protein amounts.
% Accounts for normalisation factor J (i.e. reverse normalisation of RhoP =
% 1 used in the paper) and conversion from concentrations to N.



% Specified Parameters

totalT = 100; % Total simulation time
bins = 100; % Discretization bins
boundPos = 50; % Initial boundary position in bin #
L = 134.6/2; % Length of system in um
psi = 0.174; % S/V ration from Goehring et al (2014).

J = 79 * 0.60221413; % Normalisation factor equal [PAR-2] in N / um^3

ratioAB = 1.56;

Bconc = 79; % in nM
Bconc = Bconc * 0.60221413; %  nM to N / um^3
Aconc = Bconc * ratioAB; % Aconc calculated from PAR2/6 ratio

DA = 0.28;
DB = 0.15;
alpha = 1;
beta = 2;
koffA = 0.0054; % A > cytoA
konA = 0.00858; % cytoA > A
koffB = 0.0073; % B > cytoB
konB = 0.0474; % cytoB > B
kcA = 0.19; % A + alpha*B > Acyto + alpha*B
kcB = 2; % B + beta*A > Bcyto + beta*A

% Parameter adjustments and calculations

h = L/bins; % Bin length
w = 1; % Bin width
Si = h*w; % Area of bin
totalS = Si * bins; % Total surface area
totalV = totalS / psi; % Total volume (assuming 

dA = DA/h^2;
dB = DB/h^2;
konA = konA * Si / totalV; % Adjust for concentration to N conversion
konB = konB * Si / totalV; % Adjust for concentration to N conversion
kcA = kcA / J^alpha / Si^alpha; % Adjust for concentration to N conversion
kcB = kcB / J^beta / Si^beta; % Adjust for concentration to N conversion

% Initial conditions

% Plot interval definition and counter
k = 1;
plotInt = totalT/100;
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

A(1:boundPos,1) = 1 * floor(konA * A0 / (konA*bins+koffA)); % Init distribution 
B(boundPos+1:end,1) = 1 * floor(konB * B0 / (konB*bins+koffB));


% Define initial values for species A, cytoA, B, cytoB
cytoA = A0 - sum(A(:,1));
cytoB = B0 - sum(B(:,1));

nowA = A(:,1);
nowB = B(:,1);

% Simulation
while t < totalT
    
    % Select 2 random numbers
    randA = rand;
    randB = rand;
    
    % Calculate exchange probabilities
    offA = nowA .* (koffA + kcA * nowB.^alpha);
    offB = nowB .* (koffB + kcB * nowA.^beta); 
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
    difAR = cumsum(difA(1:end-1));
    difAL = cumsum(difA(2:end));
    difBR = cumsum(difB(1:end-1));
    difBL = cumsum(difB(2:end));

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

    % Choose which function
    if randB < offAP
%         display(['-A '])  %  Display action
        probMatrix = offAsum/a0;
        j = find(probMatrix >= randB, 1, 'first');
        nowA(j) = nowA(j) - 1;
        cytoA = cytoA + 1;
       
    elseif randB < offAP + offBP
%         display(['-B '])  %  Display action
        probMatrix = offBsum/a0 + offAP;
        j = find(probMatrix >= randB, 1, 'first');
        nowB(j) = nowB(j) - 1;
        cytoB = cytoB + 1;
        
    elseif randB < offAP + offBP + onAP
%         display(['+A '])  %  Display action
        probMatrix = onAsum/a0 + offAP + offBP;
        j = find(probMatrix >= randB, 1, 'first');
        nowA(j) = nowA(j) + 1;
        cytoA = cytoA - 1;
        
    elseif randB < offAP + offBP + onAP + onBP
%         display(['+B '])  %  Display action
        probMatrix = onBsum/a0 + offAP + offBP + onAP;
        j = find(probMatrix >= randB, 1, 'first');
        nowB(j) = nowB(j) + 1;
        cytoB = cytoB - 1;
        
    elseif randB < offAP + offBP + onAP + onBP + difARP
%         display(['A Right '])  %  Display action
        probMatrix = difAR/a0 + offAP + offBP + onAP + onBP;
        j = find(probMatrix >= randB, 1, 'first');
        nowA(j) = nowA(j) - 1;
        nowA(j+1) = nowA(j+1) + 1;
        
    elseif randB < offAP + offBP + onAP + onBP + difARP + difALP
%         display(['A Left '])  %  Display action
        probMatrix = difAL/a0 + offAP + offBP + onAP + onBP + difARP;
        j = find(probMatrix >= randB, 1, 'first');
        j = j+1;
        nowA(j) = nowA(j) - 1;
        nowA(j-1) = nowA(j-1) + 1;
        
    elseif randB < offAP + offBP + onAP + onBP + difARP + difALP + difBRP
%         display(['B Right '])  %  Display action
        probMatrix = difBR/a0 + offAP + offBP + onAP + onBP + difARP + difALP;
        j = find(probMatrix >= randB, 1, 'first');
        nowB(j) = nowB(j) - 1;
        nowB(j+1) = nowB(j+1) + 1;
        
    elseif randB >= offAP + offBP + onAP + onBP + difARP + difALP + difBRP
%         display(['B Left '])  %  Display action
        probMatrix = difBL/a0 + offAP + offBP + onAP + onBP + difARP + difALP + difBRP;
        j = find(probMatrix >= randB, 1, 'first');
        j = j+1;
        nowB(j) = nowB(j) - 1;
        nowB(j-1) = nowB(j-1) + 1;    
        
    end
    
    % Update time
    t = t + tau; 
    
    % Output only select timepoints
    if t > plotInts(k)
        k = k + 1;
        A(:,k) = nowA;
        B(:,k) = nowB;
        plotTime(k) = t;
        disp(['simulation time: ', num2str(t)]);
    end

end

% Plotting
close all

xpos = linspace(0, 1, bins);
endA = A(:,end);
endB = B(:,end);

subplot(1,2,1), plot(xpos, endA,'r-'), hold on, plot(xpos', endB,'g-'), xlim([0 1]), axis square

I = cat(3,A'/max(max(A)),B'/max(max(B)),zeros(size(A')));

subplot(1,2,2), imagesc(I), axis square

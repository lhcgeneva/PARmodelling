% This creates example profiles for different diffusivities
%% WP
[u01, ~] = wave_pin(0.1, 100000, 15, 1/9, 1/9, 1, 1, 0.067/9);
[u0025, ~] = wave_pin(0.025, 100000, 15, 1/9, 1/9, 1, 1, 0.067/9);
[u02, ~] = wave_pin(0.2, 100000, 15, 1/9, 1/9, 1, 1, 0.067/9);
dlmwrite('WP.csv', [u0025(end, :); u01(end, :);u02(end, :)]);
%% Different sizes
[u15, ~] = wave_pin(0.1, 100000, 15, 1/9, 1/9, 1, 1, 0.067/9);
[u5, ~] = wave_pin(0.1, 100000, 5, 1/9, 1/9, 1, 1, 0.067/9);
[u25, ~] = wave_pin(0.1, 100000, 25, 1/9, 1/9, 1, 1, 0.067/9);
dlmwrite('WP_sizes.csv', [u5(end, :); u15(end, :);u25(end, :)]);
%% Write WP
%% OT
[uOT_01, v] = otsuji(15, 0.1, 100000, 1, 0.7, 1, 1.5);
[uOT_02, v] = otsuji(15, 0.2, 100000, 1, 0.7, 1, 1.5);
[uOT_0025, v] = otsuji(15, 0.025, 100000, 1, 0.7, 1, 1.5);
%% Plot OT
figure; hold on;
plot(uOT_01(end, :)');
plot(uOT_02(end, :)');
plot(uOT_0025(end, :)')
%% Write OT
dlmwrite('OT.csv', [uOT_0025(end, :); uOT_01(end, :);uOT_02(end, :)]);
%% Different sizes
[uOT_5, v] = otsuji(5, 0.1, 100000, 1, 0.7, 1, 1.5);
[uOT_15, v] = otsuji(15, 0.1, 100000, 1, 0.7, 1, 1.5);
[uOT_25, v] = otsuji(25, 0.1, 100000, 1, 0.7, 1, 1.5);
dlmwrite('OT_sizes.csv', [uOT_5(end, :); uOT_15(end, :);uOT_25(end, :)]);
%% GOR
[~, uGOR01, ~, ~] = goryachev(15, 0.1, 10000000, 0.67/101, 0.33/101, 0.01, 5);
[~, uGOR02, ~, ~] = goryachev(15, 0.2, 10000000, 0.67/101, 0.33/101, 0.01, 5);
[~, uGOR0025, ~, ~] = goryachev(15, 0.025, 10000000, 0.67/101, 0.33/101, 0.01, 5);
dlmwrite('GOR.csv', [uGOR0025(end, :); uGOR01(end, :);uGOR02(end, :)]);
%% Different sizes
[~, uGOR5, ~, ~] = goryachev(5, 0.1, 10000000, 0.67/101, 0.33/101, 0.01, 5);
[~, uGOR15, ~, ~] = goryachev(15, 0.1, 10000000, 0.67/101, 0.33/101, 0.01, 5);
[~, uGOR25, ~, ~] = goryachev(25, 0.1, 10000000, 0.67/101, 0.33/101, 0.01, 5);
dlmwrite('GOR_sizes.csv', [uGOR5(end, :); uGOR15(end, :);uGOR25(end, :)]);


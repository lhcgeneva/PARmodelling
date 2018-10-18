Da = logspace(-2, 1, 8);
grad5 = cell(1, length(Da));
for i = 1:length(Da)
    grad5{i} = wave_pin(Da(i), 40, 1);
end
%% Write different profiles to file to be plotted in python
ints = cell2mat(cellfun(@(x) x(end, :), grad5, 'Uni', 0)');
csvwrite('/Users/hubatsl/Desktop/interplay-cell-size/Theory/ints_wp.csv', ints);
%% Simulate for different Ds to calculate slopes
Da = logspace(-2, -0.45, 8);
grad5 = cell(1, length(Da));
for i = 1:length(Da)
    grad5{i} = wave_pin(Da(i), 40, 1);
end
%% Fit D and slope/gradient length
% figure(1)
% g = cellfun(@(x) x(end, :), grad5', 'Uni', 0);
% plot(cell2mat(g)')    
figure(2); hold on;
s = cellfun(@(x) max(abs(diff(x(end, :)))), grad5);
f = fit(Da(1:end)', 1./s(1:end)', 'a*sqrt(x)');
plot(Da(1:end), f(Da(1:end)), 'o');
plot(Da(1:end), 1./s(1:end));
%%
csvwrite('/Users/hubatsl/Desktop/interplay-cell-size/Theory/fit.csv', f(linspace(0, max(Da), 100)));
csvwrite('/Users/hubatsl/Desktop/interplay-cell-size/Theory/wave_pin.csv', s');
csvwrite('/Users/hubatsl/Desktop/interplay-cell-size/Theory/Da.csv', Da');
%%
% figure
% plot(cell2mat(grad2fold')')  
Da = 4*2*(0.1:0.1:0.92);
[s, ind] = cellfun(@(x) max(abs(diff(x))), grad2fold, 'Uni', 0);
s = cellfun(@(i, m, x) m/x(i), ind, s, grad2fold);
% plot(0.09*1./sqrt(1:35));
f = fit(Da', s', 'a*sqrt(x)');
figure; hold on;
plot(Da, s);
% plot(Da, 0.1*sqrt(Da));
plot(Da, f(Da));
%%
figure
plot(diff(grad(1, :)))
%% L ~ sqrt(Dt), t=1/delta
Da = 2*(0.1:0.1:0.92);
grad = cell(1, length(Da));
for i = 1:length(Da)
    grad{i} = wave_pin(Da(i), 10);
end
%%
figure; hold on;
for i = 1:length(Da)
    plot(grad{i}(end, :));
end
%% 2x length, 4x D
Da = 4*2*(0.1:0.1:0.92);
grad2fold = cell(1, length(Da));
for i = 1:length(Da)
    grad2fold{i} = wave_pin(Da(i), 20);
end
%%
figure; hold on;
for i = 1:length(Da)
    plot(grad2fold{i}(end, :));
end
%% 4x length, 16x D
Da = 4*4*2*(0.1:0.1:0.92);
grad2fold = cell(1, length(Da));
for i = 1:length(Da)
    grad2fold{i} = wave_pin(Da(i), 40);
end
%%
figure; hold on;
for i = 1:length(Da)
    plot(grad2fold{i}(end, :));
end
%% 8x length, 64x D
Da = 4*4*4*2*(0.1:0.1:0.92);
grad2fold = cell(1, length(Da));
for i = 1:length(Da)
    grad2fold{i} = wave_pin(Da(i), 80);
end
%%
figure; hold on;
for i = 1:length(Da)
    plot(grad2fold{i}(end, :));
end
%% 16x length, 256x D
Da = 4*4*4*4*2*(0.1:0.1:0.92);
grad2fold = cell(1, length(Da));
for i = 1:length(Da)
    grad2fold{i} = wave_pin(Da(i), 160);
end
%%
figure; hold on;
for i = 1:length(Da)
    plot(grad2fold{i}(end, :));
end
%% 20x length, 400x D
Da = 400*2*(0.1:0.1:0.92);
grad20fold = cell(1, length(Da));
for i = 1:length(Da)
    grad20fold{i} = wave_pin(Da(i), 200);
end
%%
figure; hold on;
for i = 1:length(Da)
    plot(grad20fold{i}(end, :));
end

%% Changing off rates
delta = logspace(-2, 1.5, 24);
off_change = cell(1, length(delta));
for i = 1:length(delta)
    off_change{i} = wave_pin(0.1, 10, delta(i));
end

%%
figure; hold on;
for i = 1:length(delta)
    plot(off_change{i}(end, :));
end
%% Write different off rate profiles to csv for plotting in python
intsoff = cell2mat(cellfun(@(x) x(end, :), off_change, 'Uni', 0)');
csvwrite('/Users/hubatsl/Desktop/interplay-cell-size/Theory/ints_wp_off.csv', intsoff);
%%
figure(6); hold on;
s = cellfun(@(x) max(abs(diff(x(end, :)))), off_change);
f = fit(delta(1:end)', s(1:end)', 'a*sqrt(x)');
plot(delta, f(delta), 'o');
plot(delta, s);
csvwrite('fit_off.csv', f(delta));
csvwrite('wave_pin_off.csv', s');
csvwrite('delta_off.csv', delta');
%% Breakdown size vs D
% Breakdown happens between 4th and 5th D, at 0.9, 3.6, 14.4 (40), 57.6 (80), and 360, calculating
% the proportionality constant from L_crit = A*sqrt(D) gives A =
% L_crit/sqrt(D) = 10/sqrt(0.9) = 10.5409
figure; hold on;
plot([0.9, 3.6, 14.4, 57.6, 230.4, 360], [10, 20, 40, 80, 160, 200], '.', 'MarkerSize', 20);
plot(10.5409*sqrt(0:1:360), '--', 'LineWidth', 2);
yticks([0:40:200])
set(gca,'linewidth',1.0)
set(gca, 'FontSize', 11.5);
xlabel('D [$\mu m^2/s$]', 'Interpreter', 'Latex');
ylabel('L$_{crit}\;[\mu m$]', 'Interpreter', 'Latex');
set(gca,'TickDir','out');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 3.75 3.1]);
% This plot is done in python Par_RD_PDE.ipynb
% print('/Users/hubatsl/Dropbox (Lars DMS)/LRI/PAR_Size/Theory/Lcrit_vsD_WavePin', '-dpdf');
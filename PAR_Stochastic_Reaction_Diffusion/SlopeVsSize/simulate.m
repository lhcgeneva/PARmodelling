% Not taking into account P1 cells
sizes = 10:5:100;
s = Simulations({{{'DB', 'DA', 'koffB', 'koffA', 'L', 'psi', 'bins'}},...
                {{fDP2(2*sizes), fDP6(2*sizes), fOP2(2*sizes), fOP6(2*sizes),...
                  sizes,S_to_V('Circumference', {2*sizes, 15/27}), 2*round(sizes)}}}, 20);
s.run_simulations();
save s3.mat
%% Taking into account P1 cells
sizes = 10:5:100;
s = Simulations({{{'DB', 'DA', 'koffB', 'koffA', 'L', 'psi', 'bins'}},...
                {{fDP2(2*sizes), fDP6(2*sizes), fOP2(2*sizes), fOP6(2*sizes),...
                  sizes, S_to_V('Circumference', {2*sizes, 15/27}), 2*round(sizes)}}}, 20);
s.run_simulations();
save s4.mat
%% Taking into account P1 cells and using ratio 1/1 for ellipsoid
sizes = 20:0.5:22;
StoV = 4*pi./(2*sizes/(2*pi)*4/3*pi);
s = Simulations({{{'DB', 'DA', 'koffB', 'koffA', 'L', 'psi', 'bins'}},...
                {{fDP2(2*sizes), fDP6(2*sizes), fOP2(2*sizes), fOP6(2*sizes),...
                  sizes, StoV, 2*round(sizes)}}}, 20);
s.run_simulations();
save s5.mat
%% 06/03/2018 Run OffRates.m before running this section, to access fDP2 etc.
sizes = [16, 18, 20, 22, 30, 134.6/2];
s = Simulations({{{'DB', 'DA', 'koffB', 'koffA', 'L', 'psi', 'bins'}},...
                {{fDP2(2*sizes), fDP6(2*sizes), fOP2(2*sizes), fOP6(2*sizes),...
                  sizes,S_to_V('Circumference', {2*sizes, 15/27}), round(sizes)}}}, 8);
s.run_simulations();
% save s_06_03_2017.mat
%% 02/04/2018 Also using P1 cells! Run OffRates.m before running this section, to access fDP2 etc.
sizes = [16, 18, 20, 22, 30, 134.6/2];
s = Simulations({{{'DB', 'DA', 'koffB', 'koffA', 'L', 'psi', 'bins'}},...
                {{fDP2(2*sizes), fDP6(2*sizes), fOP2(2*sizes), fOP6(2*sizes),...
                  sizes,S_to_V('Circumference', {2*sizes, 15/27}), round(sizes)}}}, 8);
s.run_simulations();
save s_02_04_2018.mat
%% Assign data from deterministic model for plotting
load /Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/s_02_04_2018.mat;
s.param_table_summary.L = s.param_table_summary.L*2 % to display full length;
a16 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Anterior16.csv');
b16 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Posterior16.csv');
a18 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Anterior18.csv');
b18 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Posterior18.csv');
a20 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Anterior20.csv');
b20 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Posterior20.csv');
a22 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Anterior22.csv');
b22 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Posterior22.csv');
a30 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Anterior30.csv');
b30 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Posterior30.csv');
a67 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Anterior67.csv');
b67 = csvread('/Users/hubatsl/Desktop/interplay-cell-size/MatsMets/Posterior67.csv');
DetSimus{1}{1} = a16;
DetSimus{1}{2} = b16;
DetSimus{2}{1} = a18;
DetSimus{2}{2} = b18;
DetSimus{3}{1} = a20;
DetSimus{3}{2} = b20;
DetSimus{4}{1} = a22;
DetSimus{4}{2} = b22;
DetSimus{5}{1} = a30;
DetSimus{5}{2} = b30;
DetSimus{6}{1} = a67;
DetSimus{6}{2} = b67;
s.average_datasets();
s.plot_Simu({'L', 'DA', 'DB', 'koffA', 'koffB'}, DetSimus);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 2*4.0 2*3.1];
% For presentation
allAxesInFigure = findall(gcf,'type','axes');
for i = 1:length(allAxesInFigure)
    allAxesInFigure(i).Color = 'black';
    allAxesInFigure(i).XColor = [247/255 214/255 31/255];
    allAxesInFigure(i).YColor = [247/255 214/255 31/255];
end

allAxesInFigure = findall(gcf,'type','line');
for i = 1:length(allAxesInFigure)
    allAxesInFigure(i).LineWidth=2;
end

allAxesInFigure = findall(gcf,'type','text')
for i = 1:length(allAxesInFigure)
    allAxesInFigure(i).Color=[247/255 214/255 31/255];
end

g = gcf;
% g.PaperPosition = [0 0 1.8*4.0 1.7*3.1];
g.InvertHardcopy = 'off';
g.Color = 'black';
saveas(gcf,'/Users/hubatsl/Desktop/VivaTalk/StochVsDet.png')
% print('/Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/StochVsDet', '-dpng', '-r600');
%% Compare noise in stochastic model and real data
% load 16_05_17_green_red
m = s.SimulationCell{end, 6}.B(:, end);
% m = m(5:end);
% m = (m-min(m))/(mean(m-min(m)));
m = (m-min(m))/(max(m)-min(m));
b = 5;
r = repmat(b, floor(length(m)/b), 1);
b = mat2cell(m, [r; length(m) - sum(r)], 1);
b = nanmean(cellfun(@std, b));
figure;
plot(m);
text(10, max(m)/2, ['\sigma_{mov} = ', num2str(b, 3)]);
make_figure_pretty([0 Inf 0 Inf], 'x [\mum]', 'particle # [a.u.]');
print('/Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/SimuNoise', '-dpng', '-r600');
%%
d = Mo_Segs_P0_15_12_05_red(7).Side1{end-10}.fitData;
n = floor(1/0.124);
r = repmat(n, floor(length(d)/n), 1);
n = mat2cell(d, 1, [r; length(d)-sum(r)]);
n = cellfun(@mean, n);
% n = (n-min(n))/(mean(n-min(n)));
n = (n-min(n))/(max(n)-min(n));
n = n(1:60);

b = 5;
r = repmat(b, floor(length(n)/b), 1);
b = mat2cell(n, 1, [r; length(n) - sum(r)]);
b = nanmean(cellfun(@std, b));

figure;
plot(n);
text(10, double(max(n)/2), ['\sigma_{mov} = ', num2str(b, 3)]);
make_figure_pretty([0 Inf 0 Inf], 'x [\mum]', 'particle # [a.u.]');
print('/Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/RealNoise', '-dpng', '-r600');
%% Kymograph for large and small cells
figure(1); hold on; % 18 microns
A = s.SimulationCell{2}.A;
B = s.SimulationCell{2}.B;
[X, Y] = meshgrid(linspace(0, 18, 18), linspace(0, 1000, 11));
I = cat(3,A'/max(max(A)),B'/max(max(B)),zeros(size(A')));
Z1 = I(2:end, :, 2);
Z2 = I(2:end, :, 1); 
% Blue layer: 
p1 = surf(X,Y,zeros(size(Z1)),'AlphaData',Z1, 'FaceAlpha','interp',...
          'FaceColor','blue', 'edgecolor','none');
p2 = surf(X,Y,zeros(size(Z2)),'AlphaData',Z2, 'FaceAlpha','interp',...
          'FaceColor','red', 'edgecolor','none');
axis([0, max(max(X)), 0, max(max(Y)), 0, 1]);
make_figure_pretty([0, max(max(X)), 0, max(max(Y)), 0, 1], 'x [\mum]','Time [s]');
view(2) 
print('/Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/StochKymo18', '-dpng', '-r600');

figure(2); hold on; %67.3 microns
A = s.SimulationCell{end}.A;
B = s.SimulationCell{end}.B;
[X, Y] = meshgrid(linspace(0, 67.3, 67), linspace(0, 1000, 11));
I = cat(3,A'/max(max(A)),B'/max(max(B)),zeros(size(A')));
Z1 = I(2:end, :, 2);
Z2 = I(2:end, :, 1); 
% Blue layer: 
p1 = surf(X,Y,zeros(size(Z1)),'AlphaData',Z1, 'FaceAlpha','interp',...
          'FaceColor','blue', 'edgecolor','none');
p2 = surf(X,Y,zeros(size(Z2)),'AlphaData',Z2, 'FaceAlpha','interp',...
          'FaceColor','red', 'edgecolor','none');
axis([0, max(max(X)), 0, max(max(Y)), 0, 1]);
make_figure_pretty([0, max(max(X)), 0, max(max(Y)), 0, 1], 'x [\mum]','Time [s]');
view(2) 
print('/Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/StochKymo67_3', '-dpng', '-r600');

figure(3); hold on; % 20 microns
A = s.SimulationCell{3}.A;
B = s.SimulationCell{3}.B;
[X, Y] = meshgrid(linspace(0, 20, 20), linspace(0, 1000, 11));
I = cat(3,A'/max(max(A)),B'/max(max(B)),zeros(size(A')));
Z1 = I(2:end, :, 2);
Z2 = I(2:end, :, 1); 
p1 = surf(X,Y,zeros(size(Z1)),'AlphaData',Z1, 'FaceAlpha','interp',...
          'FaceColor','blue', 'edgecolor','none');
p2 = surf(X,Y,zeros(size(Z2)),'AlphaData',Z2, 'FaceAlpha','interp',...
          'FaceColor','red', 'edgecolor','none');
axis([0, max(max(X)), 0, max(max(Y)), 0, 1]);
make_figure_pretty([0, max(max(X)), 0, max(max(Y)), 0, 1], 'x [\mum]','Time [s]');
view(2) 
print('/Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/StochKymo20', '-dpng', '-r600');

figure(4); hold on;% 30 microns
A = s.SimulationCell{5}.A;
B = s.SimulationCell{5}.B;
[X, Y] = meshgrid(linspace(0, 30, 30), linspace(0, 1000, 11));
I = cat(3,A'/max(max(A)),B'/max(max(B)),zeros(size(A')));
Z1 = I(2:end, :, 2);
Z2 = I(2:end, :, 1); 
p1 = surf(X,Y,zeros(size(Z1)),'AlphaData',Z1, 'FaceAlpha','interp',...
          'FaceColor','blue', 'edgecolor','none');
p2 = surf(X,Y,zeros(size(Z2)),'AlphaData',Z2, 'FaceAlpha','interp',...
          'FaceColor','red', 'edgecolor','none');
view(2) 
make_figure_pretty([0, max(max(X)), 0, max(max(Y)), 0, 1], 'x [\mum]','Time [s]');
print('/Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/StochKymo30', '-dpng', '-r600');
%%
% sizes = 30:10:150;
% s = Simulations({{{'DB', 'DA', 'koffB', 'koffA', 'L', 'psi', 'bins'}},...
%                 {{fDP2(2*sizes), linspace(0.5, 0.01, length(sizes))', fOP2(2*sizes), fOP6(2*sizes),...
%                   sizes,S_to_V('Circumference', {2*sizes, 15/27}), 2*round(sizes)}}}, 20);
% s.run_simulations();
% save trial_higherDifferenceDA1.mat
% % Keep koffA/koffB the same, but incorporate their experimental change into D
% sizes = 30:10:150;
% L2 = sqrt(fDP2(2*sizes)./fOP2(2*sizes));
% L6 = sqrt(fDP6(2*sizes)./fOP6(2*sizes));
% D2 = L2.^2*0.0073;
% D6 = L6.^2*0.0054;
%
% s = Simulations({{{'DB', 'DA', 'L', 'psi', 'bins'}},...
%                 {{D2, D6, sizes, S_to_V('Circumference',...
%                  {2*sizes, 15/27}), 2*round(sizes)}}}, 20);
% s.run_simulations();
% save trial_keepKoffChangeD1.mat
%%

%% Don't change system size, change one diffusion coefficient, then both
% sizes = 50;
% s = Simulations({{{'DA'}, {'L', 'psi', 'bins'}},...
%                 {{logspace(-3, 0, 10)'},...
%                  {sizes, S_to_V('Circumference', {2*sizes, 15/27}), 2*round(sizes)}}}, 20);
% s.run_simulations();
% save trial_higherDifferenceDA2.mat

% % Do the same for P, should be different
% s = Simulations({{{'DB'}, {'L', 'psi', 'bins'}},...
%                 {{logspace(-3, 0, 10)'},...
%                  {sizes, S_to_V('Circumference', {2*sizes, 15/27}), 2*round(sizes)}}}, 20);
% s.run_simulations();
% save trial_higherDifferenceDA3.mat
% 
% % Change both diffusion coefficients
% s = Simulations({{{'DA', 'DB'}, {'L', 'psi', 'bins'}},...
%                 {{logspace(-3, 0, 10)', linspace(0.7, 0.001, 12)'},...
%                  {sizes, S_to_V('Circumference', {2*sizes, 15/27}), 2*round(sizes)}}}, 20);
% s.run_simulations();
% save trial_higherDifferenceDA4.mat
% 
% % Change size and D
% sizes = 30:20:150;
% s = Simulations({{{'DA', 'DB'}, {'L', 'psi', 'bins'}},...
%                 {{logspace(-3, 0, 10)', linspace(0.7, 0.001, 12)'},...
%                  {sizes, S_to_V('Circumference', {2*sizes, 15/27}), 2*round(sizes)}}}, 20);
% s.run_simulations();
% save trial_higherDifferenceDA5.mat

% load trial_higherDifferenceDA3.mat;
% m_fit = eval_and_plot(s, 'DA', 0.5);
% eval_and_plot(s, 'DB', 0.5);

% load trial_higherDifferenceDA3.mat;
% m_fit = eval_and_plot(s, 'DA', 0.5);
m_fit = abs(eval_and_plot(s, 'DA', 0.5));

% function m_fit = eval_and_plot(s, species, pixelsize)
%     s.average_datasets()
%     m_all = zeros(1, length(s.AvSimus));
%     m_ind = zeros(1, length(s.AvSimus));
%     m_mp = zeros(1, length(s.AvSimus));
%     m_fit = zeros(1, length(s.AvSimus));
%     if species == 'DA'
%         n = 1;
%     elseif species == 'DB'
%         n = 2;
%     end
%     for i = 1:length(s.AvSimus)
%         int = s.AvSimus{i}{n}(:, end);
%         filt = sgolayfilt(int, 3, 31);
%         int = int/max(filt);
%         filt = filt/max(filt);
% %         filt(filt>0.9*max(filt)) = nan;
% %         filt(filt<0.1*max(filt)) = nan;
%         diff_filt = abs(diff(filt));
%         [m_all(i), m_ind(i)] = max(diff_filt);
%         % slope at midpoint:
%         [~, ind] = min(abs(filt-(max(filt)-min(filt))/2));
%         dist_side = round(5/pixelsize);
%         try 
%             m_mp(i) = diff_filt(ind);
%             fit_y = int(ind-dist_side:ind+dist_side);
%             fit_x = 0:pixelsize:2*(dist_side)*pixelsize;
%             f = fit(fit_x', fit_y, 'poly1');
%             m_fit(i) = f.p1;
%         catch
%             m_fit(i) = nan;
%         end
%     end
%     x = eval(['s.param_table_summary.', species]);
% %     f = fit(x(6:end), m_fit(6:end)', 'A./sqrt(x)');
%     figure(2); hold on;
% %     plot(x, m_all)
%     plot(x, abs(m_fit));
% %     plot(x, f(x));
% end
%%
% load s3
load s4
pixelsize = 0.5;
m_fit = zeros(size(s.SimulationCell));
Di = zeros(1, 19);

for j = 3:19
for i = 1:20
    int = s.SimulationCell{j, i}.B(:, end);
    filt = sgolayfilt(int, 3, 31);
    int1 = int/max(filt);
    int = s.SimulationCell{j, i}.A(:, end);
    filt = sgolayfilt(int, 3, 31);
    int2 = flipud(int/max(filt));
%     filt = filt/max(filt);
%     [~, ind] = min(abs(filt-(max(filt)-min(filt))/2));
%     dist_side = round(3/pixelsize);
%     fit_y = int(ind-dist_side:ind+dist_side);
%     fit_x = 0:pixelsize:2*(dist_side)*pixelsize;
%     f = fit(fit_x', fit_y, 'poly1');
%     m_fit(j, i) = f.p1;
    f = Fit(int1, 'off', 1, 'err');
    mB_fit(j, i) = f.curve.s;
    f = Fit(int2, 'off', 1, 'err');
    mA_fit(j, i) = f.curve.s;
end
end
%%
plot(s.param_table_summary.L, 1./abs(mean(m_fit, 2)))
% axis([0, 100, 0, 20])
% plot(s.param_table_summary.L, 1./abs(mean(m_fit, 2)))
%% 
figure; hold on;
x = 2*s.param_table_summary.L;
y1 = mean(mA_fit, 2);
y1std = std(mA_fit, 0, 2);
y2 = mean(mB_fit, 2);
y2std = std(mB_fit, 0, 2);
% errorbar(x(3:end), y1(3:end), y1std(3:end), 'r', 'LineWidth', 2);
plot(x(3:end), y1(3:end), 'r', 'LineWidth', 2);
plot(x(3:end), y1(3:end)+ y1std(3:end), 'r', 'LineWidth', 2);
plot(x(3:end), y1(3:end)- y1std(3:end), 'r', 'LineWidth', 2);
legend({'PAR-6 simulated'}, 'Location', 'South');
c = [76/255 134/255 198/255];
% c = 'b'; % For thesis/paper
% errorbar(x(3:end), y2(3:end), y2std(3:end), 'Color',   c, 'LineWidth', 2);
% plot(x(3:end), y2(3:end), 'Color', c, 'LineWidth', 2);
% plot(x(3:end), y2(3:end)+ y2std(3:end), 'Color', c, 'LineWidth', 2);
% plot(x(3:end), y2(3:end)- y2std(3:end), 'Color', c, 'LineWidth', 2);
% legend({'PAR-2 simulated'}, 'Location', 'South');
make_figure_pretty([0, 200, 0, 25], 'x [\mu m]', '\lambda [\mu m]');
% For presentation
% l = legend({'PAR-6 simulated', 'PAR-2 simulated'}, 'Location', 'SouthWest');
% l.FontSize=12;
% legend boxoff
% l.Color = [0 0 0];
% ax = gca;
% ax.Color = 'black';
% ax.XColor = [247/255 214/255 31/255];
% ax.YColor = [247/255 214/255 31/255];
% g = gcf;
% g.PaperPosition = [0 0 1.5*4.0 1.5*3.1];
% g.InvertHardcopy = 'off';
% g.Color = 'black';
% saveas(gcf,'/Users/hubatsl/Desktop/VivaTalk/SimusGradientP6.png')
% print('/Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/SimusGradientP6', '-dpng');
print(['/Users/lhcge/Dropbox (Lars DMS)/LRI/PAR_Size/PAR_Size_MS_Data/Data_Figure4_Kinetics/',...
        'SimusGradientP6'], '-depsc');
%%
figure; hold on;
x = 2*s.param_table_summary.L;
y1 = mean(mA_fit, 2);
y1std = std(mA_fit, 0, 2);
y2 = mean(mB_fit, 2);
y2std = std(mB_fit, 0, 2);
plot(x(3:end), y1(3:end),  'r', 'LineWidth', 2);
plot(x(3:end), y2(3:end),  'b', 'LineWidth', 2);
make_figure_pretty([0, 200, 0, 25], 'x [\mu m]', '\lambda [\mu m]');
size_range = 0:200;
p2 = sqrt(fDP2(size_range)./fOP2(size_range));
p6 = sqrt(fDP6(size_range)./fOP6(size_range));
% p6art = sqrt(linspace(0.5, 0.01, length(size_range))'./fOP6(size_range));
norm_factor = mean([y2(end)/p2(end), y1(end)/p6(end)]);
plot(size_range, y2(end)/p2(end)*p2);
plot(size_range, y2(end)/p2(end)*p6);
%%
s.plot_Simu({'DA', 'DB', 'koffA', 'koffB'});
% print('/Users/hubatsl/Desktop/interplay-cell-size/TheoryApplied/SimusSize', '-dpng');
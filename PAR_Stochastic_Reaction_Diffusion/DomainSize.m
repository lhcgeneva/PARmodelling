numSims = 16;
prec = 10;
%%
ratios = linspace(0.8,2.2,prec);
s = Simulations({{{'Bconc', 'ratioAB'}},...
                {{79./ratios, ratios}}}, numSims);
s.run_simulations();
%%
s.average_datasets;
s.plot_Simu({'ratioAB', 'Bconc'})
%%
pixelsize = 1;
m_fit = zeros(size(s.SimulationCell));
mA_fit = zeros(prec, numSims);
mB_fit = zeros(prec, numSims);
dsA_fit = zeros(prec, numSims);
dsB_fit = zeros(prec, numSims);
for j = 1:prec
for i = 1:numSims
    int = s.SimulationCell{j, i}.B(:, end);
    filt = sgolayfilt(int, 3, 7);
    int1 = int/max(filt);
    int = s.SimulationCell{j, i}.A(:, end);
    filt = sgolayfilt(int, 3, 7);
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
    dsB_fit(j, i) = f.curve.c;
    f = Fit(int2, 'off', 1, 'err');
    mA_fit(j, i) = f.curve.s;
    dsA_fit(j, i) = f.curve.c;
end
end
%%
figure;
yyaxis right; hold on;
y1std = std(mB_fit, 0, 2);
plot(cellfun(@(x) x.params.Bconc, s.SimulationCell(:, 1)), mean(mB_fit'), 'LineWidth', 2);
plot(cellfun(@(x) x.params.Bconc, s.SimulationCell(:, 1)), mean(mB_fit')+y1std', '-', 'LineWidth', 1);
plot(cellfun(@(x) x.params.Bconc, s.SimulationCell(:, 1)), mean(mB_fit')-y1std', '-', 'LineWidth', 1);
axis([0, 120, 0, 20]);
ylabel('\lambda [\mu m]');
yyaxis left; hold on;
y1std = std((120-dsB_fit)./120, 0, 2);
plot(cellfun(@(x) x.params.Bconc, s.SimulationCell(:, 1)), (120-mean(dsB_fit'))/120, 'LineWidth', 2);
plot(cellfun(@(x) x.params.Bconc, s.SimulationCell(:, 1)), (120-mean(dsB_fit'))/120+y1std', '-', 'LineWidth', 1);
plot(cellfun(@(x) x.params.Bconc, s.SimulationCell(:, 1)), (120-mean(dsB_fit'))/120-y1std', '-', 'LineWidth', 1);
axis([0, 120, 0.1, 0.7]);
make_figure_pretty([0, 120, 0.1, 0.7], 'C_P [a.u.]', 'Domain Size');
legend off
print('/Users/lhcge/Dropbox (Lars DMS)/LRI/PAR_Size/PAR_Size_MS_Data/Data_Figure3_Theory2/SimusSize', '-dpng');
savefig('/Users/lhcge/Dropbox (Lars DMS)/LRI/PAR_Size/PAR_Size_MS_Data/Data_Figure3_Theory2/SimusSize');
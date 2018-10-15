% All artefact of cell size? This is no cell size change, only D changes:
load trial_higherDifferenceDA3.mat;
m_fit = zeros(size(s.SimulationCell));
Di = zeros(1, 10);

for j = 2:10
    Di(j) = s.SimulationCell{j, 1}.params.DB;
for i = 1:20
    int = s.SimulationCell{j, i}.B(:, end);
    filt = sgolayfilt(int, 3, 31);
    int = int/max(filt);
    filt = filt/max(filt);
    [~, ind] = min(abs(filt-(max(filt)-min(filt))/2));
    dist_side = round(4/pixelsize);
    try 
        fit_y = int(ind-dist_side:ind+dist_side);
        fit_x = 0:pixelsize:2*(dist_side)*pixelsize;
        f = fit(fit_x', fit_y, 'poly1');
        m_fit(j, i) = f.p1;
    catch
        m_fit(j, i) = nan;
    end
end
end

plot(Di, abs(mean(m_fit, 2)))
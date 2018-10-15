% clear all
% load tempproper.mat
s.average_datasets();
m = cellfun(@(x) [max(abs(diff(sgolayfilt(x{1}(:, end),3, 11)))),...
                  max(abs(diff(sgolayfilt(x{2}(:, end),3, 11))))],...
                  s.AvSimus, 'UniformOutput', false);
m = mean(cell2mat(m), 2);
plot(s.param_table_summary.L, m)
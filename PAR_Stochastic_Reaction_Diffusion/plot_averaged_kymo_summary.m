% s is output of param_script, cell containing simulations. script plots
% average of all simulations for each parameter set investigated.
dimensions = size(s{1});
for j = 1:dimensions(2)
    figure(j)
    for i = 1:length(s)
        A_t = [s{i}];
        A_tt = {A_t.A};
        A = mean(cat(3, A_tt{:}), 3);
        B_t = [s{i}];
        B_tt = {B_t.B};
        B = mean(cat(3, B_tt{:}), 3);
        I = cat(3,A'/max(max(A)),B'/max(max(B)),zeros(size(A')));
        I(:, :, 3) = 0.5;
        subplot(3, 5, i);
        c = imfuse(A.', B.','falsecolor','Scaling','joint','ColorChannels','red-cyan');
        imshow(c);
        axis square;
    end
end
%%
% for i = 1:length(S.SimulationCell)
    s = S.SimulationCell{1}
    A_tt = {s.A}
    B_tt = {s.B};
    A = mean(cat(3, A_tt{:}), 3);
    B = mean(cat(3, B_tt{:}), 3);
    I = cat(3,A'/max(max(A)),B'/max(max(B)),zeros(size(A')));
    I(:, :, 3) = 0.5;
    figure
    c = imfuse(A.', B.','falsecolor','Scaling','joint','ColorChannels','red-cyan')
    imshow(c);
    axis square;
    tightfig;
% end
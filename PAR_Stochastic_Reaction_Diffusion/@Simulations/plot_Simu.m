function plot_Simu( Simu, show_params, DetSimus )
%PLOT_SIMU Plots all Simulations (last timepoint) using tight subplots
%   Also puts the parameters for each performed simulation on the plot that
%   are contained in show_params (e.g. show_params = {'kOn', 'kOff',
%   'kcB'})

numSubPlots = length(Simu.AvSimus);
dimNumR = ceil(sqrt(numSubPlots));
if dimNumR * (dimNumR - 1) >= numSubPlots
    dimNumC = dimNumR -1;
else
    dimNumC = dimNumR;
end
% dimNumR = dimNumR-1; % temporarily, to only have 16 graphs
figure;
ts =  tight_subplot(dimNumR, dimNumC, 0.08, 0.04, 0.04);
counter = 1;
for i = 1 : dimNumR
    for j = 1 : dimNumC
        axes(ts(counter));
        hold on;
        l = length(Simu.AvSimus{counter}{2}(:, end));
        x = linspace(0.5, l - 0.5, l);
        plot(x, Simu.AvSimus{counter}{2}(:, end));
        plot(x, Simu.AvSimus{counter}{1}(:, end));
        if nargin > 2
            no = mean(Simu.AvSimus{counter}{2}(end-4:end, end));
            A = DetSimus{counter}{1}(:, end);
            B = DetSimus{counter}{2}(:, end);
            m = max(B); % define separately to be able to reuse in B=no/m*B
            A = no/m*A;
            B = no/m*B;
            x = linspace(0, length(Simu.AvSimus{counter}{1}(:, end)), length(B));
            plot(x, A, 'r--');
            plot(x, B, 'b--');
        end
        %Show legend only in first graph
%         if counter == 1
%             legend('Species A', 'Species B', 'Location', 'BEST');
%         end
%         legend('boxoff');
        set(gca, 'FontName', 'Palatino');
        set(gca,'TickDir','out');
        set(gca, 'FontSize', 9, 'LineWidth', 1)
        box off;
        if nargin > 1
            %Get text position, adjust axes
            m = max(max([Simu.AvSimus{counter}{1}(:, end), ...
                    Simu.AvSimus{counter}{2}(:, end)]));%max of A & P
            l = length(Simu.AvSimus{counter}{1}(:, end));%numel in x
            h = m+0.1*m*(length(show_params)+1); %Total height of graph
            axis([0 l 0 h]);
            for ii = 1 : length(show_params)
                text(l/10, h-0.11*m*ii, [show_params{ii}, ' = ', ...
                    num2str(Simu.param_table_summary.(show_params{ii})(counter), 3)],...
                    'FontName', 'Palatino', 'FontSize', 8);
%                 text(10, 10,'Check whether A and B are in correct order!!!');
            end
        end
        counter = counter + 1;
        if counter > numSubPlots; break; end;
    end
end
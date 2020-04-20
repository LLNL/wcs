% Source: Matlab user forum discussion
% https://www.mathworks.com/matlabcentral/answers/224604-implementing-gillespie-s-algorithm
% by Harley Day on 11 Feb 2019
lineStyles = ['-', '--', ':', '-.'];
markers = ['+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'];
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];

seeds = [59 67 73 83 97 103 109 127 137 149 157];
path = 'result-trace/toplot/';

Y0 = 1000;
ssa = 0;
cms_idx = 1; % colori-marker-style index
FigW=10;
FigH=8;

for s=1:length(seeds)
    close all;
    seed = seeds(s);
    initPop = [ Y0 Y0 ];
    filename = [ path 'eq38.Y-' num2str(Y0) '.m-' num2str(ssa) '.seed-' num2str(seed) ];
    
    tY = load([ filename '.txt' ]);
    t = tY(:,1);
    Y1 = tY(:,2);
    Y2 = tY(:,3);
    clear tY;
    
    h1 = figure;
%    clf;
%    set(h1,'defaulttextinterpreter','latex',...
%                'PaperUnits','centimeters','PaperSize',[FigW FigH],...
%                'PaperPosition',[0,0,FigW,FigH],'Units','centimeters',...
%                'Position',[1,01,FigW,FigH]);

    hold on;
    c1 = colors(mod(cms_idx-1, length(colors))+1);
    c2 = colors(mod(cms_idx, length(colors))+1);
    h1 = stairs(t, Y1, ...
             'LineWidth', 0.4, ...
             'Color', c1, ...
             'MarkerSize', 0.5, ...
             'Marker', markers(mod(cms_idx-1, length(markers))+1), ...
             'MarkerFaceColor', c1, ...
             'MarkerEdgeColor', c1);
    
    h1 = stairs(t, Y2, ...
             'LineWidth', 0.4, ...
             'Color', c2, ...
             'MarkerSize', 0.5, ...
             'Marker', markers(mod(cms_idx, length(markers))+1), ...
             'MarkerFaceColor', c2, ...
             'MarkerEdgeColor', c2);
    grid on;
    hold off;
    xlim([0 30]);
    xlabel('Minutes');
    %ylabel('Number of Y molecules', 'interpreter', 'latex');
    ylabel('Number of molecules');
    
    lgn1{1} = ['$Y1, Y1_0$ = [' num2str(initPop(1)) ']'];
    lgn1{2} = ['$Y2, Y2_0$ = [' num2str(initPop(2)) ']'];
    legend(lgn1, 'interpreter', 'latex', 'Location', 'NorthWest');
    
    print([filename '-t-Y1Y2.pdf'], '-dpdf');
    
    h2 = figure;
%    clf;
%    set(h2,'defaulttextinterpreter','latex',...
%                'PaperUnits','centimeters','PaperSize',[FigW FigH],...
%                'PaperPosition',[0,0,FigW,FigH],'Units','centimeters',...
%                'Position',[1,10,FigW,FigH]);
    hold on;
    h2 = plot(Y1, Y2, ...
             'LineWidth', 0.4, ...
             'Color', c1, ...
             'MarkerSize', 0.5, ...
             'Marker', markers(mod(cms_idx-1, length(markers))+1), ...
             'MarkerFaceColor', c1, ...
             'MarkerEdgeColor', c1);
    grid on;
    hold off;
    xlabel('$Y1$', 'interpreter', 'latex');
    ylabel('$Y2$', 'interpreter', 'latex');

    xl = xlim;
    yl = ylim;
    xlim([min([xl yl]) max([xl yl])]);
    ylim([min([xl yl]) max([xl yl])]);
    
    lgn2{1} = ['$Y_0$ = [' num2str(initPop(1)) ' ' num2str(initPop(2)) ']'];
    legend(lgn2(1), 'interpreter', 'latex' );
    
    print([filename '-Y1-Y2.pdf'], '-dpdf');
end

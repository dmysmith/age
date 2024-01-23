%% Define function for plotting
function plot_rois(volmat, x, groups, ylab)
    axes(figure('Position', [10 10 900 600]));
    hold on;
    handles = {};

    cats = unique(groups');
    
    for g = 1:length(cats)
        group = cats{g};
        idx = (strcmp(groups, group));
        % xplot = x(idx,:);
        % splot = volmat(idx,:);
        xplot = x;
        splot = mean(volmat(idx,:),1);
        ci = [splot + (2*std(volmat(idx,:),1)) flip(splot - (2*std(volmat(idx,:),1)))];
        handles{end+1} = plot(xplot(:), splot(:),'DisplayName', group, 'LineWidth',2);
    end
    legend('Location','northwest','Interpreter', 'none');
    % set_colors;
    % colororder(mycolors_7);
    xlabel('Age (Years)'); 
    ylabel(ylab);
end
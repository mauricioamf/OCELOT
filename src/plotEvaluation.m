function plotEvaluation(metrics_cutoff, selectedWTReference, showPlots)
    if isequal(showPlots, 0) || isequal(string(showPlots), "0") || isempty(metrics_cutoff)
        return;
    end

    selectedReferenceIdx = find(contains(metrics_cutoff.Condition, selectedWTReference), 1, 'first');
    if (isequal(showPlots, 1) || strcmpi(string(showPlots), 'cutoff')) && ~isempty(selectedReferenceIdx)
        figure;
        tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

        nexttile;
        hold on;
        plot(metrics_cutoff.fpRate{selectedReferenceIdx,:}, ...
             metrics_cutoff.tpRate{selectedReferenceIdx,:}, ...
             'LineWidth', 1.5, ...
             'DisplayName', sprintf('Model (AUC = %.3f, n_{ess}=%d)', ...
             metrics_cutoff.AUC(selectedReferenceIdx), ...
             metrics_cutoff.n_ess_covered(selectedReferenceIdx)));
        plot([0,1],[0,1],'--','Color',[0.5 0.5 0.5],'DisplayName','Random (AUC = 0.5)');
        hold off;
        xlabel('False Positive Rate'); ylabel('True Positive Rate');
        legend('Location','southeast'); grid on;

        nexttile;
        confusionchart(metrics_cutoff.C{selectedReferenceIdx,:}, {'Non-essential','Essential'});
    end

    if isequal(showPlots, 2) || strcmpi(string(showPlots), 'all')
        figure;
        tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

        nexttile; histogram(metrics_cutoff.AUC); xlabel('AUC');
        nexttile; histogram(metrics_cutoff.Accuracy); xlabel('Accuracy');
        nexttile; histogram(metrics_cutoff.Sensitivity); xlabel('Sensitivity');
        nexttile; histogram(metrics_cutoff.Specificity); xlabel('Specificity');
    end
end
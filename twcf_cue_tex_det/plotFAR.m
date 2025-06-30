
p = twcf_analysisParams; 
cNames = {'sawFigure','sawShape','sawFigureNS'}; 

s.saveFigs = 1; 
s.figType = 'png'; 

%% Compile FAR by line length 
for iC = 1:numel(cNames)
    for iS = 1:30
        FARs.(cNames{iC})(iS,:,:) = a(iS).sdt.(cNames{iC}).det_pfa;
    end
end

%% Plot FAR by line length
cNamesLabel = {'saw figure','saw shape'};
for iC = 1:2
    figure
    hold on
    figureStyle

    ylabel(sprintf('False alarm rate\n%s',cNamesLabel{iC}))
    xlabel('Line length')
    xCorrect = [-0.2 0 0.2];
    ax = gca;
    ax.FontSize = 18;
    xlim([0.5 7.5])

    for iL = 1:7
        for iAtt = 1:3
            y = FARs.(cNames{iC})(:,iL,iAtt);
            bar(iL+xCorrect(iAtt),mean(y),0.2, 'FaceColor',p.style.attColorsMuted(iAtt,:),...
                'EdgeColor','none','FaceAlpha',1)
            errb = errorbar(iL+xCorrect(iAtt), mean(y), std(y)/sqrt(size(y,1))  );
            errb.LineWidth = 2;
            errb.Color = 'k';
            errb.CapSize = 0;
        end
    end

    if s.saveFigs
        figTitle = sprintf('group_FARs_C%d',iC);
        saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, s.figType))
    end
end



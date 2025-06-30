%% For email with Jan Brascamp April 2025


%% make pds
mus = [-1.5 1.5];
muScaling = 1.2; 
sigma = 1.5; 
sigmaScaling = 0.6;
for iX = 1:2
    for iA = 1:2
        switch iA
            case 1
                pd_signal(iX,iA) = makedist('Normal','mu',mus(iX)*muScaling,'sigma',sigma*sigmaScaling); 
            case 2
                pd_signal(iX,iA) = makedist('Normal','mu',mus(iX),'sigma',sigma); 
        end
    end
end


%% plot
figure
set(gcf,'Position',[100 100 500 200])
figureStyle
for iX = 1:2
    for iA = 1:2 % 1 attended, 2 unattended
        pd = pd_signal(iX,iA);
        pl(iX,iA) = plot(pd);
        switch iA
            case 1
                pl(iX,iA).LineStyle = '-'; 
            case 2
                pl(iX,iA).LineStyle = '--'; 
        end
  
    switch iX
        case 1
            pl(iX,iA).Color = [57 106 177]/255; 
        case 2
            pl(iX,iA).Color = [62,150,81]/255; 
    end
    end
end
xlim([-6 6])
ylim([0 0.5])
xline(0,'-k','saw stimulus','LabelOrientation','horizontal')
xlabel('Internal signal strength')
xticks('')
yticks('')
legend([pl(1,1) pl(1,2)],{'valid','invalid'},'FontSize',9)





function kt_SDTfigs_criterion

%% Setup 
p = twcf_analysisParams; 

%% Make probability distributions
% Matched stimulus strengths
criterion = 0.9; 
mus = [0 0.8 2.5];
sigma = 1; 

for iX = 1:numel(mus)
    val_mu = mus(iX);
    val_sigma = sigma;
    pd_signal(iX) = makedist('Normal','mu',val_mu,'sigma',val_sigma);
    pd_cdf(iX) = 1-sum(cdf('Normal',criterion,val_mu,val_sigma));
end

sds = -2:2;
for iX = 1:numel(mus)
    for iS = 1:numel(sds)
        c(iX,iS) = mus(iX)+sds(iS);
    end
end

%% plot PDs
figure
set(gcf,'Position',[100 100 600 200])
figureStyle
iA = 1; 
pColors = [0.5 0.5 0.5; 
    p.style.attColorsMuted(1,:);
    p.style.attColorsMuted(3,:)]; 
for iX = 1:3
    pd = pd_signal(iX);
    pl(iX) = plot(pd);
    pl(iX).LineWidth = 1.5;

    xIdx = find(pl(iX,iA).XData>criterion);
    pa(iX) = area(pl(iX).XData(xIdx),pl(iX).YData(xIdx));
    pa(iX).FaceAlpha = 0.4;
    pa(iX).EdgeColor = 'none';

    switch iX
        case 1
            pl(iX,iA).LineStyle = '--';
        case 2
            pl(iX,iA).LineStyle = '-';
    end
    pl(iX).Color = pColors(iX,:);
    pa(iX).FaceColor = pColors(iX,:);

end

xticks(c(2,:))
xticklabels(sds)

ax = gca; 
ax.YAxis.Visible = 'off'; % remove y-axis
xlim([-3 6])
ylim([0 0.5]); 
xline(criterion,'-k','"saw stimulus"','LabelOrientation','horizontal')
xlabel('Internal signal strength','FontSize',14)

yticks('')
ylabel([])
legend([pl(1,1) pl(1,2)],{'valid','invalid'},'FontSize',9,'Box','off')






%% How do multiplicative vs. additive effects of attention
% differentially map onto d'

%% Simulate d'
means = linspace(0,5,20);

attScaling = 1:0.1:5;

for iM = 1:numel(means)
    for iA = 1:numel(attScaling)
        d_noAtt(iM,iA) = means(iM); 
        d(iM,iA) = means(iM)*attScaling(iA); % mulitplicative model 
        d_add(iM,iA) = means(iM)+attScaling(iA); % additive model 
    end
end


%% 
cmap = cmocean('balance');
cmap_red = cmocean('amp');

figure
set(gcf,'Position',[0 0 350 500])
set(0, 'DefaultFigureRenderer', 'painters');
nSubplots = 4;

ax(1) = subplot (nSubplots,1,1); 
hold on 
imagesc(attScaling,means,d_noAtt)
colormap(ax(1),cmap_red)
c = colorbar;
c.Label.String = sprintf('\\itd'''); 
title('No attention')
xticks([])
yticks([0:5])
xlim([1 5])
ylim([0 5])

ax(2) = subplot (nSubplots,1,2);
hold on 
imagesc(attScaling,means,d-d_noAtt)
colormap(ax(2),cmap_red)
c = colorbar;
c.Label.String = sprintf('\\itd'''); 
title('Attention mulitiplicatively scales mean')
xticks(1:5)
yticks(0:5)
xlim([1 5])
ylim([0 5])

ax(3) = subplot (nSubplots,1,3);
hold on 
imagesc(attScaling,means,d_add)
colormap(ax(3),cmap_red)
c = colorbar;
c.Label.String = sprintf('\\itd'''); 
title('Attention additively scales mean')
xticks(1:5)
yticks(0:5)
xlim([1 5])
ylim([0 5])

ax(4) = subplot (nSubplots,1,4);
hold on 
imagesc(attScaling,means,d-d_add)
colormap(ax(4),cmap)
c = colorbar; 
c.Label.String = sprintf('\\Delta \\itd'); 
% Get color axis limits
caxis_limits = clim;
% Adjust color axis limits
max_limit = max(abs(caxis_limits));
clim([-max_limit max_limit]);

xlabel(sprintf('\\alpha Attentional scaling'))
title('Difference between multiplicative and additive scaling')
xticks(1:5)
yticks(0:5)
xlim([1 5])
ylim([0 5])
ylabel(sprintf('\\mu_{unattended}'))

%% Let's just take a 1d slice 

attScalingIdx = 6; % 1.5 
% meansIdx = ; % 1.8421 all means

figure

ax(1) = subplot (3,1,1);
hold on 
plot(means,d(:,attScalingIdx))
plot(means,d_add(:,attScalingIdx))
xlabel('Mean')
ylabel('d''')
title('Multiplicative mean scaling')

ax(2) = subplot (3,1,2);
hold on 
plot(means,d_add(:,attScalingIdx))
xlabel('Mean')
ylabel('d''')
title('Additive mean scaling')

ax(3) = subplot (3,1,3);
hold on 
plot(means,d_add(:,attScalingIdx)-d(:,attScalingIdx))
xlabel('Mean')
ylabel('d''')
title('Multiplicative - additive')
yline(0)


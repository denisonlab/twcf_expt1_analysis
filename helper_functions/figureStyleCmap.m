function figureStyleCmap()
% TWCF FOHO figure styling 
% adjusts fig axis and text styling
% Imagesc heatmap figure styling

hold on 
box off
set(gca,'TickDir','in');
set(gca, 'Layer', 'Top');
set(gca, 'Color', 'w');
% set(gca, 'color', 'none')
ax = gca;
ax.LineWidth = 1;
ax.XColor = 'black';
ax.YColor = 'black';
% ax.FontSize = 18; % 12
% ax.FontWeight = 'light'; 

smlFont = 12;
bigFont = 16; % 24 

ax.FontSize = bigFont;
ax.FontName = 'Helvetica'; 
 
ax.XAxis.FontSize = smlFont;
ax.YAxis.FontSize = smlFont;

ax.XLabel.FontSize = bigFont;
ax.YLabel.FontSize = bigFont;
function figureStyle()
% TWCF FOHO figure styling 
% adjusts fig axis and text styling

hold on 
box off
set(gca,'TickDir','out');
set(gca, 'Layer', 'Top');
set(gca, 'Color', 'w');
% set(gca, 'color', 'none')
ax = gca;
ax.LineWidth = 1.5; % 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
% ax.FontSize = 18; % 12
% ax.FontWeight = 'light'; 

smlFont = 14;
bigFont = 18; % 24 

ax.FontSize = bigFont;
ax.FontName = 'Helvetica'; 
 
ax.XAxis.FontSize = smlFont;
ax.YAxis.FontSize = smlFont;

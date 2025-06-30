function twcf_drawBracket(x1,x2,y)
% Draws bracket indicating statistical test

yl = ylim;
xl = xlim;

y = y + 0.05; % correction to better align with significance indicator 
bracketHeight = (yl(2)-yl(1))/45; % 0.007 
bracketColor = 'black'; 

lineWidth = 1; 

% x = [x1 x1;... % vertical left 
%     x1 x2;... % horizontal 
%     x2 x2]; % vertical right 
% y = [max(yl)*y-bracketHeight max(yl)*y;...
%     max(yl)*y max(yl)*y;... 
%     max(yl)*y max(yl)*y-bracketHeight]; 


% 4 point line 
x = [x1; x1; x2; x2]; 
y = [max(yl)*y-bracketHeight; max(yl)*y; max(yl)*y; max(yl)*y-bracketHeight]; 


plot(x,y,'LineWidth',lineWidth,'Color',bracketColor,'LineJoin','miter')
 

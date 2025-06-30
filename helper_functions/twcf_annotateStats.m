function txt = twcf_annotateStats(x,y,stars)
% Draws n.s. and significance stars
% Inputs: 
%   x: x position
%   y: y position, max(fh.YLim) for top of figure 
%   stars: string, ns, *, **, *** 
% Ouputs: 
%   txt: text handle 

if contains(stars,'*')
    txt = text(x,y,stars,'EdgeColor','none',...
    'FontSize',20,'HorizontalAlignment','center','VerticalAlignment','Bottom','FontName','Helvetica'); % Times, make bigger 
elseif strcmp(stars,'ns')
    txt = text(x,y,'n.s.','EdgeColor','none',...
        'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','Bottom','FontName','Helvetica');
else % any alternative text 
    txt = text(x,y,stars,'EdgeColor','none',...
        'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','Bottom','FontName','Helvetica');
end



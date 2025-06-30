function twcf_annotateN(a)

xl = xlim;
yl = ylim;
nStr = sprintf('n = %d',size(a,2));
nStrTxt = text(0.945*xl(2),yl(1)+0.01,nStr,'HorizontalAlignment','right','VerticalAlignment','bottom');
nStrTxt.FontSize = 14;
nStrTxt.FontName = 'Helvetica-Light';



% Compile line lengths 

for i = 1:15
    lineLength_inDeg_postcue(i,:) = unique(dataAll(i).lineLength_inDeg_postcue);
end

figure
figureStyle
hold on 
for iS = 1:15
    for iL = 1:7
        scatter(iL,lineLength_inDeg_postcue(iS,iL))
    end
end

xlabel('r')
ylabel('Line length (deg)')
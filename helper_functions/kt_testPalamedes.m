
StimLevels = [.01 .03 .05 .07 .09 .11]; 
NumPos = [59 53 68 83 92 99];
OutOfNum = [100 100 100 100 100 100]; 

% psychometric function type to fit 
PF = @PAL_Logistic; 

paramsFree = [1 1 0 0]; % alpha beta gamma lambda 1-->free, 0-->fixed 
searchGrid.alpha = [0.01:0.001:0.11]; % limits shoudl be lowest stim intensity and highest stim intensity
searchGrid.beta = logspace(0,3,101); % slope range, throw wide net 
searchGrid.gamma = 0.5; 
searchGrid.lambda = 0.02; 

[paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum,searchGrid, paramsFree,PF);


figure
PropCorrectData = NumPos./OutOfNum;
StimLevelsFine = [min(StimLevels):(max(StimLevels)-min(StimLevels))./1000:max(StimLevels)];
Fit = PF(paramsValues, StimLevelsFine);
plot(StimLevelsFine,Fit,'g-','linewidth',2);
hold on;
plot(StimLevels, PropCorrectData,'k.','markersize',40);
set(gca, 'fontsize',12);
axis([0 .12 .4 1]);


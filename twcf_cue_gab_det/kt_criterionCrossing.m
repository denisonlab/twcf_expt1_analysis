

%% 

% mu = [0,0.1,0.2,0.4,0.6,1,2,4];
mu = [0:0.1:4];
sigma_att = 1; % 0.75
sigma_unatt = 1; 

k = 1.5;  % 1.5
mu_attScaling = 1.25; % 1 for no change, 1.25, 2

for i = 1:size(mu,2)
    dprime_att(i) = (mu(i)*mu_attScaling)/sigma_att; 
    dprime_unatt(i) = mu(i)/sigma_unatt; 

    c_att(i) = (k-(0.5*mu(i)*mu_attScaling))/sigma_att; 
    c_unatt(i) = (k-(0.5*mu(i)))/sigma_unatt; 
end


%% plot
idx = 1:size(mu,2); 

figure
set(gcf,'Position',[100 100 200 250])
subplot 211
figureStyle
hold on 
plot(mu(idx),dprime_att(idx),'-')
plot(mu(idx),dprime_unatt(idx),'-')
ylabel('{\itd''}')
xticks()
legend('att','unatt','FontSize',8,'Location','northwest','Box','off')

subplot 212
figureStyle
hold on 
plot(mu(idx),c_att(idx),'-')
plot(mu(idx),c_unatt(idx),'-')
ylabel('Criterion')

yline(0,'-k')
yline(1.25,'--k','1.25','LabelOrientation','horizontal')

sgtitle(sprintf('sigma scaling = %0.2f\nmean scaling = %0.2f',sigma_att,mu_attScaling),'FontSize',10,'FontWeight','normal')

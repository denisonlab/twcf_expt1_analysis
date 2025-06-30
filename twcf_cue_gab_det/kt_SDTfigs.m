function kt_SDTfigs 

%% If d' is mu/sigma, assuming mu_noise = 0. And attention reduces sigma by alpha. 
% d'_attended = mu/(sigma*alpha); 
% d'_unattended = mu/sigma; 
 
% In order to match performance, then 
% d'_unattended = mu*beta / sigma, where beta = 1/alpha

% In order to match visibility...

%% Setup 
p = twcf_analysisParams; 

%% Make probability distributions
% Matched stimulus strengths
criterion = 0.8; 
mus = [0 1 3/2];
sigma = 1; 
sigmaScaling = 2/3;
muScaling = 3/2; 
for iX = 1:numel(mus)
    for iA = 1:3 % attended, unattended
        switch iA
            case 1 % variance reduction and mean scaling
                val_mu = mus(iX)*muScaling; 
                val_sigma = sigma*sigmaScaling; 
            case 2 
                val_mu = mus(iX); 
                val_sigma = sigma;
        end
        pd_signal(iX,iA) = makedist('Normal','mu',val_mu,'sigma',val_sigma); 
        pd_cdf(iX,iA) = 1-sum(cdf('Normal',criterion,val_mu,val_sigma));
    end
end

% matched stimulus strength
figure
set(gcf,'Position',[100 100 500 200])
figureStyle
for iX = 1:2
    for iA = 1:2 % 1 attended, 2 unattended
        pd = pd_signal(iX,iA);
        pl(iX,iA) = plot(pd);
        pl(iX,iA).LineWidth = 1.5;

        xIdx = find(pl(iX,iA).XData>criterion);
        pa(iX,iA) = area(pl(iX,iA).XData(xIdx),pl(iX,iA).YData(xIdx));
        pa(iX,iA).FaceAlpha = 0.4;
        pa(iX,iA).EdgeColor = 'none';

        switch iX
            case 1
                pl(iX,iA).LineStyle = '--';
            case 2
                pl(iX,iA).LineStyle = '-';
        end

        switch iA
            case 1
                pl(iX,iA).Color = p.style.attColorsMuted(3,:);
                pa(iX,iA).FaceColor = p.style.attColorsMuted(3,:); 
            case 2
                pl(iX,iA).Color = p.style.attColorsMuted(1,:);
                pa(iX,iA).FaceColor = p.style.attColorsMuted(1,:); 
        end
    end
end
ax = gca; 
ax.YAxis.Visible = 'off'; % remove y-axis
xlim([-3 4])
ylim([0 0.7]); 
xline(0.8,'-k','saw stimulus','LabelOrientation','horizontal')
xlabel('Internal signal strength','FontSize',14)
xticks('')
yticks('')
ylabel([])
legend([pl(1,1) pl(1,2)],{'valid','invalid'},'FontSize',9,'Box','off')
text(-3,0.67,'At matched stimulus strength')

% --- Bar ---
figure
set(gcf,'Position',[100 100 300 200])
subplot 121 
figureStyle
hold on
d_att = (1*muScaling)/(sigma*sigmaScaling); 
d_unatt = (1)/(sigma); 
b(1) = bar(1,d_att);
b(2) = bar(2,d_unatt); 
b(1).FaceColor = p.style.attColorsMuted(3,:); 
b(1).EdgeColor = 'none'; 
b(2).FaceColor = p.style.attColorsMuted(1,:); 
b(2).EdgeColor = 'none'; 
ylabel('{\it d''}')
ylim([0 2.5])
yticks([0 2.5])
xlim([0.3 2.7])
xticks([1 2])
xticklabels({'Attended','Unattended'})

subplot 122 
figureStyle
hold on 
pd_attended_total = pd_cdf(1,1)+pd_cdf(2,1); 
pd_unattended_total = pd_cdf(1,2)+pd_cdf(2,2); 
b(1) = bar(1,pd_attended_total/2);
b(2) = bar(2,pd_unattended_total/2); 
b(1).FaceColor = p.style.attColorsMuted(3,:); 
b(1).EdgeColor = 'none'; 
b(2).FaceColor = p.style.attColorsMuted(1,:); 
b(2).EdgeColor = 'none'; 
ylabel('p("saw stimulus")')
ylim([0 1])
yticks([0 1])
xlim([0.3 2.7])
xticks([1 2])
xticklabels({'Attended','Unattended'})

%% matched performance
criterion = 0.8; 
mus = [0 1 3/2];
sigma = 1; 
sigmaScaling = 2/3;
muScaling = 3/2; 

d_att = (1/muScaling)/(sigma*sigmaScaling); 
d_unatt = (1)/(sigma); 

for iX = 1:numel(mus)
    for iA = 1:3 % attended, unattended
        switch iA
            case 1 % sigma reduction and mean reduction
                val_mu = mus(iX)/muScaling;
                val_sigma = sigma*sigmaScaling;
            case 2 
                val_mu = mus(iX);
                val_sigma = sigma;
        end
        pd_signal(iX,iA) = makedist('Normal','mu',val_mu,'sigma',val_sigma);
        pd_cdf(iX,iA) = 1-sum(cdf('Normal',criterion,val_mu,val_sigma));
    end
end

figure
set(gcf,'Position',[100 100 500 200])
figureStyle
for iX = [1 2]
    for iA = 1:2 % 1 attended, 2 unattended
        pd = pd_signal(iX,iA);
        pl(iX,iA) = plot(pd);
        pl(iX,iA).LineWidth = 1.5;

        xIdx = find(pl(iX,iA).XData>criterion);
        pa(iX,iA) = area(pl(iX,iA).XData(xIdx),pl(iX,iA).YData(xIdx));
        pa(iX,iA).FaceAlpha = 0.4;
        pa(iX,iA).EdgeColor = 'none';

        switch iX
            case 1
                pl(iX,iA).LineStyle = '--';
            case 2
                pl(iX,iA).LineStyle = '-';
        end

        switch iA
            case 1
                pl(iX,iA).Color = p.style.attColorsMuted(3,:);
                pa(iX,iA).FaceColor = p.style.attColorsMuted(3,:); 
            case 2
                pl(iX,iA).Color = p.style.attColorsMuted(1,:);
                pa(iX,iA).FaceColor = p.style.attColorsMuted(1,:); 
        end
    end
end
ax = gca; 
ax.YAxis.Visible = 'off'; % remove y-axis
xlim([-3 4])
ylim([0 0.7])
xline(0.8,'-k','saw stimulus','LabelOrientation','horizontal')
xlabel('Internal signal strength','FontSize',14)
xticks('')
yticks('')
ylabel([])
legend([pl(1,1) pl(1,2)],{'valid','invalid'},'FontSize',9,'Box','off')
text(-3,0.67,'At matched performance')

% --- Bar ---
figure
set(gcf,'Position',[100 100 300 200])
subplot 121 
figureStyle
hold on 
b(1) = bar(1,d_att);
b(2) = bar(2,d_unatt); 
b(1).FaceColor = p.style.attColorsMuted(3,:); 
b(1).EdgeColor = 'none'; 
b(2).FaceColor = p.style.attColorsMuted(1,:); 
b(2).EdgeColor = 'none'; 
ylabel('{\it d''}')
ylim([0 2.5])
yticks([0 2.5])
xlim([0.3 2.7])
xticks([1 2])
xticklabels({'Attended','Unattended'})

subplot 122 
figureStyle
hold on 
pd_attended_total = pd_cdf(1,1)+pd_cdf(2,1); 
pd_unattended_total = pd_cdf(2,1)+pd_cdf(2,2); 
b(1) = bar(1,pd_attended_total/2);
b(2) = bar(2,pd_unattended_total/2); 
b(1).FaceColor = p.style.attColorsMuted(3,:); 
b(1).EdgeColor = 'none'; 
b(2).FaceColor = p.style.attColorsMuted(1,:); 
b(2).EdgeColor = 'none'; 
ylabel('p("saw stimulus")')
ylim([0 1])
yticks([0 1])
xlim([0.3 2.7])
xticks([1 2])
xticklabels({'Attended','Unattended'})

%% matched visibility
clear pd_signal pd_cdf
criterion = 0.8; 
mus = [0 0.4625 1];
sigma = 1; 
sigmaScaling = 2/3;
muScaling = 3/2; 

% -- attended noise distribution
pd_signal(1,1) = makedist('Normal','mu',0,'sigma',1*sigmaScaling);
pd_cdf_noise = 1-sum(cdf('Normal',criterion,0,1*sigmaScaling));

% attended signal distribution
mu_att = 1.1; % 0.7363*muScaling
pd_signal(2,1) = makedist('Normal','mu',mu_att,'sigma',1*sigmaScaling);
pd_cdf_signal = 1-sum(cdf('Normal',criterion,mu_att,sigmaScaling));
pd_attended_total = pd_cdf_noise + pd_cdf_signal;  % 0.9682, 0.5797

% -- unattended noise distribution
pd_signal(1,2) = makedist('Normal','mu',0,'sigma',1);
pd_cdf_noise_un = 1-sum(cdf('Normal',criterion,0,1));

% unattended signal distribution
pd_signal(2,2) = makedist('Normal','mu',1,'sigma',1);
pd_cdf_signal_un = 1-sum(cdf('Normal',criterion,1,1)); % 0.4625
pd_unattended_total = pd_cdf_noise_un + pd_cdf_signal_un;  % 0.7911

d_att = (mu_att)/(sigma*sigmaScaling); 
d_unatt = (1)/(sigma); % 0.4625

figure
set(gcf,'Position',[100 100 500 200])
figureStyle
for iX = [1 2]
    for iA = 1:2 % 1 attended, 2 unattended
        pd = pd_signal(iX,iA);
        pl(iX,iA) = plot(pd);
        pl(iX,iA).LineWidth = 1.5;
        
        xIdx = find(pl(iX,iA).XData>criterion);
        pa(iX,iA) = area(pl(iX,iA).XData(xIdx),pl(iX,iA).YData(xIdx));
        pa(iX,iA).FaceAlpha = 0.4;
        pa(iX,iA).EdgeColor = 'none'; 

        switch iX
            case 1
                pl(iX,iA).LineStyle = '--';
            case 2
                pl(iX,iA).LineStyle = '-';
        end

        switch iA
            case 1
                pl(iX,iA).Color = p.style.attColorsMuted(3,:);
                pa(iX,iA).FaceColor = p.style.attColorsMuted(3,:); 
            case 2
                pl(iX,iA).Color = p.style.attColorsMuted(1,:);
                pa(iX,iA).FaceColor = p.style.attColorsMuted(1,:); 
        end
  
    end
end
ax = gca; 
ax.YAxis.Visible = 'off'; % remove y-axis
xlim([-3 4])
ylim([0 0.7])
xline(0.8,'-k','saw stimulus','LabelOrientation','horizontal')
xlabel('Internal signal strength','FontSize',14)
xticks('')
yticks('')
ylabel([])
legend([pl(1,1) pl(1,2)],{'valid','invalid'},'FontSize',9,'Box','off')
text(-3,0.67,'At matched visibility')

% --- Bar ---
figure
set(gcf,'Position',[100 100 300 200])
subplot 121 
figureStyle
hold on 
b(1) = bar(1,d_att);
b(2) = bar(2,d_unatt); 
b(1).FaceColor = p.style.attColorsMuted(3,:); 
b(1).EdgeColor = 'none'; 
b(2).FaceColor = p.style.attColorsMuted(1,:); 
b(2).EdgeColor = 'none'; 
ylabel('{\it d''}')
ylim([0 2.5])
yticks([0 2.5])
xlim([0.3 2.7])
xticks([1 2])
xticklabels({'Attended','Unattended'})

subplot 122 
figureStyle
hold on 
b(1) = bar(1,pd_attended_total/2);
b(2) = bar(2,pd_unattended_total/2); 
b(1).FaceColor = p.style.attColorsMuted(3,:); 
b(1).EdgeColor = 'none'; 
b(2).FaceColor = p.style.attColorsMuted(1,:); 
b(2).EdgeColor = 'none'; 
ylabel('p("saw stimulus")')
ylim([0 1])
yticks([0 1])
xlim([0.3 2.7])
xticks([1 2])
xticklabels({'Attended','Unattended'})

%% d' by p(seen)
figure
set(gcf,'Position',[100 100 300 200])
subplot 121 
figureStyle
hold on 
b(1) = bar(1,d_att);
b(2) = bar(2,d_unatt); 
b(1).FaceColor = p.style.attColorsMuted(3,:); 
b(1).EdgeColor = 'none'; 
b(2).FaceColor = p.style.attColorsMuted(1,:); 
b(2).EdgeColor = 'none'; 
ylabel('{\it d''}')
ylim([0 2.5])
yticks([0 2.5])
xlim([0.3 2.7])
xticks([1 2])
xticklabels({'Attended','Unattended'})

subplot 122 
figureStyle
hold on 
b(1) = bar(1,pd_attended_total);
b(2) = bar(2,pd_unattended_total); 
b(1).FaceColor = p.style.attColorsMuted(3,:); 
b(1).EdgeColor = 'none'; 
b(2).FaceColor = p.style.attColorsMuted(1,:); 
b(2).EdgeColor = 'none'; 
ylabel('p("saw stimulus")')
ylim([0 1])
yticks([0 1])
xlim([0.3 2.7])
xticks([1 2])
xticklabels({'Attended','Unattended'})









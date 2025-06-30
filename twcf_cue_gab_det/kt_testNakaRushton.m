

%% NakaRushton doing? 
x = data_x(1):0.001:data_x(end); % linear not log space 
n = 1; 
c50 = data_x(end)/4; 
M = 0.1; 
Rmax = 5; % not bounded by 1, infinite? 
% whereas Weibull and Gumbel are 

slopes = 0.1:0.5:10; 

for iS = 1:numel(slopes)
    for i = 1:numel(x)
        n = slopes(iS); 
        y(i,iS) = (Rmax*x(i).^n)/(x(i).^n+c50^n)+M; 
        % (Rmax*x(i)).^n)/(x(i).^slopes(iS)+log(c50)^slopes(iS))+M;
    end
end

%% 
figure
hold on 
for iS = 1:numel(slopes)
    plot(x,y(:,iS))
end
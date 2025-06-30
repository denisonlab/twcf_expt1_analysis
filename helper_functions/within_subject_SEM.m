function SEM = within_subject_SEM(M)
% SEM = within_subject_SEM(M)
% 
% M is a i x j matrix, where rows i denote subject and columns j denote 
% within-subject condition.
% 
% SEM is the within-subject SEM from the data set M, as detailed in 
% 
% Morey, R. D. (2008). Confidence intervals from normalized data: 
%    A correction to Cousineau (2005). Tutorials in Quantitative Methods 
%    for Psychology, 4, 61–64.

grand_mean = mean(mean(M));
subject_means = mean(M,2);

n_subs  = size(M,1);
n_cells = size(M,2);

for i=1:n_cells
    norm_M(:,i) = M(:,i) - subject_means + grand_mean;
end

var_norm_M = var(norm_M) * (n_cells / (n_cells-1));
std_norm_M = sqrt(var_norm_M);
SEM = std_norm_M / sqrt(n_subs);
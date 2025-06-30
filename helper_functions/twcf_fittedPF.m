function eq = twcf_fittedPF(alpha, beta, gamma, lambda, x, p)
% x are stimulus instensities 

switch p.fit.PFtype
    case 'gumbel' % log weibull 
        F = 1 - exp(-10.^(beta*(x-alpha))); 
    case 'weibull'
        F = 1 - exp(-(x/alpha).^beta); 
end

eq = gamma + (1-gamma-lambda) * F; 
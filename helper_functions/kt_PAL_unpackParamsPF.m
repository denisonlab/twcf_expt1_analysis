%
%PAL_unpackParamsPF     Re-format parameters of PF input
%
%syntax: [alpha beta gamma lambda] = PAL_unpackParamsPF(params)
%
%Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)
%Modified: kt to unpack Naka-rushton parameters

function [alpha, beta, gamma, lambda, M] = kt_PAL_unpackParamsPF(params)

gamma = 0;
lambda = 0;

if isstruct(params)
    alpha = params.alpha;
    beta = params.beta;
    if isfield(params,'gamma')
        gamma = params.gamma;
        if isfield(params,'lambda')
            lambda = params.lambda;
            if isfield(params,'M')
                M = params.M; 
            end
        end
    end
else
    alpha = params(1);
    beta = params(2);
    if length(params) > 2
        gamma = params(3);
        if length(params) > 3
            lambda = params(4);
            if length(params) > 4
                M = params(5);
            end
        end
    end
end
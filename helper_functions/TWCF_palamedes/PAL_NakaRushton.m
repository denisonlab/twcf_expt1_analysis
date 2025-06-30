%
%PAL_NakaRushton  Evaluation of Naka-Rushton Psychometric Function
%
%   syntax: y = PAL_NakaRushton(params, x)
%
%   y = PAL_NakaRushton(params, x), where 'params' contains the four
%   parameters of a Psychometric Funtion (i.e., [alpha beta gamma lambda]),
%   returns the Psychometric Function evaluated at values in 'x'. 'x' is
%   array of any size.
%
%   x = PAL_NakaRushton(params, y, 'Inverse') returns the x-value at 
%   which the Psychometric Function evaluates to y.
%
%   dydx = PAL_NakaRushton(params, x, 'Derivative') returns the
%   derivative (slope of tangent line) of the Psychometric Function
%   evaluated at x.
%
%   'params' need not have four entries. A two element vector will be
%   interpreted as [alpha beta], a three element vector as [alpha beta
%   gamma]. Missing elements in 'params' will be assigned a value of 0.
% 
%   Is it possible to have 5 entries? For attentional scaling M_v (valid) and M_i (invalid)?
%
%   This example returns the function value at threshold when gamma 
%   ('guess-rate') and lambda ('lapse-rate') both equal 0:
%       
%   y = PAL_NakaRushton([1 2 0 0], 1)
%
%   y = 0.5000
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.0.2, 1.1.1, 1.2.0, 1.6.3 (see History.m)

function y = PAL_NakaRushton(params,x,varargin)

% [alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params);
% c50, 

% if nargin<1
%     alpha = []; % c50
%     lambda = []; % 1-Rmax
%     beta = []; % n
%     gamma = []; % M 
% else
%     alpha = params.alpha; % threshold (e.g. contrast at 50% performance)
%     lambda = params.lambda; % asymptote (1-lapse rate)
%     beta = params.beta; % slope parameter 
%     gamma = params.gamma; % performance at lowest stimulus strength (i.e. guess rate)
% end

[alpha, beta, gamma, lambda] = kt_PAL_unpackParamsPF(params); 

y = ((1-lambda).*x.^beta)./(x.^beta+alpha.^beta)+gamma; 

% if ~isempty(varargin)
%     if strncmpi(varargin{1}, 'Inverse',3)
%         c = (x - gamma)./(1 - gamma - lambda);
%         c = (1 - c)./c;
%         y = alpha - log(c)./beta;    
%     end
%     if strncmpi(varargin{1}, 'Derivative',3)
%         y = (1 - gamma - lambda).*(1+exp(-1*(beta).*(x-alpha))).^-2.*(exp(-1*(beta).*(x-alpha))).*beta;
%     end
% else
%     y = gamma + (1 - gamma - lambda).*(1./(1+exp(-1*(beta).*(x-alpha))));
% end
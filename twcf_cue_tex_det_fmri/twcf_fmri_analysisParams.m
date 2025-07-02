function p = twcf_analysisParams
% p = twcf_analysisParams for fmri_cue_tex_det

%% Figure styling
p.style.sz             = 100; % scatterplot marker size 
p.style.szSml          = 30; % scatterplot marker size (subjects) 
p.style.alpha          = 0.6; % scatterplot marker transparency 
p.style.xBuffer        = 0.01; % buffer for x limits on both sides 
p.style.xOffset        = 0.15; % offset for group error bars 
p.style.attColors       = [197 77 51; % invalid (red) 
                         81 160 215; % valid (blue)
                         113 112 183]/255; % difference (purple)
%                          0 0 0; % neutral (black)
% p.style.attColors      = [0.8 0 0;  % invalid
%                           0 0 0;    % neutral 
%                           0 1 0];   % valid
p.style.attColorsMuted = [243/255 193/255 193/255;
                          169/255 169/255 169/255;
                          194/255 255/255 193/255];
p.style.greyColor      = [0.9 0.9 0.9]; 
p.style.labelColors    = [0.9 0.9 0.9]; % light grey 
                        % [133/255 174/255 227/255];  % light blue 
                        % 55/255  64/255  139/255]; % dark blue 
p.style.figHeight     = 700; 
p.style.errCapSize    = 0; 

%% Palamedes PF and search grid 
p.fit.PFtype          = 'weibull'; % 'gumbel' (for log stim strength), 'weibull' (for linear stim strength) 
p.fit.paramsFreeA     = [1 1 0 1]; % alpha threshold, beta slope, gamma guess-rate, lambda lapse-rate, 1-->free, 0-->fixed, guesss-rate is fixed 
p.fit.paramsFreeC     = [1 1 1 1];

switch p.fit.PFtype 
    case 'gumbel' % log space 
        p.fit.PF               = @PAL_Gumbel; % type of psychometric function type to fit: @PAL_Logistic, @PAL_Gumbel, @PAL_Weibull
        p.fit.searchGrid.alpha = log10(.05:.05:3); % threshold 
    case 'weibull' % linear space 
        p.fit.PF               = @PAL_Weibull; 
        p.fit.searchGrid.alpha = .05:.05:3; % threshold 
end
p.fit.betaCorrection    = 0; % 0.001; % prevents step function 
p.fit.searchGrid.beta   = 10.^(-1:.1:1-p.fit.betaCorrection); % slope 
p.fit.searchGrid.lambda = 0:.01:.1; % lapse-rate before -:.001:.1

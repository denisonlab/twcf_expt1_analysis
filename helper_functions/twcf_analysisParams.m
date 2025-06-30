function p = twcf_analysisParams
% p = twcf_analysisParams

%% Figure styling
p.style.sz              = 100; % scatterplot marker size 
p.style.szSml           = 30; % scatterplot marker size (subjects) 
p.style.alpha           = 0.6; % scatterplot marker transparency 
p.style.xBuffer         = 0.05; % buffer for x limits on both sides % for line length 
% p.style.xBufferLLDis    = 0.1; % buffer for x limits for texture dis expt 
p.style.xBufferContrast = 0.015; % 0.005; % buffer for x limits for contrast 
p.style.xOffset         = 0.15; % offset for group error bars 
% p.style.attColors      = [0.8 0 0;  % invalid
%                           0 0 0;    % neutral 
%                           0 1 0];   % valid
% colorblind safe palette 
p.style.attColors       = [197 77 51; % invalid (red) 
                         0 0 0; % neutral (black) 
                         81 160 215; % valid (blue)
                         113 112 183]/255; % difference (purple)
% p.style.attColors      = [214/255 66/255 38/255; % invalid (red) 
%                         0/255 0/255 0/255; % neutral (black) 
%                         38/255 163/255 221/255]; % valid (blue) 
p.style.attColorsMuted  = [214 131 112;
                           97 97 97;
                           133 188 227]/255; 
p.style.attColorsMutedLight = [214 162 147;
    128 128 128;
    166 204 246;
    178 178 219]/255; 
p.style.attColorsExtraLight = [230 192 170;
                                195 195 195;
                                199 230 246]/255; 
p.style.attColorsKnit(:,:,1) = [226 181 175; % light
                                214 214 214;
                                175 202 206]/255; 
p.style.attColorsKnit(:,:,2) = [173 98 92; % dark 
                            76 76 76;
                            85 132 160]/255;
p.style.greens(:,1) = [201 219 198]/255; % stimulus vis (light green) 
p.style.greens(:,2) = [69 96 27]/255; % feature vis (dark green) 

% p.style.attColorsMuted = [243/255 193/255 193/255;
%                           169/255 169/255 169/255;
%                           141/255 176/255 255/255];
p.style.greyColor      = [0.9 0.9 0.9]; 
p.style.labelColors    = [0.9 0.9 0.9]; % light grey 
                        % [133/255 174/255 227/255];  % light blue 
                        % 55/255  64/255  139/255]; % dark blue 

p.style.figHeight     = 700; 
p.style.errCapSize    = 0; 
p.style.fitLineWidth  = 2; % fit line width 

% Figure styling 
p.style.textTickSize  = 14; % 14
p.style.textAxisSize  = 18; % 22, 18 

p.style.xtickangle    = 30; 

% Figure axes sizing 
p.size.AUC = [0.4 0.2074 0.45 0.7176]; 
p.size.AUCbySite = [0.3 0.2074 0.5 0.7176]; 
p.size.rect = [ 0 0 380 280]; % 360 260
p.size.rect2 = [ 0 0 200 260];

% Significance stars styling 
p.sigStarMainEffect = 0.95; % percent of ylim to plot the significance label for a main effect

%% Palamedes PF and search grid 
p.fit.PFtype          = 'weibull'; % 'nakarushton'; % 'gumbel' (for log stim strength), 'weibull' (for linear stim strength), 'nakarushton'
p.fit.paramsFreeA     = [1 1 0 1]; % alpha threshold, beta slope, gamma guess-rate, lambda lapse-rate, 1-->free, 0-->fixed, guesss-rate is fixed 
p.fit.paramsFreeC     = [1 1 1 1];

switch p.fit.PFtype 
    case 'gumbel' % log space 
        p.fit.PF                = @PAL_Gumbel; % type of psychometric function type to fit: @PAL_Logistic, @PAL_Gumbel, @PAL_Weibull
        p.fit.searchGrid.alpha  = log10(.05:.05:3); % threshold
        p.fit.betaCorrection    = 0; % 0.001; % prevents step function
        p.fit.searchGrid.beta   = 10.^(-1:.1:1-p.fit.betaCorrection); % slope
        p.fit.searchGrid.lambda = 0:.01:.1; % lapse-rate before -:.001:.1

    case 'weibull' % linear space 
        p.fit.PF                = @PAL_Weibull;
        p.fit.searchGrid.alpha  = .05:.05:3; % threshold
        p.fit.betaCorrection    = 0; % 0.001; % prevents step function
        p.fit.searchGrid.beta   = 10.^(-1:.1:1-p.fit.betaCorrection); % slope
        p.fit.searchGrid.lambda = 0:.01:.1; % lapse-rate before -:.001:.1

    case 'nakarushton' % will name these w the same convention as gumbel, weibulls; but note that these are not the standard naming conventions for NakaRushtons
        p.fit.PF                = @PAL_NakaRushton; 
        p.fit.searchGrid.alpha  = .05:.05:3; % threshold (e.g. contrast at 50% performance) alpha equivalent 
        p.fit.betaCorrection    = 0; % 0.001; % prevents step function
        p.fit.searchGrid.beta   = 10.^(-1:.1:1-p.fit.betaCorrection); % slope parameter, is this reasonable range? plot some nakarushtons to check? 
        % p.fit.gamma           = []; % performance at lowest stimulus strength (i.e. guess rate); gamma, 50% for A, free for C?  
        p.fit.searchGrid.lambda = 0.5:.01:1; % asymptote (1-lapse rate) related to lambda 
end

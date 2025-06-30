function baseDir = twcf_baseDir(user)

% baseDir = TWCF_BASEDIR(user)
%   Input: 
%       - user (string)
%           - eg 'karen' 'emily'
%   Output:
%       - baseDir (string)
%         returns base directory for user 

%% 
switch user
    case 'karen'
        baseDir = '/Users/kantian/Dropbox/github/TWCF_FOHO';
    case 'emily'
        baseDir = '/Users/emilyrussell/Documents/GitHub/TWCF_FOHO'; % Emily add here 
    otherwise 
        error('user not defined')
end
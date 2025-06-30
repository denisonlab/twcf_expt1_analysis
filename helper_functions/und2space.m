function a = und2space(a)

% function a = und2space(a)
%
% in a string a, replaces underscores with spaces

a(strfind(a,'_')) = ' ';
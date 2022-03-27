% Modify the string below so that it points to the parent directory of 
% your CmdStan installation
% TODO
%  o some basic checking
%  o some way to manage fileseparators?
function d = stan_home()

if ispc
 d = 'C:\cmdstan';
elseif ismac
d = '/Users/yanndufour/Documents/Github/cmdstan';
end
%d = 'C:\Users\brian\Downloads\cmdstan';

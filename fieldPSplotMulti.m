function fieldPSplotMulti(Fname,tOffset,tLimit,pathCell,scale)
%
% fieldPSplotMulti(Fname,tOffset,tLimit,pathCell,scale)
%
% Plots field power spectra from several directories, on the same graph.
% Useful for comparing between resolutions or values of s
%
%  Fname = field name, eg. E2, B2, Phi2, or dotPhi
%  tOffset = offset to improve scaling (not working yet, use any number)
%  tLimit = times between which to plot
%  inPath = cell array of paths to files, including final '/'
% e.g. inpath = {'1024/resolution/','2048/resolution/','4096/resolution/'};
%  scale = 'log' (default) or 'linear'
%
% Version 1.1 2013.6.11 MBH

if ~exist('scale','var')
    scale = 'log';
end

clf

for n=1:numel(pathCell)

inpath = pathCell{n};
h = subplot(1,1,1);
fieldPSplot(Fname,'',1,tOffset,tLimit,inpath,h);
set(h,'YScale',scale)

end

function [k,t]=fieldPSinfoRead(path,id)
%
% [k,t]=fieldPSinfoRead(path,id)
%
% Version 1.0 2013.6.10 MBH

fileName=[path 'resInfo' id '.txt'];
disp(['Loading field PS infoFile: ' fileName])
fid=fopen(fileName,'r');

% temp=fgetl(fid);   %Read comment
% nk=fscanf(fid,'%i',[1 1]);
% temp=fgetl(fid);   %Finish line

temp=fgetl(fid);   % Discard comment line
[k nk] = fscanf(fid,'%f,'); % Read and count the k-values

temp=fgetl(fid);   

temp=fgetl(fid);   % Discard comment line   
temp=fgetl(fid);   % Discard comment line   

[stats nt] = fscanf(fid,'%f,%f,%f,%f,%f,%f,%f,%f,%f',inf);

stats = reshape(stats,9,nt/9)';

t = stats(:,1);

fclose(fid);
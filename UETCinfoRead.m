function [k,r]=UETCinfoRead(path,id)
fileName=[path 'info' id '.dat'];
disp(['Loading infoFile: ' fileName])
fid=fopen(fileName,'r');

temp=fgetl(fid);   %Read comment
nk=fscanf(fid,'%i',[1 1]);
temp=fgetl(fid);   %Finsih line

temp=fgetl(fid);   
k=fscanf(fid,'%f',[1 nk]);
temp=fgetl(fid);   

temp=fgetl(fid);   
nt=fscanf(fid,'%i',[1 1]);
temp=fgetl(fid);   

temp=fgetl(fid);   
r=fscanf(fid,'%f',[1 nt])';
temp=fgetl(fid); 

fclose(fid);

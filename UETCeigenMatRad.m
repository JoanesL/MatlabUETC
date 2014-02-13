%Function to process UETCeigen output to ensure
%radiation era eigenvectors have same sign as
%corresponding matter era ones (for rad->mat 
%transition interpolation)
%
%Run UETCeigen for both Mat and Rad eras, then
%run UETCeigenMatRad.
%
%Usage: UETCeigenMatRad(outPath, number, inPath)
%
%outpath = same outPath used int UETCeigen
%          (ie. path to UETCeigen data relative
%           to inPath)
% number = ID number(s) of UETCeigen set
%
%Option input:
%
% inPath = overide gpath for UETCeigen input path
%
% JO: Ni = size of the interpolated matrix. Remember, scalar 2Ni x 2Ni
%
% JO: Eval. How many eigenvalues
%
%Can do mutiple ID numbers at once.

function UETCeigenMatRad(outPath, number, inPath, Ni, Eval)

if nargin==0; 
  help UETCeigenMatRad
  return
end

if ~exist('Ni','var'); Ni=''; end %JO

global gpath

if ~exist('inPath','var'); inPath=''; end 

if prod(size(inPath))>0; 
  path=inPath; 
else
  if prod(size(gpath))>0
    path=gpath;
  else
    disp(['Please set gpath global variable to default path'...
	  ' or specify path in fucntion call'])
    return
  end
end

path = [path outPath];

v=1:Eval;
t=1:Eval;
[V,T]=meshgrid(v,t);

for b=number

%=======
%SCALARS
%=======  

%Load scalar eigenvectors
%disp(['Loading: ' path 'UETCeigenScalarMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'])
%M=load([path 'UETCeigenScalarMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
%disp(['Loading: ' path 'UETCeigenScalarRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'])
%R=load([path 'UETCeigenScalarRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);

disp(['Loading: ' path 'UETCeigenScalarMat_' num2str(b,'%2.2d') '.dat'])
M=load([path 'UETCeigenScalarMat_' num2str(b,'%2.2d') '.dat']);
disp(['Loading: ' path 'UETCeigenScalarRad_' num2str(b,'%2.2d') '.dat'])
R=load([path 'UETCeigenScalarRad_' num2str(b,'%2.2d') '.dat']);

%Check for same sign normalization
dkt=diff(M(1,2:(end+1)/2));
for i=2:Eval
 signPhi=sum(dkt.*M(i,2:(end+1)/2-1).*R(i,2:(end+1)/2-1));
 signPsi=sum(dkt.*M(i,(end+3)/2:end-1).*R(i,(end+3)/2:end-1));
 sign=signPhi+signPsi;
 if signPhi<0 || signPsi<0
  disp(['S' num2str(i-1) ': ' num2str(signPhi) '   ' num2str(signPsi) '   =   '  num2str(sign)]) 
 end
 if sign<0
 R(i,2:end)=-R(i,2:end);
 end
end

%Save scalar eigenvectors for radiation
disp('Resaving scalar radiation and matter era data')
save([path 'UETCeigenScalarRad_' num2str(b,'%2.2d') '.dat'],'R','-ascii')
save([path 'UETCeigenScalarMat_' num2str(b,'%2.2d') '.dat'],'M','-ascii')

%=======
%VECTORS
%=======

%Load vector eigenvectors
%disp(['Loading: ' path 'UETCeigenVectorMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'])
%M=load([path 'UETCeigenVectorMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
%disp(['Loading: ' path 'UETCeigenVectorRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'])
%R=load([path 'UETCeigenVectorRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);

disp(['Loading: ' path 'UETCeigenVectorMat_' num2str(b,'%2.2d') '.dat'])
M=load([path 'UETCeigenVectorMat_' num2str(b,'%2.2d') '.dat']);
disp(['Loading: ' path 'UETCeigenVectorRad_' num2str(b,'%2.2d') '.dat'])
R=load([path 'UETCeigenVectorRad_' num2str(b,'%2.2d') '.dat']);

%Check for same sign normalization
dkt=diff(M(1,2:end));
for i=2:Eval
 sign=sum(dkt.*M(i,2:end-1).*R(i,2:end-1));
 if sign<0
  disp(['V' num2str(i-1) ': ' num2str(sign)]) 
  R(i,2:end)=-R(i,2:end);
 end
end

%Save vector eigenvectors for radiation
disp('Resaving vector radiation era data')
save([path 'UETCeigenVectorRad_' num2str(b,'%2.2d') '.dat'],'R','-ascii')
save([path 'UETCeigenVectorMat_' num2str(b,'%2.2d') '.dat'],'M','-ascii')


%=======
%TENSORS
%=======

%Load tensor eigenvectors
%disp(['Loading: ' path 'UETCeigenTensorMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'])
%M=load([path 'UETCeigenTensorMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
%disp(['Loading: ' path 'UETCeigenTensorRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'])
%R=load([path 'UETCeigenTensorRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
disp(['Loading: ' path 'UETCeigenTensorMat_' num2str(b,'%2.2d') '.dat'])
M=load([path 'UETCeigenTensorMat_' num2str(b,'%2.2d') '.dat']);
disp(['Loading: ' path 'UETCeigenTensorRad_' num2str(b,'%2.2d') '.dat'])
R=load([path 'UETCeigenTensorRad_' num2str(b,'%2.2d') '.dat']);

%Check for same sign normalization
dkt=diff(M(1,2:end));
for i=2:Eval
 sign=sum(dkt.*M(i,2:end-1).*R(i,2:end-1));
 if sign<0
  disp(['T' num2str(i-1) ': ' num2str(sign)]) 
  R(i,2:end)=-R(i,2:end);
 end
end

%Save Tensor eigenvectors for radiation
disp('Resaving tensor radiation era data')
save([path 'UETCeigenTensorRad_' num2str(b,'%2.2d') '.dat'],'R','-ascii')
save([path 'UETCeigenTensorMat_' num2str(b,'%2.2d') '.dat'],'M','-ascii')


end

%Joanes (2012/11/28)
%Function to Plot Matter vs Radiation eigenvectors
%
%Usage=UETCeigenCheck(inPath, type, ID,Evalmin, eigenvalueLimit,Ni)
%
%inPath = path to folder which contains eigen data
%ID = number corresponding to data file
%Evalmin = starting point for plotting
%eigenvalueLimit = eigenvectors will be plotted up to this value
%Optional:
%Ni = Size of the interpolated matrix
%BefAf =    -0 check eigenvectors Before interpolation 
%           -1 check After interpolation


function UETCeigenCheck(path,ID,Evalmin,lim,Ni,BefAf)

if nargin==0; 
  help UETCeigenCheck
  return
end

if ~exist('Ni','var'); Ni=''; end


%JO added type, for 'Scalar', 'Vector' and 'Tensor'
%JO added ID in order to get correctly UETCeigen data
b=ID;



if BefAf == 0 
       
    
        %JO ...JO  before interpolation
    %Sm=load([path 'UETCeigenScalarBeforeMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    %Sr=load([path 'UETCeigenScalarBeforeRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    %Vm=load([path 'UETCeigenVectorBeforeMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    %Vr=load([path 'UETCeigenVectorBeforeRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    %Tm=load([path 'UETCeigenTensorBeforeMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    %Tr=load([path 'UETCeigenTensorBeforeRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    Sm=load([path 'UETCeigenScalarBeforeMat_' num2str(b,'%2.2d') '.dat']);
    Sr=load([path 'UETCeigenScalarBeforeRad_' num2str(b,'%2.2d') '.dat']);
    Vm=load([path 'UETCeigenVectorBeforeMat_' num2str(b,'%2.2d') '.dat']);
    Vr=load([path 'UETCeigenVectorBeforeRad_' num2str(b,'%2.2d') '.dat']);
    Tm=load([path 'UETCeigenTensorBeforeMat_' num2str(b,'%2.2d') '.dat']);
    Tr=load([path 'UETCeigenTensorBeforeRad_' num2str(b,'%2.2d') '.dat']);

elseif BefAf == 1    
    
    
    
    %JO ...JO added after interpolation
    Sm=load([path 'UETCeigenScalarMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    Sr=load([path 'UETCeigenScalarRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    Vm=load([path 'UETCeigenVectorMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    Vr=load([path 'UETCeigenVectorRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    Tm=load([path 'UETCeigenTensorMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    Tr=load([path 'UETCeigenTensorRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat']);
    
    Ni=8*Ni;
end    
    
    %%%%%%%%%%JO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
  disp(['Size scalar' num2str(size(Sm))])
    figure
    scalar=multiPlot([2 4])
    figure
    vector=multiPlot([3 2])
    figure
    tensor=multiPlot([3 3])
    
    
    for i=Evalmin:lim
    ktmphi=Sm(1,2:Ni+1);
    ktrphi=Sr(1,2:Ni+1);
    
    axes(scalar(1))
    semilogx(ktmphi,Sm(i+1,2:Ni+1),'g',ktrphi,Sr(i+1,2:Ni+1),'r');
    title(['Scalar Phi:MAT->GREEN, RAD->RED' num2str(i)],'FontWeight','bold')
    
    axes(scalar(3))
    semilogx(ktmphi,Sm(i,2:Ni+1),'g',ktrphi,Sr(i+1,2:Ni+1),'r');
    title(['MAT' num2str(i-1) 'vs RAD' num2str(i) 'Phi'],'FontWeight','bold')
    
    axes(scalar(5))
    semilogx(ktmphi,Sm(i+1,2:Ni+1),'g',ktrphi,Sr(i,2:Ni+1),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i-1) 'Phi'],'FontWeight','bold')
    
    axes(scalar(7))
    semilogx(ktmphi,Sm(i+1,2:Ni+1),'g',ktrphi,Sr(i+3,2:Ni+1),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i+2) 'Phi'] ,'FontWeight','bold')
    
    ktmpsi=Sm(1,Ni+2:end);
    ktrpsi=Sr(1,Ni+2:end);
    axes(scalar(2))
    semilogx(ktmpsi,Sm(i+1,Ni+2:end),'g',ktrpsi,Sr(i+1,Ni+2:end),'r');
    title(['Scalar Psi:MAT->GREEN, RAD->RED' num2str(i)],'FontWeight','bold')
     
    axes(scalar(4))
    semilogx(ktmpsi,Sm(i,Ni+2:end),'g',ktrpsi,Sr(i+1,Ni+2:end),'r');
    title(['MAT' num2str(i-1) 'vs RAD' num2str(i) 'Psi'],'FontWeight','bold')

    axes(scalar(6))
    semilogx(ktmpsi,Sm(i+1,Ni+2:end),'g',ktrpsi,Sr(i,Ni+2:end),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i-1) 'Psi'] ,'FontWeight','bold')
    
    axes(scalar(8))
    semilogx(ktmpsi,Sm(i+1,Ni+2:end),'g',ktrpsi,Sr(i+3,Ni+2:end),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i+2) 'Psi'] ,'FontWeight','bold')
    
    
    waitforbuttonpress
    
    end
    for i=Evalmin:lim
    ktvtm=Vm(1,2:end);
    ktvtr=Vr(1,2:end);
    
    axes(vector(1))
    semilogx(ktvtm,Vm(i+1,2:end),'g',ktvtr,Vr(i+1,2:end),'r');
    title(['Vector:MAT->GREEN, RAD->RED' num2str(i)],'FontWeight','bold')
    
    axes(vector(2))
    semilogx(ktvtm,Vm(2,2:end),'g');%,ktvtr,Vr(i+1,2:end),'r');
    title(['Vector:MAT->GREEN, RAD->RED' num2str(i)],'FontWeight','bold')
    
    axes(vector(3))
    semilogx(ktvtm,Vm(i,2:end),'g',ktvtr,Vr(i+1,2:end),'r');
    title(['Rad' num2str(i-1) 'vs Mat' num2str(i)],'FontWeight','bold')
    
    axes(vector(4))
    semilogx(ktvtm,Vm(i+1,2:end),'g',ktvtr,Vr(i,2:end),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i-1)],'FontWeight','bold')
    
    axes(vector(5))
    semilogx(ktvtm,Vm(i+1,2:end),'g',ktvtr,Vr(i+2,2:end),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i+1)],'FontWeight','bold')
    
    axes(vector(6))
    semilogx(ktvtm,Vm(i+1,2:end),'g',ktvtr,Vr(i+3,2:end),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i+2)],'FontWeight','bold')
     waitforbuttonpress
     end
    for i=Evalmin:lim
      axes(tensor(1))
    semilogx(ktvtm,Tm(i+1,2:end),'g',ktvtr,Tr(i+1,2:end),'r');
       title(['Tensor:MAT->GREEN, RAD->RED' num2str(i)],'FontWeight','bold')
    
    axes(tensor(3))
    semilogx(ktvtm,Tm(i,2:end),'g',ktvtr,Tr(i+1,2:end),'r');
    title(['MAT' num2str(i-1) 'vs RAD' num2str(i)],'FontWeight','bold')
    
    axes(tensor(2))
    semilogx(ktvtm,Tm(i+1,2:end),'g',ktvtr,Tr(i+2,2:end),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i+1)],'FontWeight','bold')
    
    axes(tensor(2))
    semilogx(ktvtm,Tm(i+1,2:end),'g',ktvtr,Tr(i+4,2:end),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i+3)],'FontWeight','bold')
    
    axes(tensor(4))
    semilogx(ktvtm,Tm(i+1,2:end),'g',ktvtr,Tr(i,2:end),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i-1)],'FontWeight','bold')
    
     axes(tensor(5))
    semilogx(ktvtm,Tm(i+1,2:end),'g',ktvtr,Tr(i+3,2:end),'r');
    title(['MAT' num2str(i) 'vs RAD' num2str(i+2)],'FontWeight','bold')
    
     axes(tensor(6))
    semilogx(ktvtm,Tm(i+3,2:end),'g',ktvtr,Tr(i+1,2:end),'r');
    title(['MAT' num2str(i+2) 'vs RAD' num2str(i)],'FontWeight','bold')
    
    
    waitforbuttonpress
    
    end
   
end


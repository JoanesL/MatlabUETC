%JO 2012
%Function to write eigenvectors in a suitable way to be inserted in
%CMBeasy. It reads original UETC eigenvectors and eigenvalues. It makes the
%correct links between them, if some of them are mixed. After that checks
%the relative sign. Finally performs the interpolation. 
   %Usage: UETCoutputCMBeasy(id,inPath,outPath,Ni,Evalmin,Eval,interp,plot)
%
%id = 'UETCeigenScalarBeforeMat512_' id '.dat'
%           e.g: id=1
%inpath = folder which contains data
%outpath = output path relative to inpath
%Ni = size of the matrix
%Evalmin = min to plot
%Eval = number of eigenvalues to insert in CMBeasy 
%interp = 
%           -1 interpolate
%           -else do not interpolate
%plot = 
%           -1 Cross dot product of eigenvectors before interpolation
%           -2 Cross dot product of eigenvectors before and after interpolation
%           -0 nothing
%sign =     -1 check relative sign
%           -0 do nothing
function UETCoutputCMBeasy(id,inPath,outPath,Ni,Evalmin,Evaldim,interp,plot,negativeOrder)

if nargin==0; 
  help UETCoutputCMBeasy
  return
end

if prod(size(inPath))>0; 
  path=inPath; 
else,inPath,outPath
  if prod(size(gpath))>0
    path=gpath;
  end
end


b=id;

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

    
    for i=1:Evaldim
    sm.Evector(i,:) = Sm(i+1,2:end);
    sr.Evector(i,:) = Sr(i+1,2:end);
    sm.Evectorphi(i,:) = Sm(i+1,2:Ni+1);
    sm.Evectorpsi(i,:) = Sm(i+1,Ni+2:end);
    sr.Evectorphi(i,:) = Sr(i+1,2:Ni+1);
    sr.Evectorpsi(i,:) = Sr(i+1,Ni+2:end);
    sm.Evalue(i) = Sm(i+1,1);
    sr.Evalue(i) = Sr(i+1,1);

    
    vm.Evector(i,:) = Vm(i+1,2:end);
    vr.Evector(i,:) = Vr(i+1,2:end);
    vm.Evalue(i) = Vm(i+1,1);
    vr.Evalue(i) = Vr(i+1,1);
    
    tm.Evector(i,:) = Tm(i+1,2:end);
    tr.Evector(i,:) = Tr(i+1,2:end);
    tm.Evalue(i) = Tm(i+1,1);
    tr.Evalue(i) = Tr(i+1,1);
    end

   
    ktmphi=Sm(1,2:Ni+1);
    ktrphi=Sr(1,2:Ni+1);
    
    ktmpsi=Sm(1,Ni+2:end);
    ktrpsi=Sr(1,Ni+2:end);
    
    ktvtm=Vm(1,2:end);
    ktvtr=Vr(1,2:end);
    
    if (negativeOrder==1)
        for i=1:Evaldim
            if (sm.Evalue(i)<0)
                disp(['scalar Negatibo! ' num2str(i) ])
                
                temp_s=sm.Evalue(i);
                
                sm.Evalue(i)=[];
                
                sm.Evalue(end+1)=temp_s;                
            end
        end
        EEEOOO=sm.Evalue
        for i=1:Evaldim
            if (vm.Evalue(i)<0)
                disp(['vector Negatibo! ' num2str(i) ])
            
                temp_v=vm.Evalue(i);
                
                vm.Evalue(i)=[];
                
                vm.Evalue(end+1)=temp_v;
             end
        end
        EEEOOO=vm.Evalue
        for i=1:Evaldim
            if (tm.Evalue(i)<0)
                disp(['tensor Negatibo! ' num2str(i) ])
            
                temp_t=tm.Evalue(i);
                
                tm.Evalue(i)=[];
                
                tm.Evalue(end+1)=temp_t;
            end
        end
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%INTERPOLATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if interp==1    
%ktI=logspace(log10(ktvtm(1)),log10(ktvtm(end)),4*Ni);
%ktIs=logspace(log10(ktvtm(1)),log10(ktvtm(end)),8*Ni);
ktI=logspace(-5,log10(ktvtm(end)),4*Ni);
ktIs=logspace(-5,log10(ktvtm(end)),8*Ni);

for i=1:Evaldim 
   sm.EvectorphiI(i,:)=interp1([0 ktmphi],[0 sm.Evectorphi(i,:)], ktIs, 'pchip');
   sr.EvectorphiI(i,:)=interp1([0 ktrphi],[0 sr.Evectorphi(i,:)], ktIs, 'pchip');
   sm.EvectorpsiI(i,:)=interp1([0 ktmpsi],[0 sm.Evectorpsi(i,:)], ktIs, 'pchip');
   sr.EvectorpsiI(i,:)=interp1([0 ktmpsi],[0 sr.Evectorpsi(i,:)], ktIs, 'pchip');
   
   vm.EvectorI(i,:)=interp1([0 ktvtm],[0 vm.Evector(i,:)], ktI, 'pchip');
   vr.EvectorI(i,:)=interp1([0 ktvtr],[0 vr.Evector(i,:)], ktI, 'pchip');
   
   tm.EvectorI(i,:)=interp1([0 ktvtm],[0 tm.Evector(i,:)], ktI, 'pchip');
   tr.EvectorI(i,:)=interp1([0 ktvtr],[0 tr.Evector(i,:)], ktI, 'pchip');
end

SM=[2*max(size(ktIs)) ktIs ktIs; (abs(sm.Evalue))' sm.EvectorphiI sm.EvectorpsiI];
SR=[2*max(size(ktIs)) ktIs ktIs; (abs(sr.Evalue))' sr.EvectorphiI sr.EvectorpsiI];
%save([path 'UETCeigenScalarMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'],'SM','-ascii')
%save([path 'UETCeigenScalarRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'],'SR','-ascii')

save([path 'UETCeigenScalarMat_' num2str(b,'%2.2d') '.dat'],'SM','-ascii')
save([path 'UETCeigenScalarRad_' num2str(b,'%2.2d') '.dat'],'SR','-ascii')

VM=[max(size(ktI)) ktI; (abs(vm.Evalue))' vm.EvectorI];
VR=[max(size(ktI)) ktI; (abs(vr.Evalue))' vr.EvectorI];
%save([path 'UETCeigenVectorMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'],'VM','-ascii')
%save([path 'UETCeigenVectorRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'],'VR','-ascii')

save([path 'UETCeigenVectorMat_' num2str(b,'%2.2d') '.dat'],'VM','-ascii')
save([path 'UETCeigenVectorRad_' num2str(b,'%2.2d') '.dat'],'VR','-ascii')

TM=[max(size(ktI)) ktI; (abs(tm.Evalue))' tm.EvectorI];
TR=[max(size(ktI)) ktI; (abs(tr.Evalue))' tr.EvectorI];
%save([path 'UETCeigenTensorMat' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'],'TM','-ascii')
%save([path 'UETCeigenTensorRad' num2str(Ni) '_' num2str(b,'%2.2d') '.dat'],'TR','-ascii')

save([path 'UETCeigenTensorMat_' num2str(b,'%2.2d') '.dat'],'TM','-ascii')
save([path 'UETCeigenTensorRad_' num2str(b,'%2.2d') '.dat'],'TR','-ascii')

end

    dktmphi=diff(ktmphi(1:end));
    dktmpsi=diff(ktmpsi(1:end));
    dktvtm=diff(ktvtm(1:end));        
    
    v=Evalmin:Evaldim;
    t=Evalmin:Evaldim;
    [V,T]=meshgrid(v,t);

    for i=1:(Evaldim+1-Evalmin)
            k=Evalmin+i-1;
            for j=1:(Evaldim+1-Evalmin)
                    s=Evalmin+j-1;
                    crosssphi(i,j)=sum(dktmphi.*sm.Evectorphi(k,1:end-1).*(sr.Evectorphi(s,1:end-1)));
                    crossspsi(i,j)=sum(dktmpsi.*sm.Evectorpsi(k,1:end-1).*(sr.Evectorpsi(s,1:end-1)));
                    %corsss(i,j)=sum(dktmpsi.*sm.Evector(k,1:end-1).*(sr.Evector(s,1:end-1)));
                    crossv(i,j)=sum(dktvtm.*vm.Evector(k,1:end-1).*(vr.Evector(s,1:end-1)));
                    crosst(i,j)=sum(dktvtm.*tm.Evector(k,1:end-1).*(tr.Evector(s,1:end-1)));
                %crosssphi(i,j)=abs(sm.Evectorphi(k,:)*(sr.Evectorphi(s,:)'));
                %crossstot(i,j)=abs(sm.Evector(k,:)*(sr.Evector(s,:)'));
                %crossv(i,j)=abs(vm.Evector(k,:)*(vr.Evector(s,:)'));
                %crosst(i,j)=abs(tm.Evector(k,:)*(tr.Evector(s,:)'));
            end
    end
       
	    %disp('Lets check if the order was correctly set')
	 
  
if plot==1    
   figure()   
   %plotc(V,T,crosss,Evalmin,Evaldim,1000,0) 
   %title('Cross Scalar','FontWeight','bold')
   plotc(V,T,crosssphi,Evalmin,Evaldim,1000,0) 
   title('Cross Scalar Phi','FontWeight','bold')  
   figure()
   plotc(V,T,crossspsi,Evalmin,Evaldim,1000,0) 
   title('Cross Scalar Psi','FontWeight','bold')
   figure()
   plotc(V,T,crossv,Evalmin,Evaldim,100,0) 
   title('Cross Vector','FontWeight','bold')
   figure()
   plotc(V,T,crosst,Evalmin,Evaldim,100,0) 
   title('Cross Tensor','FontWeight','bold')
   
elseif plot==2
        
    v=Evalmin:Evaldim;
    t=Evalmin:Evaldim;
    [V,T]=meshgrid(v,t);
    
    for i=Evalmin:Evaldim
        sm.EvectorI(i,:) = [sm.EvectorphiI(i,:) sm.EvectorpsiI(i,:)];
        sr.EvectorI(i,:) = [sr.EvectorphiI(i,:) sr.EvectorpsiI(i,:)];
    end
    
    for i=1:(Evaldim+1-Evalmin)
            k=Evalmin+i-1;
            for j=1:(Evaldim+1-Evalmin)
                    s=Evalmin+j-1;
                crossstot(i,j)=abs(sm.Evector(k,:)*(sr.Evector(s,:)'));
                crossstotI(i,j)=abs(sm.EvectorI(k,:)*(sr.EvectorI(s,:)'));
                
                crossv(i,j)=abs(vm.Evector(k,:)*(vr.Evector(s,:)'));
                crossvI(i,j)=abs(vm.EvectorI(k,:)*(vr.EvectorI(s,:)'));
                crosst(i,j)=abs(tm.Evector(k,:)*(tr.Evector(s,:)'));
                crosstI(i,j)=abs(tm.EvectorI(k,:)*(tr.EvectorI(s,:)'));
                %samet(i,j)=abs(tm.Evector(k,:)*(tm.Evector(s,:)'));
                %sametI(i,j)=abs(tm.EvectorI(k,:)*(tm.EvectorI(s,:)'));
            end
       end
    
   figure
   plotc(V,T,crossstot,Evalmin,Evaldim,30,0) 
   title('Cross Scalar, Before interpolation','FontWeight','bold')  
   figure
   plotc(V,T,crossstotI,Evalmin,Evaldim,30,1) 
   title('Cross Scalar, After interpolation','FontWeight','bold') 
   figure
   plotc(V,T,crossv,Evalmin,Evaldim,30,0) 
   title('Cross Vector, Before interpolation','FontWeight','bold')
   figure
   plotc(V,T,crossvI,Evalmin,Evaldim,30,1) 
   title('Cross Vector, After interpolation','FontWeight','bold')
   figure
   plotc(V,T,crosst,Evalmin,Evaldim,30,0) 
   title('Cross Tensor, Before interpolation','FontWeight','bold')
   figure
   plotc(V,T,crosstI,Evalmin,Evaldim,30,1) 
   title('Cross Tensor, After interpolation','FontWeight','bold')
   %figure
   %plotc(V,T,samet,Evalmin,Evaldim,30,0) 
   %title('Same Tensor, Before interpolation','FontWeight','bold')
   %figure
   %plotc(V,T,sametI,Evalmin,Evaldim,30,1) 
   %title('Same Tensor, After interpolation','FontWeight','bold')
end
    

end

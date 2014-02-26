function ETCRecons(path,Ni,id)
    b=id;

%     Sm=load([path 'UETCeigenScalarBeforeMat_' num2str(b,'%2.2d') '.dat']);
%     Sr=load([path 'UETCeigenScalarBeforeRad_' num2str(b,'%2.2d') '.dat']);
%     Vm=load([path 'UETCeigenVectorBeforeMat_' num2str(b,'%2.2d') '.dat']);
%     Vr=load([path 'UETCeigenVectorBeforeRad_' num2str(b,'%2.2d') '.dat']);
     Vm=load([path 'UETCeigenTensorBeforeMat_' num2str(b,'%2.2d') '.dat']);
     Vr=load([path 'UETCeigenTensorBeforeRad_' num2str(b,'%2.2d') '.dat']);
    
    for i=1:Ni
%     sm.Evector(i,:) = Sm(i+1,2:end);
%     sr.Evector(i,:) = Sr(i+1,2:end);
%     sm.Evectorphi(i,:) = Sm(i+1,2:Ni+1);
%     sm.Evectorpsi(i,:) = Sm(i+1,Ni+2:end);
%     sr.Evectorphi(i,:) = Sr(i+1,2:Ni+1);
%     sr.Evectorpsi(i,:) = Sr(i+1,Ni+2:end);
%     sm.Evalue(i) = Sm(i+1,1);
%     sr.Evalue(i) = Sr(i+1,1);

    
    vm.Evector(i,:) = Vm(i+1,2:end);
    vr.Evector(i,:) = Vr(i+1,2:end);
    vm.Evalue(i) = Vm(i+1,1);
    vr.Evalue(i) = Vr(i+1,1);
%     
%     tm.Evector(i,:) = Tm(i+1,2:end);
%     tr.Evector(i,:) = Tr(i+1,2:end);
%     tm.Evalue(i) = Tm(i+1,1);
%     tr.Evalue(i) = Tr(i+1,1);
    end
    
    kt=Vm(1,2:end);
    
    %Interpolation Function
    tEq=150;
    
    t1=[3:3:9];
    t2=[10:5:40];
    t3=[40:30:150];
    t=[t1 t2 t3 300];
    t_teq = t/tEq;
       
    for l=1:size(t,2)
        a(l)=(((2^(1/2) - 1) * t_teq(l) + 1)^2 - 1);
        f_neil(l)=(1+a(l))^(-1);
    end
    
    
    %Pure Rad and Pure Mat ETCs
    ETCvRad=zeros(Ni,1);
    ETCvMat=ETCvRad;
    
    for i=1:Ni
      ETCvRad=ETCvRad+vr.Evalue(i)*vr.Evector(i,:)'.^2;
      ETCvMat=ETCvMat+vm.Evalue(i)*vm.Evector(i,:)'.^2;
    end
    
    figure()
    plot(kt,ETCvRad,'b',kt,ETCvMat,'k'); hold on;

    %Construct interpolated e-vectors
    
    for l=1:size(t,2)
    
        v.EvectorInt=[];
        
        for i=1:Ni
            if (vr.Evalue(i)<0)
                vr.Evalue(i)=0;
            end
            if (vm.Evalue(i)<0)
                vm.Evalue(i)=0;
            end
            v.EvectorInt(i,:) = f_neil(l)*sqrt(vr.Evalue(i))*vr.Evector(i,:) + (1-f_neil(l))*sqrt(vm.Evalue(i))*vm.Evector(i,:);
        end
    
        ETCtemp=zeros(Ni,1);
    
        for i=1:Ni
            ETCtemp=ETCtemp + (v.EvectorInt(i,:)').^2;
        end
    
        plot(kt,ETCtemp,'r'); hold on;
        set(gca, 'Xscale','log')
        whichkt=find(kt<150 & kt>50);
      
        for j=1:size(whichkt,2)
            f_t(l,j)=(ETCtemp(whichkt(j))-ETCvMat(whichkt(j)))/(ETCvRad(whichkt(j))-ETCvMat(whichkt(j)));
        end
        
    end
    
    %s = fitoptions('Method','NonlinearLeastSquares','Lower',[0.1,-10],'Upper',[0.6,0]);
    %f = fittype('(1 + a* x)^b','options',s);
    %[c2,gof2] = fit(time',funtzion',f)
    %[c2,gof2] = fit(t_teq',f_t',f)
    
    figure()
    %plot(c2,'k');hold on;
    for j=1:size(whichkt,2)
    scatter(t_teq(:),f_t(:,j),'k','filled');hold on;
    end
    plot(t_teq,f_neil,t_teq,f_neil.^2)
    set(gca, 'Xscale','log')
    
end
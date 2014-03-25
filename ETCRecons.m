function ETCRecons(path,Ni,id)
    b=id;

%     Sm=load([path 'UETCeigenScalarBeforeMat_' num2str(b,'%2.2d') '.dat']);
%     Sr=load([path 'UETCeigenScalarBeforeRad_' num2str(b,'%2.2d') '.dat']);
    Vm=load([path 'UETCeigenVectorBeforeMat_' num2str(b,'%2.2d') '.dat']);
    Vr=load([path 'UETCeigenVectorBeforeRad_' num2str(b,'%2.2d') '.dat']);
%      Tm=load([path 'UETCeigenTensorBeforeMat_' num2str(b,'%2.2d') '.dat']);
%      Tr=load([path 'UETCeigenTensorBeforeRad_' num2str(b,'%2.2d') '.dat']);
    
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
    tEq=50;
    
    t1=[3:3:9];
    t2=[10:5:40];
    t3=[40:30:150];
    t4=[150:50:1000];
    t=[t1 t2 t3 t4];
    t_teq = t/tEq;
       
    for l=1:size(t,2)
        a(l)=(((2^(1/2) - 1) * t_teq(l) + 1)^2 - 1);
        f_neil(l)=(1+a(l))^(-1);
        f_LAH(l)=(1+(1/4)*(t_teq(l)))^(-1/2);
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
    set(gca, 'Xscale','log')    
    
    
    %Cross Correlator
    
    ETCXrm=zeros(Ni,1);
    ETCXmr=zeros(Ni,1);
    for i=1:Ni
      ETCXrm=ETCXrm+(sqrt(vr.Evalue(i))*vr.Evector(i,:)').*(sqrt(vm.Evalue(i))*vm.Evector(i,:)');
      ETCXmr=ETCXmr+(sqrt(vm.Evalue(i))*vm.Evector(i,:)').*(sqrt(vr.Evalue(i))*vr.Evector(i,:)');

    end
    
    plot(kt,ETCXrm,'r'); hold on;
    set(gca, 'Xscale','log')
    
    
    %Construct interpolated e-vectors
    
    for l=1:size(t,2)
    
        v.EvectorInt=[];
        
        for i=1:Ni
%             if (vr.Evalue(i)<0)
%                 vr.Evalue(i)=0;
%             end
%             if (vm.Evalue(i)<0)
%                 vm.Evalue(i)=0;
%             end
            v.EvectorInt(i,:) = f_neil(l)*sqrt(vr.Evalue(i))*vr.Evector(i,:) + (1-f_neil(l))*sqrt(vm.Evalue(i))*vm.Evector(i,:);
            %v.EvectorInt(i,:) = f_LAH(l)*sqrt(vr.Evalue(i))*vr.Evector(i,:) + (1-f_LAH(l))*sqrt(vm.Evalue(i))*vm.Evector(i,:);

        end
    
        ETCtemp=zeros(Ni,1);
    
        for i=1:Ni
            ETCtemp=ETCtemp + (v.EvectorInt(i,:)').^2;
        end
    
%         ETCtemp = ETCtemp  - 2*f_LAH(l)*(1-f_LAH(l))*(ETCXrm-ETCvMat);
        ETCtemp = ETCtemp  - 2*f_neil(l)*(1-f_neil(l))*(ETCXrm-ETCvMat);
        
        plot(kt,ETCtemp,'r'); hold on;
        set(gca, 'Xscale','log')
        whichkt1=find(kt<51 & kt>50)
        whichkt2=find(kt<41 & kt>40)
        whichkt3=find(kt<32 & kt>30)
        whichkt4=find(kt<21 & kt>20)
        whichkt5=find(kt<11 & kt>10)
%         
      
        %for j=1:size(whichkt,2)
            %f_t(l,j)=(ETCtemp(whichkt(j))-ETCvMat(whichkt(j)))/(ETCvRad(whichkt(j))-ETCvMat(whichkt(j)));
        %end
        f_t(l)=(ETCtemp(whichkt1)-ETCvMat(whichkt1))/(ETCvRad(whichkt1)-ETCvMat(whichkt1));
        f_t2(l)=(ETCtemp(whichkt2)-ETCvMat(whichkt2))/(ETCvRad(whichkt2)-ETCvMat(whichkt2));
        f_t3(l)=(ETCtemp(whichkt3)-ETCvMat(whichkt3))/(ETCvRad(whichkt3)-ETCvMat(whichkt3));
        f_t4(l)=(ETCtemp(whichkt4)-ETCvMat(whichkt4))/(ETCvRad(whichkt4)-ETCvMat(whichkt4));
        f_t5(l)=(ETCtemp(whichkt5)-ETCvMat(whichkt5))/(ETCvRad(whichkt5)-ETCvMat(whichkt5));
%         
        mean_ft(l)=(f_t(l) + f_t2(l) + f_t3(l) + f_t4(l) + f_t5(l))/5;
    end
    
    for l=1:size(t,2)
        stdev_ft(l)= 0;
        stdev_ft(l) = sqrt((1/5)*(mean_ft(l) - f_t(l))^2 + (mean_ft(l) - f_t2(l))^2 + (mean_ft(l) - f_t3(l))^2 + (mean_ft(l) - f_t4(l))^2 + (mean_ft(l) - f_t5(l))^2);
    end
    
    s = fitoptions('Method','NonlinearLeastSquares','Lower',[0.1,-10],'Upper',[0.6,0]);
    f = fittype('(1 + a* x)^b','options',s);
    %[c2,gof2] = fit(time',funtzion',f)
    [c2,gof2] = fit(t_teq',mean_ft',f)
    
    time_teq=0.01:0.1:100;
    func_jo=(1+(1/4)*time_teq).^(-1);
    func_da=(1+(1/4)*time_teq).^(-2);
    
    figure()
    plot(c2,'k');hold on;
    %for j=1:size(whichkt,2)
    scatter(t_teq(:),mean_ft(:),'k','filled');hold on;
    plotErrorBars(t_teq,mean_ft,stdev_ft,'b');hold on;
    %end
    plot(t_teq,f_neil,'g',t_teq,f_neil.^2,'g');hold on;
    %plot(time_teq,func_jo,'r',time_teq,func_da,'k');hold on;
    set(gca, 'Xscale','log')
    
end
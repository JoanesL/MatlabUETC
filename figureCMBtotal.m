function figureCMBtotal(pathCell,idCell, Num, colour,legendCell, old, slope)


for i=1:numel(pathCell)
    if (old==1 && i==1)
    
    fileAll=[pathCell{i} 'extraCls-AH5-00.dat'];
    Cl=load(fileAll);
    
    l{i}=Cl(:,1);
    
    TT{i}=Cl(:,2);
    TE{i}=Cl(:,3);
    EE{i}=Cl(:,4);
    BB{i}=Cl(:,5);
    else
fileScalar{i}=[pathCell{i} 'SUMscalar_001-' int2str(Num(i)) '_'  idCell{i} '.dat'];
fileVector{i}=[pathCell{i} 'SUMvector_001-' int2str(Num(i)) '_'  idCell{i} '.dat'];
fileTensor{i}=[pathCell{i} 'SUMtensor_001-' int2str(Num(i)) '_'  idCell{i} '.dat'];

S=load(fileScalar{i});
V=load(fileVector{i});
T=load(fileTensor{i});

Sca{i}=S;
Vec{i}=V;
Ten{i}=T;

l{i}=S(:,1);

%TT{i}=V(:,2);
%TE{i}=V(:,3);
%EE{i}=V(:,4);
%BB{i}=V(:,5);

TT{i}=S(:,2)+V(:,2)+T(:,2);
TE{i}=S(:,3)+V(:,3)+T(:,3);
EE{i}=S(:,4)+V(:,4)+T(:,4);
BB{i}=V(:,5)+T(:,5);

    end
end


%Power-law Fit for TT, EE and BB
if slope==1
    TT1=TT{1};
    S1=Sca{1};
    V1=Vec{1};
    T1=Ten{1};
    
    EE1=EE{1};
    BB1=BB{1};
    
    l1=l{1};
for j=1:50
    fitlimit(j)=1200+(j-1)*50;
    zein{j}=find(l1>fitlimit(j));
    fitlimitEE(j)=2200+(j-1)*25;
    zeinEE{j}=find(l1>fitlimitEE(j));
end

for j=1:50
    
    TT_TOT = polyfit(log10(l1(zein{j})),log10(TT1(zein{j},:)),1);
    TTTOT(j)=TT_TOT(1);
    TT_Sca = polyfit(log10(l1(zein{j})),log10(S1(zein{j},2)),1);
    TTSca(j)=TT_Sca(1);
    TT_Vec = polyfit(log10(l1(zein{j})),log10(V1(zein{j},2)),1);
    TTVec(j)=TT_Vec(1);
    TT_Ten = polyfit(log10(l1(zein{j})),log10(T1(zein{j},2)),1);
    TTTen(j)=TT_Ten(1);
    
    EE_TOT = polyfit(log10(l1(zeinEE{j})),log10(EE1(zeinEE{j},:)),1);
    EETOT(j)=EE_TOT(1);
    BB_TOT = polyfit(log10(l1(zein{j})),log10(BB1(zein{j},:)),1);
    BBTOT(j)=BB_TOT(1);
    
    if j==1
        Fit_TT = TT_TOT(1)*log10(l1) + TT_TOT(2);
        Fit_TTSca = TT_Sca(1)*log10(l1) + TT_Sca(2);
        Fit_TTVec = TT_Vec(1)*log10(l1) + TT_Vec(2);
        Fit_TTTen = TT_Ten(1)*log10(l1) + TT_Ten(2);
        
        Fit_EE = EE_TOT(1)*log10(l1) + EE_TOT(2);
        Fit_BB = BB_TOT(1)*log10(l1) + BB_TOT(2);
    end
end

figure()
axi=multiPlot([2 2]);
axes(axi(1))
    
    plot(l1,TT1,'k','LineWidth',2); hold on;
    plotLogNegative(l1,10.^(Fit_TT),'b',0.5);
    legend('Baseline')%,'Baseline+Sudden Low kXi')
        plot(l1,S1(:,2),'-.k','LineWidth',1); hold on;
        plotLogNegative(l1,10.^(Fit_TTSca),'b',0.5);
        plot(l1,V1(:,2),'+k','LineWidth',1); hold on;
        plotLogNegative(l1,10.^(Fit_TTVec),'b',0.5);
        plot(l1,T1(:,2),'.k','LineWidth',1); hold on;
        plotLogNegative(l1,10.^(Fit_TTTen),'b',0.5);
        
        textTot=['Tot,  \gamma=' num2str(TTTOT(1)) ];
        textS=['S,  \gamma=' num2str(TTSca(1)) ];
        textV=['V,  \gamma=' num2str(TTVec(1)) ];
        textT=['T,  \gamma=' num2str(TTTen(1)) ];
text(2000,100,textTot,'Color','k','FontSize',14)        
text(1300,60,textS,'Color','k','FontSize',14)
text(1300,25,textV,'Color','k','FontSize',14)
text(1000,7,textT,'Color','k','FontSize',14)
    title('TT Decomp', 'FontWeight','bold')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'YLim',[5 600])
    set(gca,'XLim',[800 5050])
    xlabel('l')
ylabel('l(l+1)C_l / (G\mu)^2')
%set(gca,'YMinorTick','on','LineWidth',2)
set(gca,'FontSize',14)

axes(axi(2))

    plot(l1,TT1,'k','LineWidth',2); hold on;
    plotLogNegative(l1,10.^(Fit_TT),'b',0.5);
    plot(l1,EE1,'-.k','LineWidth',2); hold on;
    plotLogNegative(l1,10.^(Fit_EE),'b',0.5);
    plot(l1,BB1,'k','LineWidth',2); hold on;
    plotLogNegative(l1,10.^(Fit_BB),'b',0.5);
    
    textTT=['TT,  \gamma=' num2str(TTTOT(1)) ];
    textEE=['EE,  \gamma=' num2str(EETOT(1)) ];
    textBB=['BB,  \gamma=' num2str(BBTOT(1)) ];
text(2000,115,textTT,'Color','k','FontSize',14)        
text(3000,0.2,textEE,'Color','k','FontSize',14)
text(1300,0.01,textBB,'Color','k','FontSize',14)

    title('TT, EE and BB', 'FontWeight','bold')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('l')
ylabel('l(l+1)C_l / (G\mu)^2')
%set(gca,'YMinorTick','on','LineWidth',2)
set(gca,'FontSize',14)
set(gca,'XLim',[800 5050])
set(gca,'YLim',[1E-3 600])

axes(axi(3))
plot(fitlimit,TTTOT,'k',fitlimit,TTSca,'-.k',fitlimit,TTVec,'+k',fitlimit,TTTen,'.k')
xlabel('l')
ylabel('Evolution of the slope')
axes(axi(4))
plot(fitlimit,TTTOT,'k',fitlimitEE,EETOT,'-.k',fitlimit,BBTOT,'k')
xlabel('l')
ylabel('Evolution of the slope')



elseif slope==2 %Low-l of BB and EE
    EE_TOT = polyfit(log10(l1(l1<4)),log10(EE1(l1<4,:)),1);
    EETOT=EE_TOT(1);
    %BB_TOT = polyfit(log10(l1(l1<4)),log10(BB1(l1<4,:)),1);
    %BBTOT=BB_TOT(1);
    
    BBTOT=(log10(BB1(1,:)/BB1(2,:))/log10(l1(1)/l1(2)));
    
    Fit_BB=BBTOT*log10(l1)+((log10(BB1(2,:))-BBTOT*log10(l1(2))));
        
    Fit_EE = EE_TOT(1)*log10(l1) + EE_TOT(2);
    %Fit_BB = BB_TOT(1)*log10(l1) + BB_TOT(2);
    
    plot(l1,EE1,'-.k','LineWidth',2); hold on;
    plotLogNegative(l1,10.^(Fit_EE),'b',0.5);
    plot(l2,BB2,'k','LineWidth',2); hold on;
    %plotLogNegative(l1,10.^(Fit_BB),'b',0.5);
    
    textEE=['EE,  \gamma=' num2str(EETOT) ];
    textBB=['BB,  \gamma=' num2str(BBTOT) ];
text(2,0.002,textEE,'Color','k','FontSize',14)
%text(3,0.01,textBB,'Color','k','FontSize',14)
    
    title('EE and BB: Low l', 'FontWeight','bold')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('l')
    ylabel('l(l+1)C_l / (G\mu)^2')
    %set(gca,'YMinorTick','on','LineWidth',2)
    set(gca,'FontSize',14)
    set(gca,'XLim',[1.75 10])
    set(gca,'YLim',[1E-3 5E-3])
else
    
    
%TT10=100-((TT1(10)*100/TT2(10)));
%TT100=100-((TT1(100)*100/TT2(100)));
%TT1000=100-((TT1(1000)*100/TT2(1000)));

%TE10=100-((TE1(10)*100/TE2(10)));
%TE100=100-((TE1(100)*100/TE2(100)));
%TE1000=100-((TE1(1000)*100/TE2(1000)));

%EE10=100-((EE1(10)*100/EE2(10)));
%EE100=100-((EE1(100)*100/EE2(100)));
%EE1000=100-((EE1(1000)*100/EE2(1000)));

%BB10=100-((BB1(10)*100/BB2(10)));
%BB100=100-((BB1(100)*100/BB2(100)));
%BB1000=100-((BB1(1000)*100/BB2(1000)));

%disp('Differences with respect s=0 case')
%disp(['TT: l=10-> ' num2str(TT10) '  l=100-> ' num2str(TT100) '  l=1000-> ' num2str(TT1000) ])
%disp(['TE: l=10-> ' num2str(TE10) '  l=100-> ' num2str(TE100) '  l=1000-> ' num2str(TE1000) ])
%disp(['EE: l=10-> ' num2str(EE10) '  l=100-> ' num2str(EE100) '  l=1000-> ' num2str(EE1000) ])
%disp(['BB: l=10-> ' num2str(BB10) '  l=100-> ' num2str(BB100) '  l=1000-> ' num2str(BB1000) ])


figure()
ax=multiPlot([2 2]);
axes(ax(1))
    
    plot(l{1},TT{1},colour{1},'LineWidth',1.5); hold on;
    plot(l{2},TT{2},colour{2},'LineWidth',1.5); hold on;
    plot(l{3},TT{3},colour{3},'LineWidth',1.5); hold on;
    legend(legendCell{1},legendCell{2},legendCell{3})
%         plot(l{1},S(:,2),'LineWidth',3); hold on;
%         plot(l{1},V(:,2),'LineWidth',3); hold on;
%         plot(l{1},T(:,2),'LineWidth',3); hold on;

    title('Temperature PS', 'FontWeight','bold')
    set(gca,'XScale','log')
    %set(gca,'YScale','log')
    xlabel('l')
ylabel('l(l+1)C_l / (G\mu)^2')
%set(gca,'YMinorTick','on','LineWidth',2)
set(gca,'FontSize',14)

%text(800,150,'S','Color','k','FontSize',14)
%text(100,90,'V','Color','k','FontSize',14)
%text(150,30,'T','Color','k','FontSize',14)

axes(ax(2))
    
    plot(l{1},TE{1},colour{1},'LineWidth',1.5); hold on;
    plot(l{2},TE{2},colour{2},'LineWidth',1.5); hold on;
    plot(l{3},TE{3},colour{3},'LineWidth',1.5); hold on;
    
    title('TE', 'FontWeight','bold')
    set(gca,'XScale','log')
    %set(gca,'YScale','log')
    xlabel('l')
ylabel('l(l+1)C_l / (G\mu)^2')
%set(gca,'YMinorTick','on','LineWidth',2)
set(gca,'FontSize',14)

axes(ax(3))
    
    plot(l{1},EE{1},colour{1},'LineWidth',1.5); hold on;
    plot(l{2},EE{2},colour{2},'LineWidth',1.5); hold on;
    plot(l{3},EE{3},colour{3},'LineWidth',1.5); hold on;
    %plot(l{1},S(:,4),'k','LineWidth',3); hold on;
%         plot(l{1},V(:,4),'r','LineWidth',3); hold on;
%         plot(l{1},T(:,4),'b','LineWidth',3); hold on;
    title('EE', 'FontWeight','bold')
    set(gca,'XScale','log')
    %set(gca,'YScale','log')
    xlabel('l')
ylabel('l(l+1)C_l / (G\mu)^2')
%set(gca,'YMinorTick','on','LineWidth',2)
set(gca,'FontSize',14)

axes(ax(4))
    
    plot(l{1},BB{1},colour{1},'LineWidth',1.5); hold on;
    plot(l{2},BB{2},colour{2},'LineWidth',1.5); hold on;
    plot(l{3},BB{3},colour{3},'LineWidth',1.5); hold on;
%         plot(l{1},V(:,5),'r','LineWidth',3); hold on;
%         plot(l{1},T(:,5),'b','LineWidth',3); hold on;
    title('BB', 'FontWeight','bold')
    set(gca,'XScale','log')
    %set(gca,'YScale','log')
    xlabel('l')
ylabel('l(l+1)C_l / (G\mu)^2')
%set(gca,'YMinorTick','on','LineWidth',2)
set(gca,'FontSize',14)

multiPlotZoom(ax);
end
%elseif i==4
%    %axes(ax(2))
%    title('Polarization: EE', 'FontWeight','bold')
%    set(gca,'XScale','log')
%    set(gca,'YScale','log')
%    xlabel('l')
%ylabel('l(l+1)C_l / (G\mu)^2')
%set(gca,'XScale','log')
%set(gca,'YScale','log')
%set(gca,'YMinorTick','on','LineWidth',2)
%set(gca,'FontSize',14)

%elseif i==5
    %axes(ax(3))

%    title('Polarization: BB', 'FontWeight','bold')
%    set(gca,'XScale','log')
%    set(gca,'YScale','log')
%    xlabel('l')
%ylabel('l(l+1)C_l / (G\mu)^2')
%set(gca,'XScale','log')
%set(gca,'YScale','log')
%set(gca,'YMinorTick','on','LineWidth',2)
%set(gca,'FontSize',14)

%end
%end

end


%set(gca,'XLim',[2 5000],'YLim',[0 5])
%figure(2); clf
%plot(l,l.^2.*(S(:,2)+V(:,2)+T(:,2))./(l+1),'k','LineWidth',2); hold on
%plot(l,l.^2.*S(:,2)./(l+1),'k'); hold on
%plot(l,l.^2.*V(:,2)./(l+1),'k')
%plot(l,l.^2.*T(:,2)./(l+1),'k')

%text(1000,1.5e5,'S','FontSize',14)
%text(1000,3e4,'V','FontSize',14)
%text(1000,9e3,'T','FontSize',14)

%betaS=2.5;
%betaV=1.9;
%betaT=1.7;

%le=[4000:100:6000];
%plot(le,le.^2.*(S(end,2).*(4000./le).^(betaS-1)+V(end,2).*(4000./le).^(betaV-1)+T(end,2).*(4000./le).^(betaT-1))./(le+1),'--k','LineWidth',2)
%plot(le,le.^2.*(S(end,2).*(4000./le).^(betaS-1))./(le+1),'--k')
%plot(le,le.^2.*(V(end,2).*(4000./le).^(betaV-1))./(le+1),'--k')
%plot(le,le.^2.*(T(end,2).*(4000./le).^(betaT-1))./(le+1),'--k')

%set(gca,'XScale','log')
%set(gca,'YScale','log')
%set(gca,'YMinorTick','on','LineWidth',2)
%set(gca,'FontSize',14)
%set(gca,'XLim',[200 6000])
%set(gca,'YLim',[3e3 3e5])

%plot([4000 4000],get(gca,'YLim'),'-','Color',[0.5 0.5 0.5])

%xlabel('l')
%ylabel('l^3 C_l / (G\mu)^2')

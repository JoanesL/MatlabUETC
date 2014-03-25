function BICEP2spec(path)

%fileTensor=[path 'LCDMr015_cl_lensed.dat'];
fileTensor=[path 'LCDMr015K0002_cl_lensed.dat'];
%fileTensorBICEP=[path 'LCDMr02_cl_lensed.dat'];
fileTensorBICEP=[path 'LCDMr02K0002_cl_lensed.dat'];
%fileTensorBICEP=[path 'B2_3yr_camb_planck_withB_uK_20140314.txt'];
fileLens=[path 'LCDMr00_cl_lensed.dat'];
fileStrings=[path 'maybe_AH5_planck/extraCls_AH5_Planck.dat'];
%fileStrings=[path 'SUM_001-256_01.dat'];
fileSL=[path 'seminew.dat'];
fileTex=[path 'SLAndTex/CMBdataTX0/SUM_001-128_00.dat'];
fileTexScalar=[path 'maybe_AH5_planck/SUMscalar_001-512_14.dat'];
fileTexVector=[path 'maybe_AH5_planck/SUMvector_001-512_14.dat'];
fileTexTensor=[path 'maybe_AH5_planck/SUMtensor_001-512_14.dat'];
fileBICEP=[path 'BICEP.txt'];
filePlHigh=[path 'Planck_highl.txt'];
filePlLow=[path 'Planck_lowl.txt'];


Tens=load(fileTensor);
TensBICEP=load(fileTensorBICEP);
Lens=load(fileLens);
Str=load(fileStrings);
SL=load(fileSL);
Tex=load(fileTex);
TexS=load(fileTexScalar);
TexV=load(fileTexVector);
TexT=load(fileTexTensor);
BIC=load(fileBICEP);

PLH=load(filePlHigh);
PLL=load(filePlLow);

lPll=PLL(:,1)


TTPLlow=PLL(:,2);
TTPLdown_temp=PLL(:,3);
TTPLup_temp=PLL(:,4);
TTPLdown=TTPLlow-TTPLdown_temp;
TTPLup=TTPLup_temp-TTPLlow;

lPlhLeft_temp=PLH(:,2);
lPlh=PLH(:,1)
lPlhRight_temp=PLH(:,3);
lPlhLeft=lPlh-lPlhLeft_temp;
lPlhRight=lPlhRight_temp-lPlh;

TTPLhigh=PLH(:,4);
TTPLhStDev=PLH(:,5);

lLeft_temp=BIC(:,1);
lCenter=BIC(:,2);
lRigth_temp=BIC(:,3);
lLeft=lCenter-lLeft_temp;
lRigth=lRigth_temp-lCenter;
llp1=lCenter(:).*(lCenter(:)+1);
TTBIC=BIC(:,4);%./llp1*2*pi;
TTSdev=BIC(:,10)%./llp1*2*pi;
BBBIC=BIC(:,7);%./llp1
BBSdev=BIC(:,13);%:llp1;

T_cmb=2.726;
Factor=(T_cmb^2)/(2*pi);

lTens=Tens(:,1);
TTTens=Tens(:,2);%*Factor;
BBTens=Tens(:,4);%*Factor;

lTensBICEP=TensBICEP(:,1);
TTTensBICEP=TensBICEP(:,2);
BBTensBICEP=TensBICEP(:,4);

lLens=Lens(:,1);
TTLens=Lens(:,2);
BBLens=Lens(:,4);

lStr=Str(:,1);
llp1=lStr(:).*(lStr(:)+1);
TTStr=Str(:,2);%*Factor*(10^-12)/(2*pi);
BBStr=Str(:,5);%*Factor*(10^-12)/(2*pi);

lSL=SL(:,1);
TTSL=SL(:,2);%*Factor;
BBSL=SL(:,5);%*Factor;

lTex=TexV(:,1);
TTTex=Tex(:,2);%*Factor;
BBTex=Tex(:,4);%*Factor;

%BBTexS=TexS(:,3);
BBTexV=TexV(:,5);
BBTexT=TexT(:,5);

plot(lTex,BBTexV,'b',lTex,BBTexT,'g')

%at l=80
BicepFactor1=BBTens(79)/(BBStr(79));%+BBLens(79))
BicepFactor2=BBTens(79)/BBSL(79);
BicepFactor3=BBTens(79)/BBTex(79);


BBStr_int=interp1(lStr,BBStr,lTens);
TTStr_int=interp1(lStr,TTStr,lTens);


%Fig 2
%Normalize String, textures and SL to equate Planck data at l=10
PlanckFactor1=TTPLlow(9)/TTStr(9);
PlanckFactor2=TTPLlow(9)/TTSL(9);
PlanckFactor3=TTPLlow(9)/TTTex(9);


%Fig 3
%f_10
% beta=0.2628;
% beta2=0.03;
% beta3=0.06;
% beta4=0.15;

% %Fig 4
beta=0.2628;
beta2=0.03;
beta3=0.06;
beta4=0.015;

%alpha -> factor to get desired f_10's
%alpha= (beta * TTTens(9))/((1-beta)*TTStr_int(9)); %Normalized at l=80 (beta * TTTens(9))/((1-beta)*TTStr_int(9)*BicepFactor1)
alphaBF = 0.166 * Factor;
alpha = (beta * TTTens(9))/((1-beta)*TTStr_int(9))
alpha2= (beta2 * TTTens(9))/((1-beta2)*TTStr_int(9));%*BicepFactor1)
alpha3= (beta3 * TTTens(9))/((1-beta3)*TTStr_int(9));%*BicepFactor1)
alpha4= (beta4 * TTTens(9))/((1-beta4)*TTStr_int(9));%*BicepFactor1)

f_10=(TTStr_int(9)*alphaBF)./(TTStr_int(9)*alphaBF+TTTens(9))
f_10=(TTStr_int(9)*alpha2)./(TTStr_int(9)*alpha2+TTTens(9))
f_10=(TTStr_int(9)*alpha3)./(TTStr_int(9)*alpha3+TTTens(9))

factorPlanck=1.6297*10^(-14);

% figure()
% plot(lStr,BBStr*PlanckFactor1,'k',lSL,BBSL*PlanckFactor2,'r',lTex,BBTex*PlanckFactor3,'-.k','LineWidth',3); hold on;%lTens,BBStr_int*factorPlanck,'-.k',lSL,BBSL*BicepFactor2,'r',lTex,BBTex*BicepFactor3,'-.k');hold on;
% 
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
% set(gca,'XLim',[10,2000])
% set(gca,'YLim',[10^-3 3])
% set(gca,'TickLength',get(gca,'TickLength')*1.5)
% set(gca,'LineWidth',1.1,'FontSize',20)
%ax=multiPlot([1 2]);


%############################
%Fig 2
%############################
% figure()
% plot(lTens,BBStr_int*PlanckFactor1,'-.b','LineWidth',2);hold on;
% plot(lSL,BBSL*PlanckFactor2,'r','LineWidth',2); hold on;
% plot(lTex,BBTex*PlanckFactor3,'--k','LineWidth',2);hold on;
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
% set(gca,'XLim',[10,2000])
% set(gca,'YLim',[0.001 3])
% set(gca,'TickLength',get(gca,'TickLength')*1.5)
% set(gca,'LineWidth',1.1,'FontSize',36)
% ylabel('l(l+1)C^{BB}_l/2\pi','fontsize',36)
% xlabel('l','fontsize',36)

grey=[0.5 0.5 0.5];

%############################
%Fig 3
%############################
%Upper panel
%############################
% %High l Planck data
% %LCDM + r=0.2 + lensing
% figure()
% plot(lTens,TTTens,'k','LineWidth',2);hold on;
% %String:
% plot(lLens,TTStr_int*alpha+TTTens,'-.','Color',grey,'LineWidth',1.5); hold on;
% plot(lLens,TTStr_int*alpha,'-.b','LineWidth',1.5); hold on; %Normalized at l=80
% plot(lLens,TTStr_int*alpha2,'-.b','LineWidth',1.5); hold on;
% plot(lLens,TTStr_int*alpha3,'-.b','LineWidth',1.5); hold on;
% plot(lLens,TTStr_int*alpha4,'-.b','LineWidth',1.5); hold on;
% plotErrorBars(lPlh,TTPLhigh,TTPLhStDev,'k');hold on;
% for i=1:size(lPlh)
%     herrorbar(lPlh(i), TTPLhigh(i), lPlhLeft(i), lPlhRight(i));hold on;
% end
% %errorbar(lPll,TTPLlow,TTPLdown,TTPLup,'ok','LineWidth',2)
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
% set(gca,'XLim',[10,2000])
% set(gca,'YLim',[10 7000])
% set(gca,'TickLength',get(gca,'TickLength')*1.5)
% set(gca,'LineWidth',1.1,'FontSize',36)
% ylabel('l(l+1)C^{TT}_l/2\pi [\muK^2]','fontsize',36)
% xlabel('l','fontsize',36)
% 
% %############################
% %Lower panel
% %############################
% green=[0 0.5 0];
% figure()
% plot(lLens,BBStr_int*alpha+BBLens,'-.','Color',green,'LineWidth',2); hold on; %Normalized at l=80
% plot(lLens,BBStr_int*alpha2+BBLens,'-.','Color',green,'LineWidth',2); hold on;
% plot(lLens,BBStr_int*alpha3+BBLens,'-.','Color',green,'LineWidth',2); hold on;
% plot(lLens,BBStr_int*alpha4+BBLens,'-.','Color',green,'LineWidth',2); hold on;
% plot(lTensBICEP,BBTensBICEP,'k','LineWidth',3);hold on;
% % %plot(lTens,BBTens,'k','LineWidth',
% % %plot(lTens,BBStr_int*BicepFactor1*alpha+BBTens,'g')
%  plotErrorBars(lCenter,BBBIC,BBSdev,'k');hold on;
%  for i=1:size(lCenter)
%      herrorbar(lCenter(i), BBBIC(i), lLeft(i), lRigth(i));hold on;
% end
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
% set(gca,'XLim',[10,2000])
% set(gca,'YLim',[10^-3 0.3])
% set(gca,'TickLength',get(gca,'TickLength')*1.5)
% set(gca,'LineWidth',1.1,'FontSize',36)
% ylabel('l(l+1)C^{BB}_l/2\pi [\muK^2]','fontsize',36)
% xlabel('l','fontsize',36)
% set(gca,'Xscale','log')
% set(gca,'XLim',[10,2000])


%############################
% Fig 4
%############################

% figure()
% %plot(lTens,BBStr_int*alphaBF+BBTens,'-.b','LineWidth',2); hold on; %Normalized at l=80
% plot(lTens,BBStr_int*alpha2+BBTens,'Color',grey,'LineWidth',2); hold on;
% plot(lTens,BBStr_int*alpha3+BBTens,'Color',grey,'LineWidth',2); hold on;
% plot(lTens,BBStr_int*alpha4+BBTens,'Color',grey,'LineWidth',2); hold on;
% plot(lTens,BBStr_int*alphaBF+BBTens,'r','LineWidth',2);hold on;
% 
% plot(lTens,BBStr_int*alpha2,'-.b','LineWidth',2); hold on;
% plot(lTens,BBStr_int*alpha3,'-.b','LineWidth',2); hold on;
% plot(lTens,BBStr_int*alpha4,'-.b','LineWidth',2); hold on;
% plot(lTens,BBStr_int*alphaBF,'r','LineWidth',2);hold on;
% 
% plot(lTens,BBTens,'k','LineWidth',2);hold on;
% % %plot(lTens,BBTens,'k','LineWidth',
% % plot(lTens,BBStr_int,'-.m',lSL,BBSL,'r',lTex,BBTex,'-.k');hold on;
% % %plot(lTens,BBStr_int*BicepFactor1*alpha+BBTens,'g')
%  plotErrorBars(lCenter,BBBIC,BBSdev,'k');hold on;
%  for i=1:size(lCenter)
%      herrorbar(lCenter(i), BBBIC(i), lLeft(i), lRigth(i));hold on;
% end
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
% set(gca,'XLim',[10,2000])
% set(gca,'YLim',[10^-3 0.3])
% set(gca,'TickLength',get(gca,'TickLength')*1.5)
% set(gca,'LineWidth',1.1,'FontSize',36)
% ylabel('l(l+1)C^{BB}_l/2\pi [\muK^2]','fontsize',36)
% xlabel('l','fontsize',36)


function f=CMBpolarone(path,n,m,f10,color,pos,num)

dir1=path;
started=0;
% for i=n
%    for imag=['r' 'i']
%     fileS=[dir1 'CMBscalar_' num2str(i,'%3.3i') imag '_' num2str(m,'%2.2d' ) '.dat'];
%     fileV=[dir1 'CMBvector_' num2str(i,'%3.3i') imag '_' num2str(m,'%2.2d' ) '.dat'];
%     fileT=[dir1 'CMBtensor_' num2str(i,'%3.3i') imag '_' num2str(m,'%2.2d' ) '.dat'];
% 
%     disp(['Loading file: ' fileS])
%     CMBeigenS=load(fileS);
%     disp(['Loading file: ' fileV])
%     CMBeigenV=load(fileV);
%     disp(['Loading file: ' fileT])
%     CMBeigenT=load(fileT);
% 
%     if started==0
%        started=1;
%         CMBs=CMBeigenS;
%         CMBv=CMBeigenV;
%         CMBt=CMBeigenT;
%     else
%         CMBs(:,2:end)=CMBs(:,2:end)+CMBeigenS(:,2:end);
%         CMBv(:,2:end)=CMBv(:,2:end)+CMBeigenV(:,2:end);
%         CMBt(:,2:end)=CMBt(:,2:end)+CMBeigenT(:,2:end);
%     end
%  end
% end
% 
% l=CMBs(:,1);
% TTs=CMBs(:,2); EEs=CMBs(:,3); TEs=CMBs(:,4);
% TTv=CMBv(:,2); EEv=CMBv(:,3); BBv=CMBv(:,4); TEv=CMBv(:,5);
% TTt=CMBt(:,2); EEt=CMBt(:,3); BBt=CMBt(:,4); TEt=CMBt(:,5);
% 
% TT=TTs+TTv+TTt;
% EE=EEs+EEv+EEt;
% BB=BBv+BBt;
% TE=TEs-TEv-TEt;

% WMAP=804.1740; %WMAP-3 TT: l(l+1)Cl/(2pi) [micro-K^2]
% norm=WMAP/TT(8)*f10
% 
% TT=TT*norm; TTs=TTs*norm; TTv=TTv*norm; TTt=TTt*norm;
% EE=EE*norm; EEs=EEs*norm; EEv=EEv*norm; EEt=EEt*norm;
% BB=BB*norm; BBv=BBv*norm; BBt=BBt*norm;
% TE=TE*norm; TEs=TEs*norm; TEv=TEv*norm; TEt=TEt*norm;

fileStr=[dir1 'maybe_AH5_planck/extraCls_AH5_Planck.dat'];
filer02=[dir1 'LCDMr00_cl_lensed.dat'];
filer02Tens=[dir1 'LCDMr02_clt.dat'];

T_cmb=2.726;
Factor=(T_cmb*10^6)^2;
%Planck best fit: (Gmu)**2/(2pi)
factorPlanck=1.6297*10^(-14);


Str=load(fileStr);
r02=load(filer02);
r02Tens=load(filer02Tens);

lStr=Str(:,1);
TTStr=Str(:,2)*Factor;
BBStr=Str(:,5)*Factor;


l=r02(:,1);
TTr02=r02(:,2);
BBr02=r02(:,4);

lTens=r02Tens(:,1);
TTr02Tens=r02Tens(:,2);
BBr02Tens=r02Tens(:,4);

if pos==1
plot(lStr,TTStr*factorPlanck,'-.b',l,TTr02,'k',lTens,TTr02Tens,'--k','LineWidth',2);
%ylabel(' l(l+1) C_l^E^E / 2\pi [\muK^2]')
if num==0
set(gca,'YLim',[0 7000])
set(gca,'Yscale','log')
% set(gca,'Yscale','log')
%      annotation(gcf,'textbox','String',{'TT'},'FontSize',24,...
% 		'FitHeightToText','off',...
% 		'LineStyle','none',...
% 		'Position',[0.1627 0.8493 0.07683 0.05083]);
end
else
plot(lStr,BBStr*factorPlanck,'-.b',l,BBr02,'k',lTens,BBr02Tens,'--k','LineWidth',2);
%ylabel(' l(l+1) C_l^B^B / 2\pi [\muK^2]')
if num==0
% xlabel('l')
%      annotation(gcf,'textbox','String',{'BB'},'FontSize',24,...
% 		'FitHeightToText','off',...
% 		'LineStyle','none',...
% 		'Position',[0.1627 0.4066 0.07683 0.05083]);
% set(gca,'YLim',[2e-5 2],'YTick',[0.0001 0.001 0.01 0.1 1])
%      annotation(gcf,'textbox',...
% 		'String',{'l(l+1) C_l / 2\pi [\muK^2]'},...
% 		'FontSize',20,...
% 		'FitHeightToText','off',...
% 		'LineStyle','none',...
% 		'Rotation',90,......
% 		'Position',[0.01308 0.26 0.1548 0.1169]);
xlabel('l')
%ylog
set(gca,'Yscale','log')
end
end

set(gca,'XLim',[2 3000])

set(gca,'TickLength',get(gca,'TickLength')*1.5)
set(gca,'XScale','log','Xtick',[10 100 1000])
set(gca,'LineWidth',2,'FontSize',20)
hold on



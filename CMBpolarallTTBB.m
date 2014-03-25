
figure

clf
for pos=1:2

a2=axes;
p=get(a2,'Position');
if pos==1
    set(a2,'Position',[p(1) p(2)+p(4)/2 p(3) p(4)/2]);
else
set (a2,'Position',[p(1) p(2) p(3) p(4)/2]);
end


CMBpolaroneTTBB('../defect_spectra/',1:128,30,0.1,'--k',pos,1);
CMBpolaroneTTBB('../defect_spectra/',1:128,0,0.14,'r',pos,1);
CMBpolaroneTTBB('../defect_spectra/',1:128,0,0.16,'-.b',pos,0);

%CMBpolaroneTTBB('/net/rebus.scratch/nab21/CMBdataRealz/',1:28,30,0.11,'--k',pos,1);
%CMBpolaroneTTBB('/scratch/kap10/CMBdata/CMBdatas3/',1:28,0,0.15,'r',pos,1);
%CMBpolaroneTTBB('/scratch/kap10/CMBdata/CMBdataTX3/',1:28,0,0.16,'-.b',pos,0);


set(a2,'Color','none')
set(a2,'Box','on')
%set(a2,'FontSize',14);
set(a2,'LineWidth',2)

if pos==1; set(a2,'XTickLabel',''); end
if pos==2
    axes('Visible','off','Position',[0 0 1 1])
end

end


0.11
0.15
0.16

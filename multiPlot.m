%MULTIPLOT
%
%Creates a number of axes in the current figure
%in a similar way to the in-built subplot command.
%However, unlike subplot, the axis-labels and titles
%will not interfer with other axes. Also they are set to
%scale with the figure size, so for example, that 
%when a figure is exported (which is usually as if 
%it were maximised) the axis labels do not look tiny.
%
%Syntax: handles=MULTIPLOT(plots)
%
%The input plots is a two element array consisting of
%the number of plots horizontally followed by the number
%vertically. The output is a matrix consisting of the
%handles of each of the axes created such that handles(1,1)
%is the handle of the top left axis, handles(2,1) is the
%one to the right of that and handles(2,2) is the one below
%this one.
%
%To switch to a particular axis, use eg: axes(handles(1,1)),
%or axes(handles(3)) to switch to the third plot counting
%from the top left.

function handles=MultiPlot(plots)

%CHECK INPUT=======================================
if max(size(plots))~=2
   error('Input array must be a two element array');
end
   
%GET FIGURE DATA===================================
figure_size=get(gcf,'Position');
figure_size=figure_size(3:4);

%GET LABELS AND TITLE BOARDERS=====================
axes('Units','pixels'); 
xlabel('fgfg'); ylabel('fhggfh'); title('fgfgf');

position=get(gca,'Position');

set(get(gca,'XLabel'),'Units','pixels');
EXTENT=get(get(gca,'XLabel'),'Extent'); 
xlabel_boarder=-EXTENT(2);

set(get(gca,'YLabel'),'Units','pixels');
EXTENT=get(get(gca,'YLabel'),'Extent'); 
ylabel_boarder=-EXTENT(1);

set(get(gca,'Title'),'Units','pixels');
EXTENT=get(get(gca,'Title'),'Extent');
title_boarder=EXTENT(2)+EXTENT(4)-position(4);

delete(gca)

%ESTABLISH SIZE OF EACH AXIS WITH BOARDERS=========
axis_size(1)=(figure_size(1)-ylabel_boarder)/plots(1);
axis_size(2)=(figure_size(2)-xlabel_boarder)/plots(2);

%CREATE AXES=======================================
handles=zeros(plots);
for i=1:plots(1)
   for j=1:plots(2)
      handles(i,j)=axes('Units','pixels','Position',[(i-1)*axis_size(1),(plots(2)-j)*axis_size(2),axis_size]+[ylabel_boarder,xlabel_boarder,-ylabel_boarder,-xlabel_boarder-title_boarder]+[ylabel_boarder/3,xlabel_boarder/2,0,0]);  
      set(handles(i,j),'Units','normalized','FontUnits','normalized')   
      set(get(handles(i,j),'XLabel'),'Units','normalized','FontUnits','normalized')
      set(get(handles(i,j),'YLabel'),'Units','normalized','FontUnits','normalized')
      set(get(handles(i,j),'Title'),'Units','normalized','FontUnits','normalized')
   end	
end
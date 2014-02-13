%MULITPLOTZOOM
%
%Zoom in on a subplot in a new figure window

function multiPlotZoom(handles)

for i=1:prod(size(handles))

    handle=handles(i);
    
    if isempty(get(handle,'ButtonDownFcn'))==1   
        
        %SET BUTTONDOWN FUNCTION=================
        set(handle,'ButtonDownFcn',['multiPlotZoom(' num2str(handle,'%2.16f') ')'])
        
    else   
        
        %ZOOM IN ON SUBPLOT======================   
        f=figure;
        a=copyobj(handle,f);
        set(a,'ButtonDownFcn',['close(' num2str(f) ')'])
        b=axes;
        set(a,'Position',get(b,'Position'))
        delete(b);
        
    end
    
end
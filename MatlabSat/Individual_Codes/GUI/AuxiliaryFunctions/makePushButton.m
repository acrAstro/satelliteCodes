function ui = makePushButton(parent,value,left,top,callback,width,height)

if nargin < 6
    width = 30;
end

if nargin < 7
    height = 20;
end

parentHeight = getHeight(parent);

ui = uicontrol('parent',parent,'style','pushbutton','string',value,...
    'position',[left,parentHeight-height-top,width,height],'fontweight',...
    'bold','Backgroundcolor','green','callback',callback);

end
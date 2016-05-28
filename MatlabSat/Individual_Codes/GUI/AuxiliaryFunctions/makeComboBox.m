function ui = makeComboBox(parent,value,left,top,callback,width,height)

if nargin < 5
    callback = @nothing;
end

if nargin < 6
    width = 150;
end

if nargin < 7
    height = 20;
end

parentHeight = getHeight(parent);

ui = uicontrol('parent',parent,'style','popupmenu','string',value,'Position',[left,parentHeight-height-top,width,height],'backgroundcolor',[0.74,0.74,0.74],'fontweight','bold','callback',callback);
end

function nothing(~,~)
end
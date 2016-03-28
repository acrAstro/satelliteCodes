function ui = makeCheckbox(parent,value,left,top,callback,width,height)

if nargin < 5
    callback = @nothing;
end
if nargin < 6
    width = 15;
end
if nargin < 7
    height = 15;
end

parentHeight = getHeight(parent);

ui = uicontrol('parent',parent,'style','checkbox','value',value,'position',[left,parentHeight-height-top,width,height],'callback',callback);
end

function nothing(~,~)
end
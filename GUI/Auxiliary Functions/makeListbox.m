function ui = makeListbox(parent,left,top,callback,width,height)

if nargin < 4
    callback = @nothing;
end
if nargin < 5
    width = 200;
end
if nargin < 6
    height = 250;
end

parentHeight = getHeight(parent);

ui = uicontrol('Parent',parent,'style','listbox','Position',[left parentHeight-height-top,width,height],'fontweight','bold','callback',callback);
end

function nothing(~,~)
end
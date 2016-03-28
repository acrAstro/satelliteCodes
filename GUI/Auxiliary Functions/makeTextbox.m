function ui = makeTextbox(parent,value,left,top,callback,width,height)

if nargin < 6
    width = 100;
end

if nargin < 7
    height = 15;
end

parentHeight = getHeight(parent);

ui = uicontrol('parent',parent,'style','edit','string',value,'position',[left,parentHeight-height-top,width,height],'callback',callback);

end
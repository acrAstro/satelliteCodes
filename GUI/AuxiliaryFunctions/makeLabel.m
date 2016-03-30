function ui = makeLabel(parent,value,left,top,width,height)

if nargin < 5
    width = 160;
end

if nargin < 6
    height = 15;
end

GL = {'units','centimeters'};
HL = {'Horizontalalignment','left'};

parentHeight = getHeight(parent);

ui = uicontrol('Parent',parent,'style','text','string',value,...
    'Position',[left,parentHeight-height-top,width,height],'fontweight',...
    'bold',GL{1},GL{2},HL{1},HL{2});
end
function ui = makePanel(parent,title,left,top,width,height)

if nargin < 5
    width = 0.2;
end

if nargin < 6
    height = 0.2;
end

bottom = 1 - top - height;

ui = uipanel('Parent',parent,'title',title,'titleposition','lefttop','position',[left,bottom,width,height]);
end
function makeComboBox(parent,name,width,height,callback)

pHeight = parent.Position(4);
pWidth = parent.Position(3);

uicontrol('parent',parent,'String',name,'HorizontalAlignment','left','width',width,'height',height,'callback',callback);
end
function RelativeMotionGUI

clear; close all; clc;

guiHeight = 500;
guiWidth = 1000;
leftBuffer = 50;
topBuffer = 200;
screenSize = get(groot,'Screensize');
screenHeight = screenSize(4);
f = figure('Visible','off','position',[leftBuffer,screenHeight-guiHeight-topBuffer,guiWidth,guiHeight]);



f.Visible = 'on';

end
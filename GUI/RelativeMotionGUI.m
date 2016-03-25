function RelativeMotionGUI

clear; close all; clc;
addpath([pwd,filesep,'Auxiliary Functions']);
guiHeight = 500;
guiWidth = 1000;
leftBuffer = 50;
topBuffer = 200;
screenSize = get(groot,'Screensize');
screenHeight = screenSize(4);
f = figure('Visible','off','position',[leftBuffer,screenHeight-guiHeight-topBuffer,guiWidth,guiHeight]);

chiefPanel = makePanel(f,'Chief Orbit',0.01,0.01,0.3,0.3);

f.Visible = 'on';

end
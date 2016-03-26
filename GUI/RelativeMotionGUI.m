function RelativeMotionGUI

clear; close all; clc;
addpath([pwd,filesep,'Auxiliary Functions']);
guiHeight = 500;
guiWidth = 1000;
leftBuffer = 50;
topBuffer = 200;
screenSize = get(groot,'Screensize');
screenHeight = screenSize(4);
f = figure('Visible','on','position',[leftBuffer,screenHeight-guiHeight-topBuffer,guiWidth,guiHeight]);

chiefPanel = makePanel(f,'Chief Orbit',0.01,0.01,0.3,0.3);

chiefElemsLabels = {};
chiefElemsLabels{1,1} = 'Semi-major Axis, a';               chiefElemsLabels{1,2} = 'km';
chiefElemsLabels{2,1} = 'Eccentricity, e';                  chiefElemsLabels{2,2} = '';
chiefElemsLabels{3,1} = 'Inclination, i';                   chiefElemsLabels{3,2} = 'deg';
chiefElemsLabels{4,1} = 'Right-Ascension,O';         chiefElemsLabels{4,2} = 'deg';
chiefElemsLabels{5,1} = 'Argument of Perigee, w';    chiefElemsLabels{5,2} = 'deg';
chiefElemsLabels{6,1} = 'Mean Anomaly, M_0';                chiefElemsLabels{6,2} = 'deg';

for ii = 1:6
    makeLabel(chiefPanel,chiefElemsLabels{ii,1},10,20+ii);
    makeLabel(chiefPanel,chiefElemsLabels{ii,2},200,20+ii);
end

f.Visible = 'on';





end
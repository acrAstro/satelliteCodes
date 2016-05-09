function RelativeMotionGUI

clear; close all; clc;
addpath([pwd,filesep,'Auxiliary Functions']);
addpath([pwd,filesep,'OrbitFunctions']);
global Req;
Req = 6378;
global mu;
mu = 3.986e5;
global J2;
J2 = 1082.6e-6;

guiHeight = 500;
guiWidth = 1000;
leftBuffer = 50;
topBuffer = 200;
screenSize = get(groot,'Screensize');
screenHeight = screenSize(4);
f = figure('Visible','on','position',[leftBuffer,screenHeight-guiHeight-topBuffer,guiWidth,guiHeight]);

chiefPanel = makePanel(f,'Chief Orbit',0.01,0.01,0.3,0.3);

defChiefSMA = 6678;
defChiefEcc = 0.005;
defChiefInc = 28.5;
defChiefRaan = 45;
defChiefArgPer = 45;
defChiefTA = 0;
chiefElemsLabels = {};
chiefElemsLabels{1,1} = 'Semi-major Axis, a:';               chiefElemsLabels{1,2} = 'km';
chiefElemsLabels{2,1} = 'Eccentricity, e:';                  chiefElemsLabels{2,2} = '';
chiefElemsLabels{3,1} = 'Inclination, i:';                   chiefElemsLabels{3,2} = 'deg';
chiefElemsLabels{4,1} = 'Right-Ascension, O:';               chiefElemsLabels{4,2} = 'deg';
chiefElemsLabels{5,1} = 'Argument of Perigee, w:';           chiefElemsLabels{5,2} = 'deg';
chiefElemsLabels{6,1} = 'True Anomaly, M_0:';                chiefElemsLabels{6,2} = 'deg';

for ii = 1:6
    makeLabel(chiefPanel,chiefElemsLabels{ii,1}, 10,10+15*ii);
    makeLabel(chiefPanel,chiefElemsLabels{ii,2},200,10+15*ii);
end
chiefOrbitEditSMA = makeTextbox(chiefPanel,defChiefSMA,145,25,@nothing,50);
chiefOrbitEditEcc = makeTextbox(chiefPanel,defChiefEcc,145,40,@nothing,50);
chiefOrbitEditInc = makeTextbox(chiefPanel,defChiefInc,145,55,@nothing,50);
chiefOrbitEditRaan = makeTextbox(chiefPanel,defChiefRaan,145,70,@nothing,50);
chiefOrbitEditArgPer = makeTextbox(chiefPanel,defChiefArgPer,145,85,@nothing,50);
chiefOrbitEditTA = makeTextbox(chiefPanel,defChiefTA,145,100,@nothing,50);

f.Visible = 'on';



    function nothing(~,~)
    end

end

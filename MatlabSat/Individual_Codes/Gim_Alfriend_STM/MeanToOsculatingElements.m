function [DJ2,osc_c] = MeanToOsculatingElements(J2,meanElems,Re,mu)

% Transformation matrix D_J2 in closed form
% between mean and osculating new set of elements
% with the perturbation by only J2
%
% input :
%    mean_c(1) = a_mean
%    mean_c(2) = theta_mean
%    mean_c(3) = i_mean
%    mean_c(4) = q1_mean
%    mean_c(5) = q2_mean
%    mean_c(6) = Omega_mean
%
% output :
%    D_J2 = I + (-J2*Re^2)*(D_lp+D_sp1+D_sp2) = 6x6 transformation matrix D_J2
%    osc_c = osculating new set of elements (6x1)

%
gamma   = -J2*Re^2;
a       = meanElems(1);
argLat  = meanElems(2); 
inc     = meanElems(3);
q1      = meanElems(4);
q2      = meanElems(5);
RAAN    = meanElems(6);
%
si      = sin(inc);
ci      = cos(inc);
s2i     = sin(2*inc);
c2i     = cos(2*inc);
sth     = sin(argLat);
cth     = cos(argLat);
s2th    = sin(2*argLat);
c2th    = cos(2*argLat);
s3th    = sin(3*argLat);
c3th    = cos(3*argLat);
s4th    = sin(4*argLat);
c4th    = cos(4*argLat);
s5th    = sin(5*argLat);
c5th    = cos(5*argLat);
%
p   = a*(1 - (q1^2+q2^2));
R   = p/(1 + q1*cth + q2*sth);
Vr  = sqrt(mu/p)*(q1*sth - q2*cth);
Vt  = sqrt(mu/p)*(1 + q1*cth + q2*sth);
%
Ttheta  = 1/(1-5*ci^2);
eta     = sqrt(1 - (q1^2+q2^2));
eps1    = sqrt(q1^2 + q2^2);
eps2    = q1*cth + q2*sth;
eps3    = q1*sth - q2*cth;
[lambda] = theta2lam(a, argLat, q1, q2);
argLatLam = argLat - lambda;

lam_q1 = (q1*Vr)/(eta*Vt) + q2/(eta*(1+eta)) - eta*R*(a+R)*(q2+sin(argLat))/(p^2);
lam_q2 = (q2*Vr)/(eta*Vt) - q1/(eta*(1+eta)) + eta*R*(a+R)*(q1+cos(argLat))/(p^2);

% Long period part,  D_lp
lamLp = (si^2/(8*a^2*eta^2*(1+eta)))*(1-10*Ttheta*ci^2)*q1*q2 ...
   +(q1*q2/(16*a^2*eta^4))*(3-55*ci^2-280*Ttheta*ci^4-400*Ttheta^2*ci^6);
%
aLp        = 0;
argLatLp   = lamLp -(si^2/(16*a^2*eta^4))*(1-10*Ttheta*ci^2)*((3+2*eta^2/(1+eta))*q1*q2 ...
   +2*q1*sth+2*q2*cth+(1/2)*(q1^2+q2^2)*s2th);
incLp      = (s2i/(32*a^2*eta^4))*(1-10*Ttheta*ci^2)*(q1^2-q2^2);
q1Lp       = -(q1*si^2/(16*a^2*eta^2))*(1-10*Ttheta*ci^2) ...
   -(q1*q2^2/(16*a^2*eta^4))*(3-55*ci^2-280*Ttheta*ci^4-400*Ttheta^2*ci^6);
q2Lp       = (q2*si^2/(16*a^2*eta^2))*(1-10*Ttheta*ci^2) ...
   +(q1^2*q2/(16*a^2*eta^4))*(3-55*ci^2-280*Ttheta*ci^4-400*Ttheta^2*ci^6);
RAANLp     = (q1*q2*ci/(8*a^2*eta^4))*(11+80*Ttheta*ci^2+200*Ttheta^2*ci^4);
%
DLP11 = -(1/a)*aLp;
DLP12 = 0;
DLP13 = 0;
DLP14 = 0;
DLP15 = 0;
DLP16 = 0;

DLP21 = -(2/a)*argLatLp;
DLP22 = -(si^2/(16*a^2*eta^4))*(1-10*Ttheta*ci^2)*(2*(q1*cth-q2*sth)+eps1*c2th);
DLP23 = (s2i/(16*a^2*eta^4))*(5*q1*q2*(11+112*Ttheta*ci^2+520*Ttheta^2*ci^4+800*Ttheta^3*ci^6) ...
   -(2*q1*q2+(2+eps2)*(q1*sth+q2*cth))*((1-10*Ttheta*ci^2)+10*Ttheta*si^2*(1+5*Ttheta*ci^2)));
DLP24 = (1/(16*a^2*eta^6))*((eta^2+4*q1^2)*(q2*(3-55*ci^2-280*Ttheta*ci^4-400*Ttheta^2*ci^6) ...
   -si^2*(1-10*Ttheta*ci^2)*(3*q2+2*sth)) - 2*si^2*(1-10*Ttheta*ci^2)*(4*q2+sth*(1+eps1))*q1*cth);
DLP25 = (1/(16*a^2*eta^6))*((eta^2+4*q2^2)*(q1*(3-55*ci^2-280*Ttheta*ci^4-400*Ttheta^2*ci^6) ...
   -si^2*(1-10*Ttheta*ci^2)*(3*q1+2*cth)) - 2*si^2*(1-10*Ttheta*ci^2)*(4*q1+cth*(1+eps1))*q2*sth);
DLP26 = 0;

DLP31 = -(2/a)*incLp;
DLP32 = 0;
DLP33 = ((q1^2-q2^2)/(16*a^2*eta^4))*(c2i*(1-10*Ttheta*ci^2)+5*Ttheta*s2i^2*(1+5*Ttheta*ci^2));
DLP34 = (q1*s2i/(16*a^2*eta^6))*(1-10*Ttheta*ci^2)*(eta^2+2*(q1^2-q2^2));
DLP35 = -(q2*s2i/(16*a^2*eta^6))*(1-10*Ttheta*ci^2)*(eta^2-2*(q1^2-q2^2));
DLP36 = 0;

DLP41 = -(2/a)*q1Lp;
DLP42 = 0;
DLP43 = -(q1*s2i/(16*a^2*eta^4))*(eta^2*((1-10*Ttheta*ci^2)+10*Ttheta*si^2*(1+5*Ttheta*ci^2)) ...
   +5*q2^2*(11+112*Ttheta*ci^2+520*Ttheta^2*ci^4+800*Ttheta^3*ci^6));
DLP44 = -(1/(16*a^2*eta^6))*(eta^2*si^2*(1-10*Ttheta*ci^2)*(eta^2+2*q1^2) ...
   +q2^2*(eta^2+4*q1^2)*(3-55*ci^2-280*Ttheta*ci^4-400*Ttheta^2*ci^6));
DLP45 = -(q1*q2/(8*a^2*eta^6))*(eta^2*si^2*(1-10*Ttheta*ci^2) ...
   +(eta^2+2*q2^2)*(3-55*ci^2-280*Ttheta*ci^4-400*Ttheta^2*ci^6));
DLP46 = 0;

DLP51 = -(2/a)*q2Lp;
DLP52 = 0;
DLP53 = (q2*s2i/(16*a^2*eta^4))*(eta^2*(1-10*Ttheta*ci^2)+10*Ttheta*eta^2*si^2*(1+5*Ttheta*ci^2) ...
   +5*q1^2*(11+112*Ttheta*ci^2+520*Ttheta^2*ci^4+800*Ttheta^3*ci^6));
DLP54 = (q1*q2/(8*a^2*eta^6))*(eta^2*si^2*(1-10*Ttheta*ci^2) ...
   +(3-55*ci^2-280*Ttheta*ci^4-400*Ttheta^2*ci^6)*(eta^2+2*q1^2));
DLP55 = (1/(16*a^2*eta^6))*(eta^2*si^2*(1-10*Ttheta*ci^2)*(eta^2+2*q2^2) ...
   +q1^2*(3-55*ci^2-280*Ttheta*ci^4-400*Ttheta^2*ci^6)*(eta^2+4*q2^2));
DLP56 = 0;

DLP61 = -(2/a)*RAANLp;
DLP62 = 0;
DLP63 = -(q1*q2*si/(8*a^2*eta^4))*((11+80*Ttheta*ci^2+200*Ttheta^2*ci^4)+160*Ttheta*ci^2*(1+5*Ttheta*ci^2)^2);
DLP64 = (q2*ci/(8*a^2*eta^6)) *(eta^2+4*q1^2)*(11+80*Ttheta*ci^2+200*Ttheta^2*ci^4);
DLP65 = (q1*ci/(8*a^2*eta^6)) *(eta^2+4*q2^2)*(11+80*Ttheta*ci^2+200*Ttheta^2*ci^4);
DLP66 = 0;

DLP = [ DLP11  DLP12  DLP13  DLP14  DLP15  DLP16;
         DLP21  DLP22  DLP23  DLP24  DLP25  DLP26;
         DLP31  DLP32  DLP33  DLP34  DLP35  DLP36;
         DLP41  DLP42  DLP43  DLP44  DLP45  DLP46;
         DLP51  DLP52  DLP53  DLP54  DLP55  DLP56;
         DLP61  DLP62  DLP63  DLP64  DLP65  DLP66 ];
       
       
% First short period part,  D_sp1
lamSp1 = (eps3*(1-3*ci^2)/(4*a^2*eta^4*(1+eta)))*((1+eps2)^2+(1+eps2)+eta^2) ...
   +(3*(1-5*ci^2)/(4*a^2*eta^4))*(argLatLam+eps3);
%
aSp1 = ((1-3*ci^2)/(2*a*eta^6))*((1+eps2)^3-eta^3);
argLatSp1 = lamSp1 - (eps3*(1-3*ci^2)/(4*a^2*eta^4*(1+eta)))*((1+eps2)^2+eta*(1+eta));
IncSp1 = 0;
q1Sp1 = -(3*q2*(1-5*ci^2)/(4*a^2*eta^4))*(argLatLam+eps3) ...
   +((1-3*ci^2)/(4*a^2*eta^4*(1+eta)))*(((1+eps2)^2+eta^2)*(q1+(1+eta)*cth)+(1+eps2)*((1+eta)*cth+q1*(eta-eps2)));
q2Sp1 = (3*q1*(1-5*ci^2)/(4*a^2*eta^4))*(argLatLam+eps3) ...
   +((1-3*ci^2)/(4*a^2*eta^4*(1+eta)))*(((1+eps2)^2+eta^2)*(q2+(1+eta)*sth)+(1+eps2)*((1+eta)*sth+q2*(eta-eps2)));
RAANSp1 = (3*ci/(2*a^2*eta^4))*(argLatLam+eps3);
%
DSP111 = -(1/a)*aSp1;
DSP112 = -(3*eps3/(2*a*eta^6))*(1-3*ci^2)*(1+eps2)^2;
DSP113 = (3*s2i/(2*a*eta^6))*((1+eps2)^3-eta^3);
DSP114 = (3*(1-3*ci^2)/(2*a*eta^8))*(2*q1*(1+eps2)^3+eta^2*(1+eps2)^2*cth-eta^3*q1);
DSP115 = (3*(1-3*ci^2)/(2*a*eta^8))*(2*q2*(1+eps2)^3+eta^2*(1+eps2)^2*sth-eta^3*q2);
DSP116 = 0;

DSP121 = -(2/a)*argLatSp1;
DSP122 = ((1-3*ci^2)/(4*a^2*eta^4*(1+eta)))*(eps2*(1+eps2-eta)-eps3^2) ...
   +(3*(1-5*ci^2)/(4*a^2*eta^4*(1+eps2)^2))*((1+eps2)^3-eta^3);
DSP123 = (3*eps3*s2i/(4*a^2*eta^4*(1+eta)))*((1+eps2)+(5+4*eta)) +(15*s2i/(4*a^2*eta^4))*(argLatLam);
DSP124 = ((1-3*ci^2)/(4*a^2*eta^6*(1+eta)^2))*(eta^2*(eps1*sth+(1+eta)*(eps2*sth+eps3*cth)) ...
   +q1*eps3*(4*(eps1+eps2)+eta*(2+5*eps2))) +(3*(1-5*ci^2)/(4*a^2*eta^6))*(4*q1*(argLatLam+eps3)+eta^2*sth) ...
   -(3*(1-5*ci^2)/(4*a^2*eta^4))*(lam_q1);
DSP125 = -((1-3*ci^2)/(4*a^2*eta^6*(1+eta)^2))*(eta^2*(eps1*cth+(1+eta)*(eps2*cth-eps3*sth)) ...
   -q2*eps3*(4*(eps1+eps2)+eta*(2+5*eps2))) +(3*(1-5*ci^2)/(4*a^2*eta^6))*(4*q2*(argLatLam+eps3)-eta^2*cth) ...
   -(3*(1-5*ci^2)/(4*a^2*eta^4))*(lam_q2);
DSP126 = 0;

DSP131 = -(2/a)*IncSp1;
DSP132 = 0;
DSP133 = 0;
DSP134 = 0;
DSP135 = 0;
DSP136 = 0;

DSP141 = -(2/a)*q1Sp1;
DSP142 = -((1-3*ci^2)/(4*a^2*eta^4))*((1+eps2)*(2*sth+eps2*sth+2*eps3*cth)+eps3*(q1+cth)+eta^2*sth) ...
   -(3*q2*(1-5*ci^2)/(4*a^2*eta^4*(1+eps2)^2))*((1+eps2)^3-eta^3);
DSP143 = (3*q1*s2i/(4*a^2*eta^2*(1+eta))) ...
   +(3*s2i/(4*a^2*eta^4))*((1+eps2)*(q1+(2+eps2)*cth)-5*q2*eps3+eta^2*cth) -(15*q2*s2i/(4*a^2*eta^4))*(argLatLam);
DSP144 = ((1-3*ci^2)/(4*a^2*eta^2*(1+eta))) +((1-3*ci^2)*q1^2*(4+5*eta)/(4*a^2*eta^6*(1+eta)^2)) ...
   +((1-3*ci^2)/(8*a^2*eta^6))*(eta^2*(5+2*(5*q1*cth+2*q2*sth)+(3+2*eps2)*c2th) ...
   +2*q1*(4*(1+eps2)*(2+eps2)*cth+(3*eta+4*eps2)*q1)) ...
   -(3*q2*(1-5*ci^2)/(4*a^2*eta^6))*(4*q1*eps3+eta^2*sth) -(3*q1*q2*(1-5*ci^2)/(a^2*eta^6))*(argLatLam) ...
   +(3*q2*(1-5*ci^2)/(4*a^2*eta^4))*(lam_q1);
DSP145 = ((1-3*ci^2)/(8*a^2*eta^6))*(eta^2*(2*(q1*sth+2*q2*cth)+(3+2*eps2)*s2th) ...
   +2*q2*(4*(1+eps2)*(2+eps2)*cth+(3*eta+4*eps2)*q1)) +((1-3*ci^2)*q1*q2*(4+5*eta)/(4*a^2*eta^6*(1+eta)^2)) ...
   -(3*(1-5*ci^2)/(4*a^2*eta^6))*(eps3*(eta^2+4*q2^2)-eta^2*q2*cth) ...
   -(3*(1-5*ci^2)/(4*a^2*eta^6))*(argLatLam*(eta^2+4*q2^2)) +(3*q2*(1-5*ci^2)/(4*a^2*eta^4))*(lam_q2);
DSP146 = 0;

DSP151 = -(2/a)*q2Sp1;
DSP152 = ((1-3*ci^2)/(4*a^2*eta^4))*((1+eps2)*(2*cth+eps2*cth-2*eps3*sth)-eps3*(q2+sth)+eta^2*cth) ...
   +(3*q1*(1-5*ci^2)/(4*a^2*eta^4*(1+eps2)^2))*((1+eps2)^3-eta^3);
DSP153 = (3*q2*s2i/(4*a^2*eta^2*(1+eta))) ...
   +(3*s2i/(4*a^2*eta^4))*((1+eps2)*(q2+(2+eps2)*sth)+5*q1*eps3+eta^2*sth) -(15*q1*s2i/(4*a^2*eta^4))*(argLatLam);
DSP154 = ((1-3*ci^2)/(8*a^2*eta^6))*(eta^2*(2*(2*q1*sth+q2*cth)+(3+2*eps2)*s2th) ...
   +2*q1*(4*(1+eps2)*(2+eps2)*sth+(3*eta+4*eps2)*q2)) +((1-3*ci^2)*q1*q2*(4+5*eta)/(4*a^2*eta^6*(1+eta)^2)) ...
   +(3*(1-5*ci^2)/(4*a^2*eta^6))*(eps3*(eta^2+4*q1^2)+eta^2*q1*sth) ...
   +(3*(1-5*ci^2)/(4*a^2*eta^6))*(argLatLam*(eta^2+4*q1^2)) -(3*q1*(1-5*ci^2)/(4*a^2*eta^4))*(lam_q1);
DSP155 = ((1-3*ci^2)/(4*a^2*eta^2*(1+eta))) +((1-3*ci^2)*q2^2*(4+5*eta)/(4*a^2*eta^6*(1+eta)^2)) ...
   +((1-3*ci^2)/(8*a^2*eta^6))*(eta^2*(5+2*(2*q1*cth+5*q2*sth)-(3+2*eps2)*c2th) ...
   +2*q2*(4*(1+eps2)*(2+eps2)*sth+(3*eta+4*eps2)*q2)) ...
   +(3*q1*(1-5*ci^2)/(4*a^2*eta^6))*(4*q2*eps3-eta^2*cth) +(3*q1*q2*(1-5*ci^2)/(a^2*eta^6))*(argLatLam) ...
   -(3*q1*(1-5*ci^2)/(4*a^2*eta^4))*(lam_q2);
DSP156 = 0;

DSP161 = -(2/a)*RAANSp1;
DSP162 = (3*ci/(2*a^2*eta^4*(1+eps2)^2))*((1+eps2)^3-eta^3);
DSP163 = -(3*eps3*si/(2*a^2*eta^4)) -(3*si/(2*a^2*eta^4))*(argLatLam);
DSP164 = (3*ci/(2*a^2*eta^6))*(4*q1*eps3+eta^2*sth) +(6*q1*ci/(a^2*eta^6))*(argLatLam) ...
   -(3*ci/(2*a^2*eta^4))*(lam_q1);
DSP165 = (3*ci/(2*a^2*eta^6))*(4*q2*eps3-eta^2*cth) +(6*q2*ci/(a^2*eta^6))*(argLatLam) ...
   -(3*ci/(2*a^2*eta^4))*(lam_q2);
DSP166 = 0;

D_sp1 = [ DSP111  DSP112  DSP113  DSP114  DSP115  DSP116;
          DSP121  DSP122  DSP123  DSP124  DSP125  DSP126;
          DSP131  DSP132  DSP133  DSP134  DSP135  DSP136;
          DSP141  DSP142  DSP143  DSP144  DSP145  DSP146;
          DSP151  DSP152  DSP153  DSP154  DSP155  DSP156;
          DSP161  DSP162  DSP163  DSP164  DSP165  DSP166 ];
       
       
% Second short period part,  D_sp2
lamSp2 = -(3*eps3*si^2*c2th/(4*a^2*eta^4*(1+eta)))*(1+eps2)*(2+eps2) ...
   -(si^2/(8*a^2*eta^2*(1+eta)))*(3*(q1*sth+q2*cth)+(q1*s3th-q2*c3th)) ...
   -((3-5*ci^2)/(8*a^2*eta^4))*(3*(q1*sth+q2*cth)+3*s2th+(q1*s3th-q2*c3th));
%
aSp2 = -(3*si^2/(2*a*eta^6))*(1+eps2)^3*c2th;
argLatSp2 = lamSp2 -(si^2/(32*a^2*eta^4*(1+eta)))*(36*q1*q2-4*(3*eta^2+5*eta-1)*(q1*sth+q2*cth) ...
   +12*eps2*q1*q2-32*(1+eta)*s2th-(eta^2+12*eta+39)*(q1*s3th-q2*c3th) ...
   +36*q1*q2*c4th-18*(q1^2-q2^2)*s4th+3*q2*(3*q1^2-q2^2)*c5th-3*q1*(q1^2-3*q2^2)*s5th);
incSp2 = -(s2i/(8*a^2*eta^4))*(3*(q1*cth-q2*sth)+3*c2th+(q1*c3th+q2*s3th));
q1Sp2 = (q2*(3-5*ci^2)/(8*a^2*eta^4))*(3*(q1*sth+q2*cth)+3*s2th+(q1*s3th-q2*c3th)) ...
   +(si^2/(8*a^2*eta^4))*(3*(eta^2-q1^2)*cth+3*q1*q2*sth-(eta^2+3*q1^2)*c3th-3*q1*q2*s3th) ...
   -(3*si^2*c2th/(16*a^2*eta^4))*(10*q1+(8+3*q1^2+q2^2)*cth+2*q1*q2*sth ...
   +6*(q1*c2th+q2*s2th)+(q1^2-q2^2)*c3th+2*q1*q2*s3th);
q2Sp2 = -(q1*(3-5*ci^2)/(8*a^2*eta^4))*(3*(q1*sth+q2*cth)+3*s2th+(q1*s3th-q2*c3th)) ...
   -(si^2/(8*a^2*eta^4))*(3*(eta^2-q2^2)*sth+3*q1*q2*cth+(eta^2+3*q2^2)*s3th+3*q1*q2*c3th) ...
   -(3*si^2*c2th/(16*a^2*eta^4))*(10*q2+(8+q1^2+3*q2^2)*sth+2*q1*q2*cth ...
   +6*(q1*s2th-q2*c2th)+(q1^2-q2^2)*s3th-2*q1*q2*c3th);
RAANSp2 = -(ci/(4*a^2*eta^4))*(3*(q1*sth+q2*cth)+3*s2th+(q1*s3th-q2*c3th));
%
DSP211 = -(1/a)*aSp2;
DSP212 = (3*si^2/(2*a*eta^6))*(1+eps2)^2*(3*eps3*c2th+2*(1+eps2)*s2th);
DSP213 = -(3*s2i*c2th/(2*a*eta^6))*(1+eps2)^3;
DSP214 = -(9*si^2*c2th/(2*a*eta^8))*(1+eps2)^2*(2*q1*(1+eps2)+eta^2*cth);
DSP215 = -(9*si^2*c2th/(2*a*eta^8))*(1+eps2)^2*(2*q2*(1+eps2)+eta^2*sth);
DSP216 = 0;

DSP221 = -(2/a)*argLatSp2;
DSP222 = -(1/(8*a^2*eta^4))*(3*(3-5*ci^2)*((q1*cth-q2*sth)+2*c2th+(q1*c3th+q2*s3th)) ...
   -si^2*(5*(q1*cth-q2*sth)+16*c2th+9*(q1*c3th+q2*s3th)));
DSP223 = -(s2i/(8*a^2*eta^4))*(10*(q1*sth+q2*cth)+7*s2th+2*(q1*s3th-q2*c3th));
DSP224 = -((3-5*ci^2)/(8*a^2*eta^6))*(4*q1*(3*s2th+q2*(3*cth-c3th))+(eta^2+4*q1^2)*(3*sth+s3th)) ...
   -(si^2*(3*sth+s3th)/(8*a^2*eta^2*(1+eta))) -(si^2/(32*a^2*eta^4*(1+eta)))*(36*q2-4*(2+3*eta)*sth ...
   -(eta*(12+eta)+39)*s3th+9*eps1*s5th+12*q2*(2*q1*cth+q2*sth)+9*q1*(q1*s3th-q2*c3th)+18*(3*q1*s4th+2*q2*c4th) ...
   -3*q1*(q1*s5th-11*q2*c5th)+24*((1+eps2)*(2+eps2)*sth+eps3*(3+2*eps2)*cth)*c2th) ...
   -(3*si^2/(32*a^2*eta^4*(1+eta)^2))*(4*sth-6*q1*s4th-q1*(q1*s5th+q2*c5th)) ...
   +(q1*si^2/(8*a^2*eta^6*(1+eta)))*(20*(1+eta)*(q1*sth+q2*cth)+32*(1+eta)*s2th+3*(4+3*eta)*(q1*s3th-q2*c3th)) ...
   -(q1*si^2*(4+5*eta)/(32*a^2*eta^6*(1+eta)^2))*(24*(q1*sth+q2*cth)+24*eps3*(1+eps2)*(2+eps2)*c2th ...
   -(27+3*eta)*(q1*s3th-q2*c3th)-18*s4th-3*(q1*s5th+q2*c5th) ...
   +12*q2*((3+eps2)*q1+3*(q1*c4th+q2*s4th)+q1*(q1*c5th+q2*s5th)));
DSP225 = -((3-5*ci^2)/(8*a^2*eta^6))*(4*q2*(3*s2th+q1*(3*sth+s3th))+(eta^2+4*q2^2)*(3*cth-c3th)) ...
   -(si^2*(3*cth-c3th)/(8*a^2*eta^2*(1+eta))) -(si^2/(32*a^2*eta^4*(1+eta)))*(36*q1-4*(2+3*eta)*cth ...
   +(eta*(12+eta)+39)*c3th+9*eps1*c5th+12*q1*(q1*cth+2*q2*sth)+9*q2*(q1*s3th-q2*c3th)+18*(2*q1*c4th+7*q2*s4th) ...
   +3*q2*(11*q1*s5th-q2*c5th)+24*(eps3*(3+2*eps2)*sth-(1+eps2)*(2+eps2)*cth)*c2th) ...
   -(3*si^2/(32*a^2*eta^4*(1+eta)^2))*(4*cth-6*q2*s4th-q2*(q1*s5th+q2*c5th)) ...
   +(q2*si^2/(8*a^2*eta^6*(1+eta)))*(20*(1+eta)*(q1*sth+q2*cth)+32*(1+eta)*s2th+3*(4+3*eta)*(q1*s3th-q2*c3th)) ...
   -(q2*si^2*(4+5*eta)/(32*a^2*eta^6*(1+eta)^2))*(24*(q1*sth+q2*cth)+24*eps3*(1+eps2)*(2+eps2)*c2th ...
   -(27+3*eta)*(q1*s3th-q2*c3th)-18*s4th-3*(q1*s5th+q2*c5th) ...
   +12*q2*((3+eps2)*q1+3*(q1*c4th+q2*s4th)+q1*(q1*c5th+q2*s5th)));
DSP226 = 0;

DSP231 = -(2/a)*incSp2;
DSP232 = (3*s2i/(8*a^2*eta^4))*((q1*sth+q2*cth)+2*s2th+(q1*s3th-q2*c3th));
DSP233 = -(c2i/(4*a^2*eta^4))*(3*(q1*cth-q2*sth)+3*c2th+(q1*c3th+q2*s3th));
DSP234 = -(s2i/(8*a^2*eta^6))*(4*q1*(3*c2th-q2*(3*sth-s3th))+(eta^2+4*q1^2)*(3*cth+c3th));
DSP235 = -(s2i/(8*a^2*eta^6))*(4*q2*(3*c2th+q1*(3*cth+c3th))-(eta^2+4*q2^2)*(3*sth-s3th));
DSP236 = 0;

DSP241 = -(2/a)*q1Sp2;
DSP242 = (3*q2*(3-5*ci^2)/(8*a^2*eta^4))*((q1*cth-q2*sth)+2*c2th+(q1*c3th+q2*s3th)) ...
   +(3*si^2/(16*a^2*eta^4))*((2*eps2*q2-9*q2*(q1*c3th+q2*s3th)+12*(q1*s4th-q2*c4th)-5*q2*(q1*c5th+q2*s5th)) ...
   +(1/2)*(4*(1+3*q1^2)*sth+40*q1*s2th+(28+17*eps1)*s3th+5*eps1*s5th));
DSP243 = -(s2i/(16*a^2*eta^4))*((36*q1*(q1*cth-q2*sth)+30*(q1*c2th-q2*s2th)-q2*(q1*s3th-q2*c3th) ...
   +9*(q1*c4th+q2*s4th)+3*q2*(q1*s5th-q2*c5th)) ...
   +(1/2)*(6*q1*(3+2*q1*cth)+12*(1-4*eps1)*cth+(28+17*eps1)*c3th+3*eps1*c5th));
DSP244 = (q2*(3-5*ci^2)/(8*a^2*eta^6))*(4*q1*(3*s2th+q2*(3*cth-c3th))+(eta^2+4*q1^2)*(3*sth+s3th)) ...
   -(si^2/(8*a^2*eta^4))*((8*q1*c3th-3*q2*(sth-s3th))+3*(5+eps2+3*c2th+3*(q1*c3th+q2*s3th))*c2th) ...
   -(3*q1*si^2/(4*a^2*eta^6))*(2*q1*((q1*cth-q2*sth)+(q1*c3th+q2*s3th)) ...
   +(9*cth-c3th+2*q1*(5+eps2)+6*(q1*c2th+q2*s2th)+2*q1*(q1*c3th+q2*s3th))*c2th);
DSP245 = ((3-5*ci^2)/(8*a^2*eta^6))*((eta^2+4*q2^2)*(3*s2th+q1*(3*sth+s3th)) ...
   +2*(eta^2+2*q2^2)*q2*(3*cth-c3th)) +(si^2/(16*a^2*eta^4))*(6*(q1*sth+2*q2*cth) ...
   -(9*q1*s3th+q2*c3th)-9*s4th-3*(q1*s5th+q2*c5th)) ...
   -(3*q2*si^2/(8*a^2*eta^6))*(2*q1*(3+2*(2*q1*cth-q2*sth)+10*c2th+3*(q1*c3th+q2*s3th)+(q1*c5th+q2*s5th)) ...
   +(8*cth+9*c3th+6*(q1*c4th+q2*s4th)-c5th));
DSP246 = 0;

DSP251 = -(2/a)*q2Sp2;
DSP252 = -(3*q1*(3-5*ci^2)/(8*a^2*eta^4))*((q1*cth-q2*sth)+2*c2th+(q1*c3th+q2*s3th)) ...
   +(3*si^2/(16*a^2*eta^4))*((2*eps2*q1+9*q1*(q1*c3th+q2*s3th)-12*(q1*c4th+q2*s4th)-5*q1*(q1*c5th+q2*s5th)) ...
   +(1/2)*(4*(1+3*q2^2)*cth+40*q2*s2th-(28+17*eps1)*c3th+5*eps1*c5th));
DSP253 = -(s2i/(16*a^2*eta^4))*((36*q1*(q1*sth+q2*cth)+30*(q1*s2th+q2*c2th)+q1*(q1*s3th-q2*c3th) ...
   +9*(q1*s4th-q2*c4th)+3*q1*(q1*s5th-q2*c5th)) ...
   -(1/2)*(6*q2*(3+2*q2*sth)+12*(1+2*eps1)*sth-(28+17*eps1)*s3th+3*eps1*s5th));
DSP254 = -((3-5*ci^2)/(8*a^2*eta^6))*((eta^2+4*q1^2)*(3*s2th+q2*(3*cth-c3th)) ...
   +2*(eta^2+2*q1^2)*q1*(3*sth+s3th)) -(si^2/(16*a^2*eta^4))*(6*(2*q1*sth+q2*cth) ...
   +(q1*s3th+9*q2*c3th)+9*s4th-3*(q1*s5th+q2*c5th)) ...
   +(3*q1*si^2/(8*a^2*eta^6))*(2*q2*(3-2*(q1*cth-2*q2*sth)-10*c2th-3*(q1*c3th+q2*s3th)+(q1*c5th+q2*s5th)) ...
   +(8*sth-9*s3th-6*(q1*s4th-q2*c4th)-s5th));
DSP255 = -(q1*(3-5*ci^2)/(8*a^2*eta^6))*((eta^2+4*q2^2)*(3*cth-c3th)+4*q2*(3*s2th+q1*(3*sth+s3th))) ...
   -(si^2/(8*a^2*eta^4))*(8*q2*s3th+3*q1*(cth+c3th)+3*(5+eps2-3*c2th-(q1*c3th-q2*s3th))*c2th) ...
   -(3*si^2*q2*c2th/(4*a^2*eta^6))*(9*sth-s3th+2*q2*(5+eps2)+6*(q1*s2th-q2*c2th)+2*q1*(q1*s3th-q2*c3th));
DSP256 = 0;

DSP261 = -(2/a)*RAANSp2;
DSP262 = -(3*ci/(4*a^2*eta^4))*((q1*cth-q2*sth)+2*c2th+(q1*c3th+q2*s3th));
DSP263 = (si/(4*a^2*eta^4))*(3*(q1*sth+q2*cth)+3*s2th+(q1*s3th-q2*c3th));
DSP264 = -(ci/(4*a^2*eta^6))*(4*q1*(3*s2th+q2*(3*cth-c3th))+(eta^2+4*q1^2)*(3*sth+s3th));
DSP265 = -(ci/(4*a^2*eta^6))*(4*q2*(3*s2th+q1*(3*sth+s3th))+(eta^2+4*q2^2)*(3*cth-c3th));
DSP266 = 0;

D_sp2 = [ DSP211  DSP212  DSP213  DSP214  DSP215  DSP216;
          DSP221  DSP222  DSP223  DSP224  DSP225  DSP226;
          DSP231  DSP232  DSP233  DSP234  DSP235  DSP236;
          DSP241  DSP242  DSP243  DSP244  DSP245  DSP246;
          DSP251  DSP252  DSP253  DSP254  DSP255  DSP256;
          DSP261  DSP262  DSP263  DSP264  DSP265  DSP266 ];

%
var_lsp = gamma*[ aLp aSp1 aSp2; argLatLp argLatSp1 argLatSp2; incLp IncSp1 incSp2;
                   q1Lp q1Sp1 q1Sp2; q2Lp q2Sp1 q2Sp2; RAANLp RAANSp1 RAANSp2 ];
% Evaluating Osculaing Elements
lamOsc = lambda + gamma*(lamLp + lamSp1 + lamSp2);
%
aOsc       = a + gamma*(aLp + aSp1 + aSp2);
argLatOsc  = argLat + gamma*(argLatLp + argLatSp1 + argLatSp2);
iOsc       = inc + gamma*(incLp + IncSp1 + incSp2);
q1Osc      = q1 + gamma*(q1Lp + q1Sp1 + q1Sp2);
q2Osc      = q2 + gamma*(q2Lp + q2Sp1 + q2Sp2);
OmegaOsc   = RAAN + gamma*(RAANLp + RAANSp1 + RAANSp2);

% Transformation Matrix D_J2
DJ2 = eye(6) + gamma*(DLP + D_sp1 + D_sp2);

% Osculating Elements from Mean elements
osc_c = [ aOsc;  argLatOsc;  iOsc;  q1Osc;  q2Osc;  OmegaOsc ];
end
clc;
clear all;
global G G0 E E0 v v0 mu mu0 kPsi kPsi0 kr kr0 ar ar0 k1 k0 k00 a a0 sigma0 rA rB rC rD sigmaA pwA pwE m0 m1...
    C1RA C2RA ZRA RE PA PR...
    PAA ZRB C1RB C2RB REB PRB PBB PRA C1RBA C2RBA...
    REC PCC C1RCBA C2RCBA C2RCB C1RCB PRC C1RC C2RC 

%%%%%%%%%%size of the model%%%%%%%%%%
rA = 4.5;
% radius of the tunnel
rB = 7.5;
% radius of bolts-grouting
rC = 10.5;
%rC=rg,radius of grouting
rD = 60.0;
% model radius considering seepage
rE = 500.0;
% model radius 

%%%%%%%%%%seepage parameter%%%%%%%%%%
m1 = 1.2;
% empirical coefficient for grout region
m0 = 1.0001;
% empirical coefficient for original surrounding rock
k1 = 1.0*10^(-9);
%permeability coefficient for grout region
k0 = 1.0*10^(-8);
%permeability coefficient for original surrounding rock

%%%%%%%%%%boundary conditions%%%%%%%%%%
pwA = 0.;
%pore pressure at r=rA
pwE = 2.5e6;
%pore pressure at before excavation
pwD = pwE;
%pore pressure for region C
gammaw = 10000.0;
%weight of water
sigmaE = 5.0e6;
%total radial stress before excavation
sigma0=sigmaE-pwE;
%total effective radial stress before excavation
Pvir=0.0.*sigma0;
%virtual support pressure
PWI=0.0;
%Water pressure when the tunnel is full of water
sigmaA=PWI+pwA+Pvir;
%sigmaA=P
%%%%%%%%%%rockbolts reinforcement parameters%%%%%%%%%%
db = 0.032;
%diameter of rockbolts
nb = 25;
%rockbolts number
S = ((2.*3.14.*rA./nb)+(2.*3.14.*rB./nb))./2;
%ring spacing of the rockbolts
DS=0.5;
%spacing of rockbolts along the tunnel excavation direction
Y = 3.14.*(0.5*db)^2;
%cross-sectional area of the rockbolt
Eb = 210e9;
%modulus of elasticity of the rockbolts
arfa = (nb.*Y)./(2.*3.14.*rA.*DS);
%rockbolts reinforcement density
sigmabf= 1e6;
%Residual axial stresses on the rockbolt

%%%%%%%%%%elastic parameters of the surrounding rock%%%%%%%%%%
E = 1500e6;
%modulus of elasticity of the grouted surrounding rock
v = 0.25;
%Poisson's ratio of the grouted surrounding rock
G = E./(2.*(1 + v));
%shear modulus of the grouted surrounding rock
mu = v./(1 - 2.*v);
%elastic constant of the grouted surrounding rock
E0 = 1200e6;
%modulus of elasticity of the original surrounding rock
v0 = 0.31;
%Poisson's ratio of the original surrounding rock
G0 = E0./(2.*(1 + v0));
%shear modulus of the original surrounding rock
mu0 = v0./(1 - 2.*v0);
%elastic constant of the original surrounding rock

beita = (arfa.*Eb)./(2.*G);
%ground-reinforcement stiffness coefficient (rockbolts reinforcement parameters)

%%%%%%%%%%plastic parameters of the surrounding rock%%%%%%%%%%
%%%before brittle-plastic destruction%%%%%%
c = 0.30e6;
%cohesion for surrounding rock reinforced by grout
phi = 3.14.*28./180;
%friction angle for surrounding rock reinforced by grout
k = tan(3.14./4 + phi./2).^2;
%Mohr每Coulomb constant for surrounding rock reinforced by grout
a = c./tan(phi);
%Mohr每Coulomb constant for surrounding rock reinforced by grout
c0 = 0.20e6;
%cohesion for surrounding rock without reinforcement
phi0 = 3.14.*20./180;
%friction angle for surrounding rock without reinforcement
k00 = tan(3.14./4 + phi0./2).^2;
%Mohr每Coulomb constant for surrounding rock without reinforcement
a0 = c0./tan(phi0);
%Mohr每Coulomb constant for surrounding rock without reinforcement

%%%after brittle-plastic destruction%%%%%%
cr = 0.10e6;
%cohesion for surrounding rock reinforced by grout in region of brittle-plastic destruction
phir = 3.14.*15./180;
%friction angle for surrounding rock reinforced by grout in region of brittle-plastic destruction
kr = tan(3.14./4 + phir./2).^2;
%Mohr每Coulomb constant for surrounding rock reinforced by grout in region of brittle-plastic destruction
ar = cr./tan(phir);
%Mohr每Coulomb constant for surrounding rock reinforced by grout in region of brittle-plastic destruction
cr0 = 0.05e6;
%cohesion for surrounding rock without reinforcement in region of brittle-plastic destruction
phir0 = 3.14.*13./180;
%friction angle for surrounding rock without reinforcement in region of brittle-plastic destruction
kr0 = tan(3.14./4 + phir0./2).^2;
%Mohr每Coulomb constant for surrounding rock without reinforcement in region of brittle-plastic destruction
ar0 = cr0./tan(phir0);
%Mohr每Coulomb constant for surrounding rock without reinforcement in region of brittle-plastic destruction

psi = 3.14*30/180;
%dilation angle for regions reinforced by grout
kPsi = tan(3.14/4 + psi/2).^2;
%dilation factor for regions reinforced by grout
psi0 = 3.14*35/180;
%dilation angle for unreinforced regions 
kPsi0 = tan(3.14/4 + psi/2).^2;
%dilation factor for unreinforced regions 

pwC=(k1.*m1.*pwA.*rD.*exp(m1.*log(rA./rC) + m0.*log(rC./rD)) ...
    - k1.*pwA.*rD.*exp(m1.*log(rA./rC) + m0.*log(rC./rD)) ...
    - exp(m1.*log(rA./rC)).*k1.*m1.*pwA.*rC + exp(m1.*log(rA./rC)).*k0.*m0.*pwD.*rC ...
    + exp(m1.*log(rA./rC)).*k1.*pwA.*rC - exp(m1.*log(rA./rC)).*k0.*pwD.*rC ...
    - k0.*m0.*pwD.*rA + k0.*pwD.*rA)./(exp(m1.*log(rA./rC)).*(k0.*m0.*rC ...
    - k0.*rC + exp(m0.*log(rC./rD)).*k1.*m1.*rD - exp(-m1.*log(rA./rC)).*k0.*m0.*rA ...
    - exp(m0.*log(rC./rD)).*k1.*rD + exp(-m1.*log(rA./rC)).*k0.*rA - k1.*m1.*rC + k1.*rC));
%pore pressure at grouting of reinforced boundary
C1W = -exp((m1.*log(rA./rC) + log((m1 - 1).*k1.*(pwA - pwC)./(gammaw.*(rC.*exp(m1.*log(rA./rC)) - rA))))./m1).*rC;
C2W = (rC.*pwA.*exp(m1.*log(rA./rC)) -pwC.*rA)./(rC.*exp(m1.*log(rA./rC)) - rA);
C1W0 = -exp((m0.*log(rC./rD) + log((m0 - 1).*k0.*(pwC - pwD)./(gammaw.*(rD.*exp(m0.*log(rC./rD)) - rC))))./m0).*rD;
C2W0 = (rD.*pwC.*exp(m0.*log(rC./rD)) - pwD.*rC)./(rD.*exp(m0.*log(rC./rD)) - rC);
%integration constants for pore pressure general solution
NA1 = sqrt(2.0).*sqrt(DS).*sqrt(G).*sqrt(S).*sqrt(1 + mu)./sqrt(2.0.*G.*mu.*DS.*S ...
    + 2.*G.*DS.*S + Eb.*Y);
NA2 = -gammaw.*S.*DS./(2.*k1.*(G.*S.*(m1 - 1).*(m1 - 3).*(1 + mu).*DS ...
    + Y.*Eb.*(m1 - 2).^2./2));
NB = -gammaw./(2.*G.*k1.*(m1 - 1).*(m1 - 3).*(1 + mu));
NC = -gammaw./(2.*G0.*k0.*(m0 - 1).*(m0 - 3).*(1 + mu0));
%constant for displacement general solution
A1 = 2.*G.*(1 + mu).*rA.^NA1.*NA1./rA + 2.*G.*mu.*rA.^NA1./rA ...
    + rA.^NA1.*NA1.*Eb.*Y./(rA.*S.*DS);
A2 = -2.*G.*(1 + mu).*rA.^(-NA1).*NA1./rA + 2.*G.*mu.*rA.^(-NA1)./rA ...
    - rA.^(-NA1).*NA1.*Eb.*Y./(rA.*S.*DS);
A3 =2.*G.*(1 + mu).*(-(-C1W./rA).^m1.*m1.*rA.*NA2 + 2.*(-C1W./rA).^m1.*rA.*NA2) ...
    + 2.*G.*mu.*(-C1W./rA).^m1.*rA.*NA2 + sigmaE - pwE + (-(-C1W/rA).^m1.*m1.*rA.*NA2 ...
    + 2.*(-C1W./rA).^m1.*rA.*NA2).*Eb.*Y./(S.*DS) - rA.*(-C1W./rA).^m1.*gammaw./((m1 - 1).*k1) + C2W;
A4 = 2.*G.*(1 + mu).*rB.^NA1.*NA1./rB + 2.*G.*mu.*rB.^NA1./rB ...
    + rB.^NA1.*NA1.*Eb.*Y./(rB.*S.*DS);
A5 = -2.*G.*(1 + mu).*rB.^(-NA1).*NA1./rB + 2.*G.*mu.*rB.^(-NA1)./rB ...
    - rB.^(-NA1).*NA1.*Eb.*Y./(rB.*S.*DS);
A6 =2.*G.*(1 + mu).*(-(-C1W./rB).^m1.*m1.*rB.*NA2 + 2.*(-C1W./rB).^m1.*rB.*NA2) ...
    + 2.*G.*mu.*(-C1W./rB).^m1.*rB.*NA2 + sigmaE - pwE + (-(-C1W./rB)^m1.*m1.*rB.*NA2 ...
    + 2.*(-C1W./rB).^m1.*rB.*NA2).*Eb.*Y./(S.*DS) - rB.*(-C1W./rB).^m1.*gammaw./((m1 - 1).*k1) + C2W;
%constant for displacement general solution in region A
B1 = -2.*G.*(1 + mu)./rB.^2 + 2.*G.*mu./rB.^2;
B2 = 2.*G.*(1 + mu) + 2.*G.*mu;
B3 = 2.*G.*(1 + mu).*(-rB.*(-C1W./rB).^m1.*m1.*NB + 2.*rB.*(-C1W./rB).^m1.*NB) ...
    + 2.*G.*mu.*rB.*(-C1W./rB).^m1.*NB + sigmaE - pwE - rB.*(-C1W./rB).^m1.*gammaw./((m1 - 1).*k1) + C2W;
B4 = -2.*G.*(1 + mu)./rC.^2 + 2.*G.*mu./rC.^2;
B5 = 2.*G.*(1 + mu) + 2.*G.*mu;
B6 = 2.*G.*(1 + mu).*(2.*rC.*(-C1W./rC).^m1.*NB - rC.*(-C1W./rC).^m1.*m1.*NB) + 2.*G.*mu.*rC.*(-C1W./rC).^m1.*NB + sigmaE ...
    - pwE - rC.*(-C1W./rC).^m1.*gammaw./((m1 - 1).*k1) + C2W;
%constant for displacement general solution in region B
C1 = -2.*G0.*(1 + mu0)./rC.^2 + 2.*G0.*mu0./rC.^2;
C2 = 2.*G0.*(1 + mu0) + 2.*G0.*mu0;
C3 =2.*G0.*(1 + mu0).*(-rC.*(-C1W0./rC).^m0.*m0.*NC + 2.*rC.*(-C1W0./rC).^m0.*NC) ...
    + 2.*G0.*mu0.*rC.*(-C1W0./rC).^m0.*NC + sigmaE - pwE - rC.*(-C1W0./rC).^m0.*gammaw./((m0 - 1).*k0) + C2W0;
C4 = -2.*G0.*(1 + mu0)./rD.^2 + 2.*G0.*mu0./rD.^2;
C5 = 2.*G0.*(1 + mu0) + 2.*G0.*mu0;
C6 = 2.*G0.*(1 + mu0).*(2.*rD.*(-C1W0./rD).^m0.*NC - rD.*(-C1W0./rD).^m0.*m0.*NC)...
    + 2.*G0.*mu0.*rD.*(-C1W0./rD).^m0.*NC + sigmaE - pwE - rD.*(-C1W0./rD).^m0.*gammaw./((m0 - 1).*k0) + C2W0;
%constant for displacement general solution in region C
D1 = 2.*G0.*(1 + mu0) + 2.*G0.*mu0;
D2 = -2.*G0.*(1 + mu0)./rD.^2 + 2.*G0.*mu0./rD.^2;
D3 = sigmaE;
D4 = 2.*G0.*(1 + mu0) + 2.*G0.*mu0;
D5 = -2.*G0.*(1 + mu0)./rE.^2 + 2.*G0.*mu0./rE.^2;
D6 = sigmaE;
%constant for displacement general solution in region D
E1 = -A4./(A1.*A5 - A2.*A4);
E2 = A1./(A1.*A5 - A2.*A4);
E3 = -(A1.*A6 - A3.*A4)./(A1.*A5 - A2.*A4);
E4 = A5./(A1.*A5 - A2.*A4);
E5 = -A2./(A1.*A5 - A2.*A4);
E6 = (A2.*A6 - A3.*A5)./(A1.*A5 - A2.*A4);
%constant for displacement general solution in region A
F1 = B1./(B1.*B5 - B2.*B4);
F2 = -B4./(B1.*B5 - B2.*B4);
F3 = -(B1.*B6 - B3.*B4)./(B1.*B5 - B2.*B4);
F4 = -B2./(B1.*B5 - B2.*B4);
F5 = B5./(B1.*B5 - B2.*B4);
F6 = (B2.*B6 - B3.*B5)./(B1.*B5 - B2.*B4);
%constant for displacement general solution in region B
G1 = C1./(C1.*C5 - C2.*C4);
G2 = -C4./(C1.*C5 - C2.*C4);
G3 = -(C1.*C6 - C3.*C4)./(C1.*C5 - C2.*C4);
G4 = -C2./(C1.*C5 - C2.*C4);
G5 = C5./(C1.*C5 - C2.*C4);
G6 = (C2.*C6 - C3.*C5)./(C1.*C5 - C2.*C4);
%constant for displacement general solution in region C
H1 = D1./(D1.*D5 - D2.*D4);
H2 = -D4./(D1.*D5 - D2.*D4);
H3 = -(D1.*D6 - D3.*D4)./(D1.*D5 - D2.*D4);
H4 = -D2./(D1.*D5 - D2.*D4);
H5 = D5./(D1.*D5 - D2.*D4);
H6 = (D2.*D6 - D3.*D5)./(D1.*D5 - D2.*D4);
%constant for displacement general solution in region D
M1 = rB.^(-NA1).*E1 + rB.^NA1.*E4;
M2 = rB.^(-NA1).*E2 + rB.^NA1.*E5;
M3= E6.*rB.^NA1 + E3.*rB.^(-NA1) + (-C1W./rB).^m1.*rB.^2.*NA2;
M4 = rB.*F1 + F4./rB;
M5 = rB.*F2 + F5./rB;
M6=F3.*rB + F6./rB + rB.^2.*(-C1W/rB).^m1.*NB;
M7 = rC.*F1 + F4./rC;
M8 = rC.*F2 + F5./rC;
M9=F3.*rC + F6./rC + rC.^2.*(-C1W./rC).^m1.*NB;
M10 = rC.*G1 + G4./rC;
M11 = rC.*G2 + G5./rC;
M12=G3.*rC + G6./rC + rC.^2.*(-C1W0./rC).^m0.*NC;
M13 = rD.*G1 + G4./rD;
M14 = rD.*G2 + G5./rD;
M15=G3.*rD + G6./rD + rD.^2.*(-C1W0./rD).^m0.*NC;
M16 = H1./rD + H4.*rD;
M17 = H2./rD + H5.*rD;
M18 = H3./rD + H6.*rD;
%constant for total stresses at the interface
sigmaB = -(M1.*M10.*M14.*sigmaA - M1.*M11.*M13.*sigmaA + M1.*M11.*M17.*sigmaA ...
    + M1.*M13.*M7.*sigmaA - M1.*M17.*M7.*sigmaA - M10.*M16.*M4.*sigmaE + M10.*M14.*M3 ...
    - M10.*M14.*M6 + M10.*M15.*M4 - M10.*M18.*M4 - M11.*M13.*M3 + M11.*M13.*M6 + M11.*M17.*M3 ...
    - M11.*M17.*M6 - M12.*M13.*M4 + M12.*M17.*M4 + M13.*M3.*M7 + M13.*M4.*M9 - M13.*M6.*M7 ...
    - M17.*M3.*M7 - M17.*M4.*M9 + M17.*M6.*M7)./(M10.*M14.*M2 - M10.*M14.*M5 - M11.*M13.*M2 ...
    + M11.*M13.*M5 + M11.*M17.*M2 - M11.*M17.*M5 + M13.*M2.*M7 + M13.*M4.*M8 - M13.*M5.*M7 ...
    - M17.*M2.*M7 - M17.*M4.*M8 + M17.*M5.*M7);
sigmaC = (M1.*M13.*M8.*sigmaA - M1.*M17.*M8.*sigmaA + M10.*M16.*M2.*sigmaE ...
    - M10.*M16.*M5.*sigmaE - M10.*M15.*M2 + M10.*M15.*M5 + M10.*M18.*M2 ...
    - M10.*M18.*M5 + M12.*M13.*M2 - M12.*M13.*M5 - M12.*M17.*M2 + M12.*M17.*M5 ...
    - M13.*M2.*M9 + M13.*M3.*M8 + M13.*M5.*M9 - M13.*M6.*M8 + M17.*M2.*M9 ...
    - M17.*M3.*M8 - M17.*M5.*M9 + M17.*M6.*M8)./(M10.*M14.*M2 - M10.*M14.*M5 ...
    - M11.*M13.*M2 + M11.*M13.*M5 + M11.*M17.*M2 - M11.*M17.*M5 + M13.*M2.*M7 ...
    + M13.*M4.*M8 - M13.*M5.*M7 - M17.*M2.*M7 - M17.*M4.*M8 + M17.*M5.*M7);
sigmaD = -(M1.*M14.*M8.*sigmaA + M11.*M16.*M2.*sigmaE - M11.*M16.*M5.*sigmaE ...
    - M16.*M2.*M7.*sigmaE - M16.*M4.*M8.*sigmaE + M16.*M5.*M7.*sigmaE - M11.*M15.*M2 ...
    + M11.*M15.*M5 + M11.*M18.*M2 - M11.*M18.*M5 + M12.*M14.*M2 - M12.*M14.*M5 ...
    - M14.*M2.*M9 + M14.*M3.*M8 + M14.*M5.*M9 - M14.*M6.*M8 + M15.*M2.*M7 + M15.*M4.*M8 ...
    - M15.*M5.*M7 - M18.*M2.*M7 - M18.*M4.*M8 + M18.*M5.*M7)./(M10.*M14.*M2 - M10.*M14.*M5 ...
    - M11.*M13.*M2 + M11.*M13.*M5 + M11.*M17.*M2 - M11.*M17.*M5 + M13.*M2.*M7 + M13.*M4.*M8 ...
    - M13.*M5.*M7 - M17.*M2.*M7 - M17.*M4.*M8 + M17.*M5.*M7);
%total stresses at the interface
C1A = E1.*sigmaA + E2.*sigmaB + E3;
C2A = E4.*sigmaA + E5.*sigmaB + E6;
%integration constant of displacement general solution for region A
C1B = F1.*sigmaC + F2.*sigmaB + F3;
C2B = F4.*sigmaC + F5.*sigmaB + F6;
%integration constant of displacement general solution for region B
C1C = G1.*sigmaD + G2.*sigmaC + G3;
C2C = G4.*sigmaD + G5.*sigmaC + G6;
%integration constant of displacement general solution for region C
C2D = H4.*sigmaE + H5.*sigmaD + H6;
C1D = H1.*sigmaE + H2.*sigmaD + H3;
%integration constant of displacement general solution for region D

%%%%%%%%%%coordinate discretization%%%%%%%%%%
r1=rA:0.002:rB;
r2=rB:0.002:rC;
r3=rC:0.002:rD;
r4=rD:0.002:rE;
%%%%%%%%%%displacement solution%%%%%%%%%%
urA = r1.^NA1.*C2A + r1.^(-NA1).*C1A+(((-C1W./r1).^m1).*(r1.^2)).*NA2;
%displacement solution for region A under elasticity
urB=r2.*C1B + C2B./r2 + r2.^2.*(-C1W./r2).^m1.*NB;
%displacement solution for region B under elasticity
urC=r3.*C1C + C2C./r3 + r3.^2.*(-C1W0./r3).^m0.*NC;
%displacement solution for region C under elasticity
urD = C1D./r4 + C2D.*r4;
%displacement solution for region D under elasticity

%%%%%%%%%%first-order derivatives of displacement solutions%%%%%%%%%%
urA1 =r1.^NA1.*NA1.*C2A./r1 - r1.^(-NA1).*NA1.*C1A./r1 ...
    -(-C1W./r1).^m1.*m1.*r1.*NA2 + 2.*(-C1W./r1).^m1.*r1.*NA2;
%first-order derivatives of displacement solution A under elasticity
urB1 =C1B - C2B./r2.^2 + 2.*r2.*(-C1W./r2).^m1.*NB - r2.*(-C1W./r2).^m1.*m1.*NB;
%first-order derivatives of displacement solution B under elasticity
urC1 =C1C - C2C./r3.^2 + 2.*r3.*(-C1W0./r3).^m0.*NC - r3.*(-C1W0./r3).^m0.*m0.*NC;
%first-order derivatives of displacement solution C under elasticity
urD1 = -C1D./r4.^2 + C2D;
%first-order derivatives of displacement solution D under elasticity

%%%%%%%%%%strain solutions%%%%%%%%%%
varepsilonrA = urA1;
%radial strain of the surrounding rock in region A under elasticity
varepsilonthetaA = urA./r1;
%hoop strain of the surrounding rock in region A under elasticity
varepsilonrB= urB1;
%radial strain of the surrounding rock in region B under elasticity
varepsilonthetaB = urB./r2;
%hoop strain of the surrounding rock in region B under elasticity
varepsilonrC = urC1;
%radial strain of the surrounding rock in region C under elasticity
varepsilonthetaC = urC./r3;
%hoop strain of the surrounding rock in region C under elasticity
varepsilonrD = urD1;
%radial strain of the surrounding rock in region D under elasticity
varepsilonthetaD = urD./r4;
%hoop strain of the surrounding rock in region D under elasticity

%%%%%%%%%%pore pressure solutions%%%%%%%%%%
Apw1 = -r1.*(-C1W./r1).^m1.*gammaw./((m1 - 1).*k1) + C2W;
%pore pressure solution in the surrounding rock for region A
Apw2 = -r2.*(-C1W./r2).^m1.*gammaw./((m1 - 1).*k1) + C2W;
%pore pressure solution in the surrounding rock for region B
Bpw = -r3.*(-C1W0./r3).^m0.*gammaw./((m0 - 1).*k0) + C2W0;
%pore pressure solution in the surrounding rock for region C

%%%%%%%%%%total stresses solutions%%%%%%%%%%
sigmarA = 2.*G.*(1 + mu).* varepsilonrA + 2.*G.*mu.*varepsilonthetaA + sigmaE - pwE + Apw1;
%elastic total radial stress solution of the surrounding rock in region A
sigmathetaA = 2.*G.*(1 + mu).* varepsilonthetaA + 2.*G.*mu.*varepsilonrA + sigmaE - pwE + Apw1;
%elastic total hoop stress solution of the surrounding rock in region A
sigmarB = 2.*G.*(1 + mu).*varepsilonrB + 2.*G.*mu.*varepsilonthetaB + sigmaE - pwE + Apw2;
%elastic total radial stress solution of the surrounding rock in region B
sigmathetaB = 2.*G.*(1 + mu).* varepsilonthetaB + 2.*G.*mu.* varepsilonrB + sigmaE - pwE + Apw2;
%elastic total hoop stress solution of the surrounding rock in region B
sigmarC = 2.*G0.*(1 + mu0).* varepsilonrC + 2.*G0.*mu0.* varepsilonthetaC + sigmaE - pwE + Bpw;
%elastic total radial stress solution of the surrounding rock in region C
sigmathetaC = 2.*G0.*(1 + mu0).* varepsilonthetaC + 2.*G0.*mu0.* varepsilonrC + sigmaE - pwE + Bpw;
%elastic total hoop stress solution of the surrounding rock in region C
sigmarD = 2.*G0.*(1 + mu0).* varepsilonrD + 2.*G0.*mu0.* varepsilonthetaD + sigmaE;
%elastic total radial stress solution of the surrounding rock in region D
sigmathetaD = 2.*G0.*(1 + mu0).* varepsilonthetaD + 2.*G0.*mu0.* varepsilonrD + sigmaE;
%elastic total hoop stress solution of the surrounding rock in region D
%%%%%%%%%%pore pressure solutions%%%%%%%%%%
sigmarAYX= sigmarA - Apw1;
%elastic effective radial stress solution of the surrounding rock in region A
sigmathetaAYX = sigmathetaA - Apw1;
%elastic effective hoop stress solution of the surrounding rock in region A
sigmarBYX = sigmarB - Apw2;
% elastic effective radial stress solution of the surrounding rock in region B 
sigmathetaBYX = sigmathetaB - Apw2;
% elastic effective hoop stress solution of the surrounding rock in region B
sigmarCYX = sigmarC - Bpw;
% elastic effective radial stress solution of the surrounding rock in region C
sigmathetaCYX = sigmathetaC - Bpw;
% elastic effective hoop stress solution of the surrounding rock in region C
sigmarDYX = sigmarD - pwE;
% elastic effective radial stress solution of the surrounding rock in region D
sigmathetaDYX = sigmathetaD - pwE;
% elastic effective hoop stress solution of the surrounding rock in region D

%%%%%%%%%%Solving for the radius of the brittle-plastic region%%%%%%%%%%
JMCA=sigmathetaAYX-k.*sigmarAYX;
MCA=(k-1).*a;
resA=JMCA-MCA;
resA=resA(:);
resAA=min(abs(resA));
r1=r1(:);
mA=length(resA);
for i=1:1:mA
    if (resA(i)==resAA|resA(i)==-resAA);
        break
    end
end

RE = r1(i);
%RE is radius of the brittle-plastic region
if RE == rA
    RE=rA
    urA=urA(:);
    urAMAX(V)=urA(1);
end
if RE==rB
    JMCB=sigmathetaBYX-k.*sigmarBYX;
    MCB=(k-1).*a;
    resB=JMCB-MCB;
    resB=resB(:);
    resBB=min(abs(resB));
    r2=r2(:);
    mB=length(resB);
    for i=1:1:mB
        if (resB(i)==resBB|resB(i)==-resBB);
            break
        end
    end
    RE=r2(i);
end
if RE==rC
    JMCC=sigmathetaCYX-k00.*sigmarCYX;
    MCC=(k00-1).*a0;
    resC=JMCC-MCC;
    resC=resC(:);
    resCC=min(abs(resC));
    r3=r3(:);
    mC=length(resC);
    for i=1:1:mC
        if (resC(i)==resCC|resC(i)==-resCC);
            break
        end
    end
    RE=r3(i);
end

%%%%%%%%%%case I%%%%%%%%%
if RE ~= rA  & ( RE<rB | RE==rB )
uRA=urA(i);
%displacement at r=Re
uRA1=urA1(i);
%first-order derivatives of displacement at r=Re
sigmarAYXR=sigmarAYX(i)+((Eb.*uRA1.*Y)./(S.*DS));
%effectice radial stress at r=Re
ZRA=sigmabf.*Y.*nb/(2.*3.14.*DS.*RE);
PR=sigmarAYXR+(ZRA./1);
PA=sigmaA-pwA+(ZRA./(rA./RE));
%boundary conditions
C1RA = ((-ZRA - PR + PA).*rA + RE.*ZRA).*(kr - 1)./((rA./RE).^kr.*RE - rA);
C2RA = (RE.*(ZRA + PR).*(rA./RE).^kr - RE.*ZRA - PA.*rA)./((rA./RE).^kr.*RE - rA);
%integration constants
%%%%%%%%%%brittle-plastic region coordinate discretization%%%%%%%%%%
rR=rA:0.1:RE;
rho=rR./RE;
sigmarRYXA= (((C1RA.*rho.^kr)./(kr - 1)) - ZRA)./rho + C2RA;
sigmathetarRYXA= sigmarRYXA.*kr + (kr - 1).*ar;
rhoR=rA./RE;
keci=1.0;
%jump factor for brittle-plasticity
[tRA,yRA]=ode45(@znxpfun1R,[1:-0.0001:rhoR],[uRA,(keci*uRA1*RE)]);
r0=tRA.*RE;
%radial coordinates of the brittle plastic region

subplot(2,2,1)
h1(1)=plot(r0,yRA(:,1).*1000,'r','LineWidth',1.5);
%yRA is the brittle-plastic displacement field in region A under case I
hold on
h1(2)=plot(r1,urA.*1000,'k','LineWidth',1.5);
hold on
plot(r2,urB.*1000,'k','LineWidth',1.5);
legend([h1(1),h1(2)],'brittle-plastic displacement','elastic displacement');
xlabel('radial distance to circular tunnel axis/m');
ylabel('radial displacement/mm');
title('Displacement field');


subplot(2,2,2)
h2(1)=plot(rR,sigmarRYXA.*1e-6,'r','LineWidth',1.5);
hold on
h2(2)=plot(r1,sigmarAYX.*1e-6,'k','LineWidth',1.5);
hold on
plot(r2,sigmarBYX.*1e-6,'k','LineWidth',1.5);
hold on 
plot(rR,sigmathetarRYXA.*1e-6,'r','LineWidth',1.5);
hold on
plot(r1,sigmathetaA.*1e-6,'k','LineWidth',1.5);
hold on
plot(r2,sigmathetaB.*1e-6,'k','LineWidth',1.5);
hold on 
legend([h2(1),h2(2)],'brittle-plastic effective stress','elastic effective stress');
xlabel('radial distance to circular tunnel axis/m');
ylabel('effective stress/MPa');
title('Effective stress field, radial (below),hoop (above)');

subplot(2,2,3)
plot(r1,Apw1.*1e-6,'b','LineWidth',1.5);
hold on
plot(r2,Apw2.*1e-6,'b','LineWidth',1.5);
hold on
plot(r3,Bpw.*1e-6,'b','LineWidth',1.5);
hold on
xlabel('radial distance to circular tunnel axis/m');
ylabel('pore pressure/MPa');
title('Pore pressure field');

length(r0);
rROCKBOLT=ones(1,length(r0));
rROCKBOLT=-1.0.*sigmabf.*rROCKBOLT.*1e-6;
subplot(2,2,4)
h3(1)=plot(r1,((Eb.*urA1)./(1.0e6)),'k','LineWidth',1.5);
hold on
%Eb.*urA1 is axial stress on the rockbolts
h3(2)=plot(r0,rROCKBOLT,'r','LineWidth',1.5);
hold on
legend([h3(1),h3(2)],'in brittle-plastic region','under elasticity')
xlabel('radial distance to circular tunnel axis/m');
ylabel('Axial stress/MPa');
title('Axial stress on the rockbolts');


%%%%%%%%%%Code used when extracting parameters to draw a graph%%%%%%%%%
% r0 = r0(:);
% r0=r0';
% r0=fliplr(r0);
% r0=r0';
% AyRA=yRA(:,1);
% AyRA=AyRA';
% AyRA=fliplr(AyRA);
% AyRA=AyRA';
% cigmaBBB=Eb.*urA1;
% cigmaBBB=cigmaBBB';
% rR=rR';
% sigmarRYXA=sigmarRYXA';
% sigmathetarRYXA=sigmathetarRYXA';
% r2=r2';
% r3=r3';
% urA=urA';
% urB=urB';
% urC=urC';
% sigmathetaAYX=sigmathetaAYX';
% sigmarAYX=sigmarAYX';
% sigmathetaBYX=sigmathetaBYX';
% sigmarBYX=sigmarBYX';
% sigmathetaCYX=sigmathetaCYX';
% sigmarCYX=sigmarCYX';
% Apw1=Apw1';
% Apw2=Apw2';
end

%%%%%%%%%%case II%%%%%%%%%%
if (rB<RE & RE<rC) | RE==rC
uRB=urB(i);
%displacement at r=Re under case II
uRB1=urB1(i);
%first-order derivatives of displacement at r=Re under case II
sigmarBYXR=sigmarBYX(i);
%effectice radial stress at r=Re under case II
PRB=sigmarBYXR;
%with the help of the plastic condition the stress at r=rB in region A is in fact the force at rho=1
REB=rB;
%boundary condition
ZRB=sigmabf.*Y.*nb/(2.*3.14.*DS.*rB);
PAA=sigmaA-pwA+(ZRB./(rA./REB));
mAR=length(sigmarAYX);
uRA1rB=urA1(mAR);
sigmarAYXR=sigmarAYX(mAR)+((Eb.*uRA1rB.*Y)./(S.*DS));
PRA=sigmarAYXR+(ZRB./1);
C1RBA = ((-ZRB - PRA + PAA).*rA + REB.*ZRB).*(kr - 1)./((rA./REB).^kr.*REB - rA);
C2RBA = (REB.*(ZRB + PRA).*(rA./REB).^kr - REB.*ZRB - PAA.*rA)./((rA./REB).^kr.*REB - rA);
%integration constants
sigmarRYXBA= (C1RBA.*(1.0).^kr./(kr - 1) - ZRB)./(1.0) + C2RBA;
PBB=sigmarRYXBA;
C1RB = (PRB*(rB/RE)^(-1 + kr) - PBB)/((rB/RE)^(-1 + kr) - 1);
C2RB = (PBB - PRB)/((rB/RE)^(-1 + kr) - 1);
%integration constants
%%%%%%%%%%brittle-plastic region coordinate discretization%%%%%%%%%%
rRB=rB:0.01:RE;
rhoB=rRB./RE;

sigmarRYXB=C1RB + C2RB.*rhoB.^(-1 + kr);
sigmathetarRYXB= sigmarRYXB.*kr + (kr - 1).*ar;
%effective brittle-plastic stress solution in region B under case II

%%%%%%%%%%brittle-plastic region coordinate discretization%%%%%%%%%%
rRA=rA:0.01:rB;
rhoA=rRA./rB;
sigmarRYXBA= (C1RBA.*(rhoA).^kr./(kr - 1) - ZRB)./(rhoA) + C2RBA;
sigmathetarRYXBA= sigmarRYXBA.*kr + (kr - 1).*ar;
%effective brittle-plastic stress solution in region A under case II

%%%%%%%%%displacement under case II%%%%%%%%%%
keci=1.0;
rhoRB=rB./RE;
[tRB,yRB]=ode45(@znxpfun1RB,[1:-0.0001:rhoRB],[uRB,(keci*uRB1*RE)]); 
uRBA=yRB(length(yRB(:,1)),1);
uRBA1=yRB(length(yRB(:,2)),2);
rhoRBA=rA./rB;
[tRBA,yRBA]=ode45(@znxpfun1RBA,[1:-0.0001:rhoRBA],[uRBA,(uRBA1*rB)]);
uRBAMAX=yRBA(length(yRBA(:,1)),1);
r0RB=tRB.*RE;
r0RBA=tRBA.*rB;

subplot(2,2,1)
h1(1)=plot (rRB,sigmarRYXB.*1e-6,'r','LineWidth',1.5);
hold on
plot (rRA,sigmarRYXBA.*1e-6,'r','LineWidth',1.5);
hold on
h1(2)=plot (rRB,sigmathetarRYXB.*1e-6,'k','LineWidth',1.5);
hold on
plot (rRA,sigmathetarRYXBA.*1e-6,'k','LineWidth',1.5);
hold on
legend([h1(1),h1(2)],'effective hoop stress','effective radial stress');
xlabel('radial distance to circular tunnel axis/m');
ylabel('effective stress/MPa');
title('Effective stress field in brittle-plastic region');

subplot(2,2,2)
plot(r0RB,yRB(:,1).*1000,'r','LineWidth',1.5);
hold on
plot(r0RBA,yRBA(:,1).*1000,'r','LineWidth',1.5);
hold on
xlabel('radial distance to circular tunnel axis/m');
ylabel('displacement/MPa');
title('Displacement field in brittle-plastic region');

subplot(2,2,3)
plot(r1,Apw1.*1e-6,'b','LineWidth',1.5);
hold on
plot(r2,Apw2.*1e-6,'b','LineWidth',1.5);
hold on
plot(r3,Bpw.*1e-6,'b','LineWidth',1.5);
hold on
xlabel('radial distance to circular tunnel axis/m');
ylabel('pore pressure/MPa');
title('Pore pressure field');

rROCKBOLT=ones(1,length(r0RBA));
rROCKBOLT=-1.0.*sigmabf.*rROCKBOLT.*1e-6;
subplot(2,2,4)
plot(r0RBA,rROCKBOLT,'r','LineWidth',1.5);
hold on
xlabel('radial distance to circular tunnel axis/m');
ylabel('Axial stress/MPa');
title('Axial stress on the rockbolts');
end

%%%%%%%%%%case III%%%%%%%%%
if (rC<RE & RE<rD) | RE==rD
uRC=urC(i);
%displacement at r=Re under case III
uRC1=urC1(i);
%first-order derivatives of displacement at r=Re under case III
sigmarCYXR=sigmarCYX(i);
%effectice radial stress at r=Re under case III
PRC=sigmarCYXR;
REC=rC;
REB=rB;
ZRB=sigmabf.*Y.*nb/(2.*3.14.*DS.*rB);
PAA=sigmaA-pwA+(ZRB./(rA./REB));
mAR=length(sigmarAYX);
uRA1rB=urA1(mAR);
sigmarAYXR=sigmarAYX(mAR)+((Eb.*uRA1rB.*Y)./(S.*DS));
PRA=sigmarAYXR+(ZRB./1);
C1RCBA = ((-ZRB - PRA + PAA).*rA + REB.*ZRB).*(kr - 1)./((rA./REB).^kr.*REB - rA);
C2RCBA = (REB.*(ZRB + PRA).*(rA./REB).^kr - REB.*ZRB - PAA.*rA)./((rA./REB).^kr.*REB - rA);
%integration constants
sigmarRYXCBA= (C1RCBA.*(1.0).^kr./(kr - 1) - ZRB)./(1.0) + C2RCBA;
PBB=sigmarRYXCBA;
mBR=length(sigmarBYX);
sigmarBYXR=sigmarBYX(mBR);
PCC=sigmarBYXR;
C1RCB = (PCC*(rB/REC)^(-1 + kr) - PBB)/((rB/REC)^(-1 + kr) - 1);
C2RCB = (PBB - PCC)/((rB/REC)^(-1 + kr) - 1);
C1RC = (PRC*(rC/RE)^(-1 + kr0) - PCC)/((rC/RE)^(-1 + kr0) - 1);
C2RC = (PCC - PRC)/((rC/RE)^(-1 + kr0) - 1);
%integration constants

%%%%%%%%%%brittle-plastic region coordinate discretization%%%%%%%%%%
rRC=rC:0.01:RE;
rhoC=rRC./RE;
sigmarRYXC=C1RC + C2RC.*rhoC.^(-1 + kr0);
sigmathetarRYXC= sigmarRYXC.*kr0 + (kr0 - 1).*ar0;
%effective brittle-plastic stress solution in region C under case III
%%%%%%%%%%brittle-plastic region coordinate discretization%%%%%%%%%%
rRB=rB:0.01:rC;
rhoB=rRB./rC;
sigmarRYXCB=C1RCB + C2RCB.*rhoB.^(-1 + kr);
sigmathetarRYXCB= sigmarRYXCB.*kr + (kr - 1).*ar;
%effective brittle-plastic stress solution in region B under case III
%%%%%%%%%%brittle-plastic region coordinate discretization%%%%%%%%%%
rRA=rA:0.01:rB;
rhoA=rRA./rB;
sigmarRYXCBA= (C1RCBA.*(rhoA).^kr./(kr - 1) - ZRB)./(rhoA) + C2RCBA;
sigmathetarRYXCBA= sigmarRYXCBA.*kr + (kr - 1).*ar;
%effective brittle-plastic stress solution in region A under case III
keci0=1.0;
rhoRC=rC./RE;
[tRC,yRC]=ode45(@znxpfun1RC,[1:-0.0001:rhoRC],[uRC,(keci0*uRC1*RE)]); 
uRCB=yRC(length(yRC(:,1)),1);
uRCB1=yRC(length(yRC(:,2)),2);
rhoRCB=rB./rC;
[tRCB,yRCB]=ode45(@znxpfun1RCB,[1:-0.0001:rhoRCB],[uRCB,(uRCB1*rC)]); 
uRCBA=yRCB(length(yRCB(:,1)),1);
uRCBA1=yRCB(length(yRCB(:,2)),2);
rhoRCBA=rA./rB;
[tRCBA,yRCBA]=ode45(@znxpfun1RCBA,[1:-0.0001:rhoRCBA],[uRCBA,(uRCBA1*rB)]); 
r0RC=tRC.*RE;
r0RCB=tRCB.*rC;
r0RCBA=tRCBA.*rB;

subplot(2,2,1)
h1(1)=plot (rRB,sigmarRYXCB.*1e-6,'r','LineWidth',1.5);
hold on
plot (rRA,sigmarRYXCBA.*1e-6,'r','LineWidth',1.5);
hold on
h1(2)=plot (rRB,sigmathetarRYXCB.*1e-6,'k','LineWidth',1.5);
hold on
plot (rRA,sigmathetarRYXCBA.*1e-6,'k','LineWidth',1.5);
hold on
plot (rRC,sigmarRYXC.*1e-6,'r','LineWidth',1.5);
hold on
plot (rRC,sigmathetarRYXC.*1e-6,'k','LineWidth',1.5);
hold on
legend([h1(1),h1(2)],'effective radial stress','effective hoop stress');
xlabel('radial distance to circular tunnel axis/m');
ylabel('effective brittle-plastic stress/MPa');
title('Effective stress field in brittle-plastic region');


subplot(2,2,2)
plot(r0RC,yRC(:,1).*1000,'r','LineWidth',1.5);
hold on
plot(r0RCB,yRCB(:,1).*1000,'r','LineWidth',1.5);
hold on
plot(r0RCBA,yRCBA(:,1).*1000,'r','LineWidth',1.5);
hold on
xlabel('radial distance to circular tunnel axis/m');
ylabel('brittle-plastic displacement/mm');
title('Displacement field in brittle-plastic region');

subplot(2,2,3)
plot(r1,Apw1.*1e-6,'b','LineWidth',1.5);
hold on
plot(r2,Apw2.*1e-6,'b','LineWidth',1.5);
hold on
plot(r3,Bpw.*1e-6,'b','LineWidth',1.5);
hold on
xlabel('radial distance to circular tunnel axis/m');
ylabel('pore pressure/MPa');
title('Pore pressure field');

rROCKBOLT=ones(1,length(r0RCBA));
rROCKBOLT=-1.0.*sigmabf.*rROCKBOLT.*1e-6;
subplot(2,2,4)
plot(r0RCBA,rROCKBOLT,'r','LineWidth',1.5);
hold on
xlabel('radial distance to circular tunnel axis/m');
ylabel('Axial stress/MPa');
title('Axial stress on the rockbolts');
end


function ydotRC=znxpfun1RC(tRC,yRC)
global v0 kPsi0 RE G0 kr0   C2RC  sigmarRC1 sigmathiRC1

sigmarRC1=C2RC*tRC^(-1 + kr0)*(-1 + kr0)/tRC;
sigmathiRC1=(C2RC*tRC^(-1 + kr0)*(-1 + kr0)/tRC)*kr0;

ydotRC=[yRC(2);...
-((-1/tRC)*yRC(2)*(kPsi0-1))-...
(RE*((v0*kPsi0-kPsi0-v0)/(2*G0))*(sigmathiRC1))-...
(RE*((1-v0+v0*kPsi0)/(2*G0))*(sigmarRC1))];
end
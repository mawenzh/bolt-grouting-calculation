function ydotRB=znxpfun1RB(tRB,yRB)
global v kPsi kr G ZRB C1RB C2RB sigmarRB1 sigmathiRB1 RE PAA REB PRB PBB

sigmarRB1=C2RB*tRB^(-1 + kr)*(-1 + kr)/tRB;
sigmathiRB1=(C2RB*tRB^(-1 + kr)*(-1 + kr)/tRB)*kr;
ydotRB=[yRB(2);...
-((-1/tRB)*yRB(2)*(kPsi-1))-...
(RE*((v*kPsi-kPsi-v)/(2*G))*(sigmathiRB1))-...
(RE*((1-v+v*kPsi)/(2*G))*(sigmarRB1))];
end
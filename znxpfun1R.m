function ydotR=znxpfun1R(t,y)
global C1RA C2RA kr ZRA kPsi G RE v 

sigmarRA1=C1RA*t^kr*kr/((kr - 1)*t^2)-(C1RA*t^kr/(kr - 1) - ZRA)/t^2;

sigmathiRA1=(C1RA*t^kr*kr/((kr - 1)*t^2)-(C1RA*t^kr/(kr - 1) - ZRA)/t^2)*kr;

ydotR=[y(2);...
-((-1/t)*y(2)*(kPsi-1))-...
(RE*((v*kPsi-kPsi-v)/(2*G))*(sigmathiRA1))-...
(RE*((1-v+v*kPsi)/(2*G))*(sigmarRA1))];
end
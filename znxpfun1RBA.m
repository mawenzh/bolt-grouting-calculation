function ydotRBA=znxpfun1RBA(tRBA,yRBA)
global v kPsi rB G  kr sigmarRBA1 sigmathiRBA1 C1RBA ZRB

sigmarRBA1=C1RBA*tRBA^kr*kr/((kr - 1)*tRBA^2)-(C1RBA*tRBA^kr/(kr - 1) - ZRB)/tRBA^2;

sigmathiRBA1=(C1RBA*tRBA^kr*kr/((kr - 1)*tRBA^2)-(C1RBA*tRBA^kr/(kr - 1) - ZRB)/tRBA^2)*kr;

ydotRBA=[yRBA(2);...
-((-1/tRBA)*yRBA(2)*(kPsi-1))-...
(rB*((v*kPsi-kPsi-v)/(2*G))*(sigmathiRBA1))-...
(rB*((1-v+v*kPsi)/(2*G))*(sigmarRBA1))];
end
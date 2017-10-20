function [gradA,gradB]= gradCalc(gamma,hs,alpha,gs,p_A,p_B,x)
gradA = gs./(gamma(:,2)+hs*alpha*p_B+gs*p_A)+hs*alpha./(gamma(:,1)+hs*alpha*p_A+gs*p_B)-hs*alpha./(gamma(:,1)+p_A*hs*alpha);
gradA=gradA.*x;
gradA=gradA';
gradB = gs./(gamma(:,1)+hs*alpha*p_A+gs*p_B)+hs*alpha./(gamma(:,2)+hs*alpha*p_B+gs*p_A)-hs*alpha./(gamma(:,2)+p_B*hs*alpha);
gradB=gradB.*x;
gradB=gradB';
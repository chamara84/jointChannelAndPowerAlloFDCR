function [funcVal] = functionEval(gamma,hs,alpha,gs,p_A,p_B,x)

funcVal = x*(log(gamma(:,2)+alpha*hs*p_B+gs*p_A))-x*log(gamma(:,2) ...
    +alpha*hs*p_B)+x*(log(gamma(:,1)+alpha*hs*p_A+gs*p_B))-x*log(gamma(:,1)+alpha*hs*p_A);
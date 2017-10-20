function [obj,feasibility]=calculateObjectiveFD(boxesUpperCorners,GammaInitlinear,hs,alpha,gainVecSUTrnsSURecv,gainVecSUTrnsPURecv,IntThreshold,SUIntAtPU,x_k,RxSensitivity,Pmax)

numCols=size(boxesUpperCorners,2);
obj=zeros(1,numCols);
feasibility = zeros(1,numCols);

for i=1:numCols
obj(1,i)= prod([(GammaInitlinear(1:n)+boxesUpperCorners(1:n,i)*hs*alpha+boxesUpperCorners(n+1:end,i)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+boxesUpperCorners(1:n,i)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+boxesUpperCorners(n+1:end,i)*hs*alpha+boxesUpperCorners(1:n,i)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+boxesUpperCorners(n+1:end,i)*hs*alpha)]);
    intAtPU=boxesUpperCorners(1:n,i)*gainVecSUTrnsPURecv(1,1)+boxesUpperCorners(n+1:end,i)*gainVecSUTrnsPURecv(2,1);
if sum(boxesUpperCorners(1:n,i))<=Pmax && sum(boxesUpperCorners(n+1:end,i))<=Pmax && all(intAtPU<=(IntThreshold-SUIntAtPU))
            if all(boxesUpperCorners(n+1:end,i)*gainVecSUTrnsSURecv./(GammaInitlinear(1:n)+boxesUpperCorners(1:n,i)*hs*alpha)>=RxSensitivity.*x_k(n+1:end,1)) && ...
                    all(boxesUpperCorners(1:n,i)*gainVecSUTrnsSURecv./(GammaInitlinear(n+1:end)+boxesUpperCorners(n+1:end,i)*hs*alpha)>=RxSensitivity.*x_k(1:n,1)) 
                
                feasibility(1,i)=1;
            else
                feasibility(1,i)=0;
            end
end
end
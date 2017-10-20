function [obj,feasibility]=feasibilityFD(boxesUpperCorners,GammaInitlinear,hs,alpha,gainVecSUTrnsSURecv,gainVecSUTrnsPURecv,IntThreshold,SUIntAtPU,x_k,RxSensitivity,Pmax)


obj= prod([(GammaInitlinear(1:n)+boxesUpperCorners(1:n)*hs*alpha+boxesUpperCorners(n+1:end)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+boxesUpperCorners(1:n)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+boxesUpperCorners(n+1:end)*hs*alpha+boxesUpperCorners(1:n)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+boxesUpperCorners(n+1:end)*hs*alpha)]);
    intAtPU=boxesUpperCorners(1:n)*gainVecSUTrnsPURecv(1,1)+boxesUpperCorners(n+1:end)*gainVecSUTrnsPURecv(2,1);
if sum(boxesUpperCorners(1:n))<=Pmax && sum(boxesUpperCorners(n+1:end))<=Pmax && all(intAtPU<=(IntThreshold-SUIntAtPU))
            if all(boxesUpperCorners(n+1:end)*gainVecSUTrnsSURecv./(GammaInitlinear(1:n)+boxesUpperCorners(1:n)*hs*alpha)>=RxSensitivity.*x_k(n+1:end,1)) && ...
                    all(boxesUpperCorners(1:n)*gainVecSUTrnsSURecv./(GammaInitlinear(n+1:end)+boxesUpperCorners(n+1:end)*hs*alpha)>=RxSensitivity.*x_k(1:n,1)) 
                
                feasibility=true;
            else
                feasibility=false;
                obj=[];
            end
end
end
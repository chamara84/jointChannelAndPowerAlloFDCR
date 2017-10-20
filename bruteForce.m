%% brute force algorithm
clear all
%clc

numChannels =2;
n=numChannels; 
SUIndex = 1;
gainVecSUTrnsSURecv = 0.6;
gainVecSUTrnsPURecv = [ 0.2947;0.2307];%rand(2,1)*0.3;
Xinit = ones(numChannels*2,1);
epsilon=0.1;
GammaInit = [0.0581 ,   0.0928;
              0.0580 ,   0.0017];
%              0.0121 ,   0.0863;
%              0.0484 ,   0.0845;
%              0.0209 ,   0.0552];%rand(2,numChannels)'*0.1;
Pmax = 1;
bestAllocation = zeros(2*n,1);
maxiterations = 100;
hs = 0.8;
alpha = 0.001;
SUIntAtPU = [ 0.0504;0.0026];%0.0492;0.0290;0.0040];
IntThreshold = ones(numChannels,1)*0.1;
RxSensitivity  = 0.001;
GammaInitlinear=[GammaInit(:,1);GammaInit(:,2)];
P_A=0:0.001:Pmax;
P_B=0:0.001:Pmax;
cbv=-inf;
%channel usage
X=dec2bin(0:2^(2*n)-1); %X_A(1,:) gives the first sequence
powerCombinations=[];
utility=[];
for channelAllocationIndex=1:2^(2*n)
    channelAllocation=X(channelAllocationIndex,:);
    powerCombinations=[];
    xVec = [];
    for index=1:2*n
        xVec = [xVec,str2num(channelAllocation(index))];
        
        if channelAllocation(index)=='1'
            if ~isempty(powerCombinations);
                powerCombinations=combvec(P_A,powerCombinations);
            else
                powerCombinations=P_A;
            end
        else
            if ~isempty(powerCombinations)
             powerCombinations=combvec(0,powerCombinations);
            else
                powerCombinations=0;
            end
        end
    end
    disp(xVec)
    for combinationIndex=1:size(powerCombinations,2)
        power=powerCombinations(:,combinationIndex);
        intAtPU=power(1:n)*gainVecSUTrnsPURecv(1,1)+power(n+1:end)*gainVecSUTrnsPURecv(2,1);
        objFuncUpperbound= prod([(GammaInitlinear(1:n)+power(1:n)*hs*alpha+power(n+1:end)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+power(1:n)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+power(n+1:end)*hs*alpha+power(1:n)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+power(n+1:end)*hs*alpha)]);
        if sum(power(1:n))<=Pmax && sum(power(n+1:end))<=Pmax && all(intAtPU<=(IntThreshold-SUIntAtPU))
            if all(power(n+1:end)*gainVecSUTrnsSURecv./(GammaInitlinear(1:n)+power(1:n)*hs*alpha)>=RxSensitivity.*xVec(n+1:end)') && ...
                    all(power(1:n)*gainVecSUTrnsSURecv./(GammaInitlinear(n+1:end)+power(n+1:end)*hs*alpha)>=RxSensitivity.*xVec(1:n)') && objFuncUpperbound>cbv
                    cbv=objFuncUpperbound;
                    bestAllocation = power;
                
            end
        end
    
    
    end
end
                
           
            
            
            
      
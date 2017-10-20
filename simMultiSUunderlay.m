function [ rateAchievedbySUs,interferenceOnPU ] = simMultiSUunderlay( numSUPairs,numPUs,gridSize,timeSlots,phi,threshold,K,PMax,SUPairPos,PUPos)
%This function simulates an underlay network of decentralized SUs and
%cetralized PUs for a given number of time slots. 

frequency = 2.4e6; %transmit frequency
vecSUs = [];
Xinit = zeros(1,numPUs);
Pinit = zeros(1,numPUs);
%GammaInit = zeros(1,numPUs);
 %% set initial values for X,P,Gamma and the lagrangian multipliers
 %Lagrange multipliers
 betaQInit = 0;
 alphaInit = ones(1,numPUs);
 lambda = zeros(1,numPUs);
 betaXInit = 0;
 betaPInit = 0;
 %X initialize
 chanSet = randi(numPUs,1,K);
 Xinit(chanSet)=1;
 %P initialize
 Pinit(chanSet) = PMax/sum(Xinit);
 
 %% matrices to  collect power values and channel allocation values
 P_mat = zeros(numSUPairs,numPUs);
 X_mat = zeros(numSUPairs,numPUs);
 PUIntOnSU = zeros(numPUs,numSUPairs);
 %% parameters to collect QoS
 interferenceOnPU = zeros(timeSlots,numPUs);
 rateAchievedbySUs = zeros(timeSlots,numSUPairs);
 
 
 PUTxPower = PMax*rand(numPUs,1); % PU transmit power
 
    [gainMatSUTrnsSURecv,gainMatSUTrnsPURecv,gainMatSURecvPUTrns ] = chanGainCalc(SUPairPos,PUPos,frequency);
    for channel=1:numPUs
        PUIntOnSU(channel,:) = gainMatSURecvPUTrns(channel,:)*PUTxPower(channel,1);
        
    end
 
 %% create SU objects
    
    
 for SUIndex=1:numSUPairs
     %Gamma initialize
     GammaInit = optGamma( Xinit,alphaInit,gainMatSUTrnsSURecv(SUIndex,:),Pinit,betaQInit,SUindex);
     vecSUs = [vecSUs,SecondaryUser(SUIndex,gainMatSUTrnsSURecv(SUIndex,:),gainMatSUTrnsPURecv(SUIndex,:),Xinit,...
         Pinit,GammaInit,alphaInit,betaQInit,betaXInit,K,maxIter,step,lambda,betaPInit,Pmax,phi)];
     
 end
 
%% simulate for time slots
for timeSlot=1:timeSlots
    PUTxPower = PMax*rand(numPUs,1);
    [gainMatSUTrnsSURecv,gainMatSUTrnsPURecv,gainMatSURecvPUTrns ] = chanGainCalc(SUPairPos,PUPos,frequency);
    for channel=1:numPUs
        PUIntOnSU(channel,:) = gainMatSURecvPUTrns(channel,:)*PUTxPower(channel,1);
        
    end
    for SUIndex=1:numSUPairs
        % update channel gains
        vecSUs(SUIndex).gainVecSUTrnsSURecv = gainMatSUTrnsSURecv(SUIndex,:);
        vecSUs(SUIndex).gainVecSUTrnsPURecv = gainMatSUTrnsPURecv(SUIndex,:);
        
    end
    
       
    
    
    %% execute the optimization of SU parameters in parallel
    parfor index = 1: numSUPairs
        vecSUs(index).optXPGamma();
    end
    
       
    
    
    %% get the optimized parameters and find if the QoS is met and if the PU interference threshold is violated
    for SUIndex=1:numSUPairs
        P_mat(SUIndex,:) = vecSUs(SUIndex).P;
        X_mat(SUIndex,:) = vecSUs(SUIndex).X;
        
    end
    
    for SUIndex=1:numSUPairs
        rateAchievedbySUs(timeSlot,SUIndex) =  sum(log2(1+gainMatSUTrnsSURecv(SUIndex,SUIndex)*P_mat(SUIndex,:).*X_mat(SUIndex,:)...
            ./(gainMatSUTrnsSURecv(1:end ~= SUIndex,SUIndex)'*(P_mat(1:end ~= SUIndex,:).*X_mat(1:end ~= SUIndex,:))+PUIntOnSU(:,SUIndex)')));
    end
    
    
    for channel=1:numPUs
        interferenceOnPU(timeSlot,channel) = gainMatSUTrnsPURecv(:,channel)'*(P_mat(:,channel).*X_mat(:,channel));
    end
    %% calculate the lagrangian parameters lambda and alpha
    
    
    lambda = lambda + step*(sum(P_mat.*gainMatSUTrnsPURecv)-threshold);
    for SUIndex=1:numSUPairs
      vecSUs(SUIndex).alpha = vecSUs(SUIndex).alpha+step*(sum(gainMatSUTrnsSURecv(1:end ~= SUIndex,SUIndex)' ...
                            *P_mat(1:end ~= SUIndex,:).*X_mat(1:end ~= SUIndex,:))+PUIntOnSU(:,SUIndex)'-vecSUs(SUIndex).Gamma);
       vecSUs(SUIndex).lambda = lambda; 
    end
    
    
end
    

end


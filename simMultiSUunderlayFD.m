function [ rateAchievedbySUs,interObservedPU,SNRdB,TotalPowerOfSUs,TotalChannelsOfSUs,rateAchievedbySUsCentral,interObservedPUCentral,SNRdBCentral,TotalPowerOfSUsCentral ...
    ] = simMultiSUunderlayFD( numSUPairs,numPUs,timeSlots,IntThreshold,RxSensitivity,PMax,SUPairPos,PUPos)
%This function simulates an underlay network of decentralized SUs and
%cetralized PUs for a given number of time slots. 

frequency = 2.4e6; %transmit frequency
vecSUs = [];
PMax = 1;

%maxtimum number of iterations
maxiterations = 100;
hs = 0.8; 
alpha = 0.001; 
 %% matrices to  collect power values and channel allocation values
 P_mat = zeros(numSUPairs,numPUs*2);
 X_mat = zeros(numSUPairs,numPUs*2);
 PUIntOnSU = zeros(numPUs,numSUPairs*2);
 %% parameters to collect QoS
 SNRdB = zeros(timeSlots,numPUs);
 interObservedPU = zeros(timeSlots,numPUs);
 rateAchievedbySUs = zeros(timeSlots,numSUPairs);
 TotalPowerOfSUs = zeros(timeSlots,numSUPairs*2);
 TotalChannelsOfSUs = zeros(timeSlots,numSUPairs*2);
 SUIntAtPU = zeros(numPUs,1);
 
 SNRdBCentral = zeros(timeSlots,numPUs);
 interObservedPUCentral = zeros(timeSlots,numPUs);
 rateAchievedbySUsCentral = zeros(timeSlots,numSUPairs);
 TotalPowerOfSUsCentral = zeros(timeSlots,numSUPairs*2);
 TotalChannelsOfSUsCentral = zeros(timeSlots,numSUPairs*2);
 SUIntAtPUCentral = zeros(numPUs,1);
 
    [gainMatSUTrnsSURecv,gainMatSUTrnsPURecv,gainMatPUTrnsSURecv,gainMatSUTransMSRecv,gainMatPUTransPURecv] = chanGainCalc(SUPairPos,PUPos,frequency,hs,alpha);
    
    PUTxPower = 5*(IntThreshold+abs(randn()*sqrt(PMax*1e-6)))./gainMatPUTransPURecv; % PU transmit power
 
    
    for channel=1:numPUs
        PUIntOnSU(channel,:) = gainMatPUTrnsSURecv(channel,:)*PUTxPower(channel,1);
        
    end
 
 %% create SU objects
    
    
 for SUIndex=1:numSUPairs
     %Gamma initialize
     GammaInit = PUIntOnSU(:,[SUIndex,SUIndex+numSUPairs]);
     vecSUs = [vecSUs,SecondaryUserFDJointAlloc(SUIndex,gainMatSUTrnsSURecv([SUIndex,SUIndex+numSUPairs],:),gainMatSUTrnsPURecv([SUIndex,SUIndex+numSUPairs],:),...
                    GammaInit,PMax,maxiterations,hs,alpha,SUIntAtPU,IntThreshold,RxSensitivity,numSUPairs)];
              
                
 end
 
%% simulate for time slots
for timeSlot=1:timeSlots
    disp('Timeslot') 
    disp(timeSlot)
    % PU and SU position generation
     %[ SUPairPos,PUPos ] = positionGenSUandPU( gridSize,numSUPairs,numPUs);
    
     [gainMatSUTrnsSURecv,gainMatSUTrnsPURecv,gainMatPUTrnsSURecv,gainMatSUTransMSRecv,gainMatPUTransPURecv] = chanGainCalc(SUPairPos,PUPos,frequency,hs,alpha);
    if timeSlot >1
        %PUTxPower = 5*(PMax*1e-3)./gainMatPUTransPURecv';
         PUTxPower = 5*(IntThreshold+abs(randn()*sqrt(PMax*1e-6)))./gainMatPUTransPURecv'; % PU transmit power
        %PUTxPower =PUTxPower';
    else
        PUTxPower = 5*(abs(randn()*sqrt(PMax*1e-6))+IntThreshold)./gainMatPUTransPURecv'; % PU transmit power
    end
        
    for channel=1:numPUs
        PUIntOnSU(channel,:) = gainMatPUTrnsSURecv(channel,:)*PUTxPower(1,channel);
        
    end
    IntAtMS=zeros(numPUs,1);
    for SUTxIndex=1:numSUPairs
        IntAtMS=IntAtMS+vecSUs(SUTxIndex).gainVecSUTrnsMSRecv(1,:)'.*vecSUs(SUTxIndex).P(1:numPUs,1);
        IntAtMS=IntAtMS+vecSUs(SUTxIndex).gainVecSUTrnsMSRecv(2,:)'.*vecSUs(SUTxIndex).P(numPUs+1:end,1);
    end
    GammaInit=PUIntOnSU;
    GammaInitCentral=PUIntOnSU;
    for SUIndex=randperm(numSUPairs)
        
        vecSUs(SUIndex).GammaInit = GammaInit(:,[SUIndex,SUIndex+numSUPairs]);
        vecSUs(SUIndex).J = IntAtMS;
        % update channel gains
        vecSUs(SUIndex).gainVecSUTrnsSURecv = gainMatSUTrnsSURecv([SUIndex,SUIndex+numSUPairs],:);
        vecSUs(SUIndex).gainVecSUTrnsMSRecv = gainMatSUTrnsPURecv([SUIndex,SUIndex+numSUPairs],:);
        %initiate X,Gamma and P
        vecSUs(SUIndex).optXPGamma();
        %centralPower = centralizedAllocation(PMax,hs,alpha,gainMatSUTrnsPURecv,2,PUIntOnSU,PUTxPower,maxiterations);
        centralPower = zeros(numPUs,2*numSUPairs);
        for SURxIndex= 1:numSUPairs
            
            GammaInit(:,SURxIndex)=GammaInit(:,SURxIndex)+vecSUs(SUIndex).gainVecSUTrnsSURecv(1,SURxIndex)*vecSUs(SUIndex).P(1:numPUs,1)...
                +vecSUs(SUIndex).gainVecSUTrnsSURecv(2,SURxIndex)*vecSUs(SUIndex).P(numPUs+1:end,1);
            GammaInit(:,SURxIndex+numSUPairs)=GammaInit(:,SURxIndex+numSUPairs)+vecSUs(SUIndex).gainVecSUTrnsSURecv(1,SURxIndex+numSUPairs)*vecSUs(SUIndex).P(1:numPUs,1)...
                +vecSUs(SUIndex).gainVecSUTrnsSURecv(2,SURxIndex+numSUPairs)*vecSUs(SUIndex).P(numPUs+1:end,1);
            
            %% central
            GammaInitCentral(:,SURxIndex)=GammaInitCentral(:,SURxIndex)+vecSUs(SUIndex).gainVecSUTrnsSURecv(1,SURxIndex)*centralPower(1:numPUs,SUIndex)...
                +vecSUs(SUIndex).gainVecSUTrnsSURecv(2,SURxIndex)*centralPower(1:numPUs,SUIndex+numSUPairs);
            GammaInitCentral(:,SURxIndex+numSUPairs)=GammaInitCentral(:,SURxIndex+numSUPairs)+vecSUs(SUIndex).gainVecSUTrnsSURecv(1,SURxIndex+numSUPairs)*centralPower(1:numPUs,SUIndex)...
                +vecSUs(SUIndex).gainVecSUTrnsSURecv(2,SURxIndex+numSUPairs)*centralPower(1:numPUs,SUIndex+numSUPairs);
            
        end
        
                     
    end
    
       
    
    
%     %% execute the optimization of SU parameters in parallel
%     parfor index = 1: numSUPairs
%          vecSUDummy = vecSUs(index);
%          vecSUDummy.optXPGamma();
%          vecSUs(index)=vecSUDummy;
%             %vecSUs(index).optXPGamma();
%     end
    
       
    
    
    %% get the optimized parameters and find if the QoS is met and if the PU interference threshold is violated
    
        for SUIndex=1:numSUPairs
            P_mat(SUIndex,:) = vecSUs(SUIndex).P;
            X_mat(SUIndex,:) = vecSUs(SUIndex).X;
%               disp('underlay');
%               disp(P_mat);
%               disp(X_mat);
%             
            
        end
        n=numPUs;
        for SUIndex=1:numSUPairs
            GammaInitlinear=[GammaInit(:,SUIndex);GammaInit(:,SUIndex+numSUPairs)];
            GammaInitlinearCtrl=[GammaInitCentral(:,SUIndex);GammaInitCentral(:,SUIndex+numSUPairs)];
            cbv = prod([(GammaInitlinear(1:n)+vecSUs(SUIndex).P(n+1:end)*vecSUs(SUIndex).gainVecSUTrnsSURecv(2,SUIndex))./(GammaInitlinear(1:n)); ...
                    (GammaInitlinear(n+1:end)+vecSUs(SUIndex).P(1:n)*vecSUs(SUIndex).gainVecSUTrnsSURecv(1,SUIndex+numSUPairs))./(GammaInitlinear(n+1:end))]);
            
                
            rateAchievedbySUs(timeSlot,SUIndex) =  log2(cbv);
            TotalPowerOfSUs(timeSlot,SUIndex) =  sum(P_mat(SUIndex,:));
            TotalChannelsOfSUs(timeSlot,SUIndex) = sum(X_mat(SUIndex,:));
            
            %% central
            cbvCentral = prod([(GammaInitlinearCtrl(1:n)+centralPower(1:numPUs,SUIndex+numSUPairs)*vecSUs(SUIndex).gainVecSUTrnsSURecv(2,SUIndex))./(GammaInitlinearCtrl(1:n)); ...
                    (GammaInitlinearCtrl(n+1:end)+centralPower(1:numPUs,SUIndex)*vecSUs(SUIndex).gainVecSUTrnsSURecv(1,SUIndex++numSUPairs))./(GammaInitlinearCtrl(n+1:end))]);    
            rateAchievedbySUsCentral(timeSlot,SUIndex) =  log2(cbvCentral);
            TotalPowerOfSUsCentral(timeSlot,SUIndex) =  sum([centralPower(1:numPUs,SUIndex);centralPower(1:numPUs,SUIndex+numSUPairs)]);
            
        end
        
        
        for channel=1:numPUs
            interObservedPU(timeSlot,channel) = gainMatSUTrnsPURecv(1:numSUPairs,channel)'*((P_mat(:,channel)))+gainMatSUTrnsPURecv(numSUPairs+1:end,channel)'*((P_mat(:,channel+numPUs)));
            SNRdB(timeSlot,channel) = 10*log10((PUTxPower(1,channel)*gainMatPUTransPURecv(channel,1))./(abs(randn()*sqrt(PMax*1e-6)) + interObservedPU(timeSlot,channel)));
            %% central
            interObservedPUCentral(timeSlot,channel) = (centralPower(channel,:))*gainMatSUTrnsPURecv(:,channel);
            SNRdBCentral(timeSlot,channel) = 10*log10((PUTxPower(1,channel)*gainMatPUTransPURecv(channel,1))./(abs(randn()*sqrt(PMax*1e-6)) +interObservedPUCentral(timeSlot,channel)));
        end
        
   
    
end
    
   

end


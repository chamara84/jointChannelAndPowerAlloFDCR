%Script to run the optimization program for different settings
clear all
clc
%% Parameters
numSUPairs = 1;
numPUs = 5;
gridSize = 100;
timeSlots = 10000;
phi = 0.1; %for unit band width
K = 1;
PMax = 1;
RxSensitivity = PMax*1e-8;
thresholdVec = 0.01*(1:9)*PMax;
% PU and SU position generation

%[ SUPairPos,PUPos ] = positionGenSUandPU( gridSize,numSUPairs,numPUs);
PUPos = [50,50,1.77679189106659,58.5173969397831;
         50,50,44.9121421483536,92.2081818024207;
         50,50,16.8921728989299,7.3033142145689;
         50,50,59.8929492399754,26.3307248905415;
         50,50,96.0564997301030,50.1707944368526];
%          50,50,47.1793905369834,99.9780722518256;
%          50,50,84.9294800093629,28.9901957945708;
%          50,50,60.0174238027374,95.8361221865986;];

SUPairPos =  [58.3939406627141,81.8854209754088;
              %13.3682967576487,60.8174078727042;
              %63.3448938939210,70.3291107874416; 
              65.3932779042565,70.7468111874852];
              %12.3163300779066,72.9020347185759];
              %56.3537951516341,34.7476745558460];
%               12.5016308214126,32.8200523927005;
%               37.8864931620977,50.4527785781561] ;   
index =1;

percentOfThresholdViolation = zeros(1,length(thresholdVec));
outageProbability=zeros(1,length(thresholdVec));
percentOfPowerConstVio = zeros(1,length(thresholdVec));
percentOfChannelConstVio = zeros(1,length(thresholdVec));
meanInterOnPUs = zeros(timeSlots,length(thresholdVec));
meanRateOfSusIndividual = zeros(numSUPairs,length(thresholdVec));
meanSNROfPUsIndividual = zeros(numPUs,length(thresholdVec));

%% central
percentOfThresholdViolationCtrl = zeros(1,length(thresholdVec));
outageProbabilityCtrl=zeros(1,length(thresholdVec));
meanInterOnPUsCtrl = zeros(timeSlots,length(thresholdVec));
meanRateOfSusIndividualCtrl = zeros(numSUPairs,length(thresholdVec));
meanSNROfPUsIndividualCtrl = zeros(numPUs,length(thresholdVec));
for threshold=thresholdVec
%function to execute optimization
[ rateAchievedbySUs,interObservedPU,SNRdB,TotalPowerOfSUs,TotalChannelsOfSUs,rateAchievedbySUsCentral,interObservedPUCentral,SNRdBCentral,TotalPowerOfSUsCentral ] = simMultiSUunderlayFD( numSUPairs,numPUs,timeSlots,threshold,RxSensitivity,PMax,SUPairPos,PUPos);


%% evaluation and plotting
meanRateOfSusIndividual(:,index) = mean(rateAchievedbySUs,1)';
meanSNROfPUsIndividual(:,index) = mean(SNRdB,1)';

percentOfPowerConstVio(1,index) = length(find(TotalPowerOfSUs>PMax))/(size(TotalPowerOfSUs,1)*size(TotalPowerOfSUs,2));

percentOfThresholdViolation(1,index)  = length(find(interObservedPU>=threshold))/(size(interObservedPU,1)*size(interObservedPU,2));
outageProbability(1,index)  = length(find(SNRdB<10*log10(5)))/(size(SNRdB,1)*size(SNRdB,2));

%%central
meanRateOfSusIndividualCtrl(:,index) = mean(rateAchievedbySUsCentral,1)';
meanSNROfPUsIndividualCtrl(:,index) = mean(SNRdBCentral,1)';
percentOfThresholdViolationCtrl(1,index)  = length(find(interObservedPUCentral>=threshold))/(size(interObservedPUCentral,1)*size(interObservedPUCentral,2));
outageProbabilityCtrl(1,index)  = length(find(SNRdBCentral<10*log10(5)))/(size(SNRdBCentral,1)*size(SNRdBCentral,2));

index=index+1;
end
meanRateOfSus = mean(meanRateOfSusIndividual,1);
meanSNROfPUs = mean(meanSNROfPUsIndividual,1);
%central
meanRateOfSusCtrl =mean(meanRateOfSusIndividualCtrl,1);
meanSNROfPUsCtrl = mean(meanSNROfPUsIndividualCtrl,1);

save(strcat('PUDownlink',num2str(phi),'.mat'),'percentOfThresholdViolation','meanRateOfSus','meanSNROfPUs','percentOfPowerConstVio');

figure(1)
plotNum(1)=plot(thresholdVec,meanSNROfPUs,'-*b');
 hold on
plotNum(2)=plot(thresholdVec,meanSNROfPUsCtrl,'-.xr');

title('mean SNR of PUs Vs. Interference threshold')
xlabel('Threshold')
ylabel('SNR (dB)')
legend (plotNum,'Decentralized','Centralized')
hold off

figure(2)
plotNum(1)=plot(thresholdVec,meanRateOfSus,'-*b');
hold on
plotNum(2)=plot(thresholdVec,meanRateOfSusCtrl,'-.xr');
% plotNum(3)=plot(timeSlots/2:timeSlots/2+timeSlots*0.1,meanRateOfSus(timeSlots/2:timeSlots/2+timeSlots*0.1,6),'--k');
% plotNum(4)=plot(timeSlots/2:timeSlots/2+timeSlots*0.1,meanRateOfSus(timeSlots/2:timeSlots/2+timeSlots*0.1,8),'-.b');

title(' mean data Rate of SUs Vs. Interference threshold');
xlabel('Threshold')
ylabel('Data rate (bps/Hz)')
legend ([plotNum(1),plotNum(2)],'Decentralized','Centralized');

hold off

figure(3)
plotNum(1)=plot(thresholdVec,outageProbability,'-*b');
hold on
plotNum(2)=plot(thresholdVec,outageProbabilityCtrl,'-.xr');
title('Percentage of PU Outage Vs. Threshold');
xlabel('Threshold');
ylabel('Percentage of PU Outage');
legend ([plotNum(1),plotNum(2)],'Decentralized','Centralized');

hold off

figure(4)
plotNum(1)=plot(thresholdVec,percentOfThresholdViolation,'-*b');
hold on
plotNum(2)=plot(thresholdVec,percentOfThresholdViolationCtrl,'-.xr');
title('Percentage of PU threshold violation Vs. Threshold');
xlabel('Threshold');
ylabel('Percentage of PU threshold violation');
legend ([plotNum(1),plotNum(2)],'Decentralized','Centralized');

hold off

figure(5)
plotNum(1)=scatter(PUPos(:,3),PUPos(:,4),'r^');
hold on
plotNum(2)=scatter(PUPos(1,1),PUPos(1,2),'o');
plotNum(3)=scatter(SUPairPos(:,1),SUPairPos(:,2),'k+');
%plotNum(4)=scatter(SUPairPos(:,3),SUPairPos(:,4),'bx');
title('Network Topology');
xlabel('X (m)');
ylabel('Y (m)');
legend ([plotNum(1),plotNum(2),plotNum(3)],'PU Receivers/MSs','PU Base station','CRs');
hold off


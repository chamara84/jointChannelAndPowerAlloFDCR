clear all
clc


PUPos = [50,50,1.77679189106659,58.5173969397831;
         50,50,44.9121421483536,92.2081818024207;
         50,50,16.8921728989299,7.3033142145689;
         50,50,59.8929492399754,26.3307248905415;
         50,50,96.0564997301030,50.1707944368526];
%          50,50,47.1793905369834,99.9780722518256;
%          50,50,84.9294800093629,28.9901957945708;
%          50,50,60.0174238027374,95.8361221865986;];

SUPairPos =  [58.3939406627141,81.8854209754088;
              13.3682967576487,60.8174078727042;
              63.3448938939210,70.3291107874416;
              37.8864931620977,50.4527785781561;
              12.5016308214126,32.8200523927005;
              81.3939406627141,37.8864931620977;
              91.3939406627141,31.8864931620977;
              77.4568123627141,96.4568123627141;
              65.3932779042565,70.7468111874852;
              12.3163300779066,72.9020347185759;
              56.3537951516341,34.7476745558460;
              50.8864931620977,40.4527785781561;
              22.5016308214126,52.8200523927005;
              86.3939406627141,30.8864931620977;
              98.3939406627141,20.8864931620977;
              97.4568123627141,90.4568123627141];

MSPos = [25.0   25.0;
         75.0    25.0;
         25.0  75.0;
         75.0   75.0];

noPred1 = load('./PredResults/NumSUs3React_update1slot.mat');
pred1 = load('./PredResults/NumSUs3Proact_update1slot.mat');
noPred2 = load('./PredResults/NumSUs5React_update1slot.mat');
pred2 = load('./PredResults/NumSUs5Proact_update1slot.mat');
noPred3 = load('./PredResults/NumSUs8React_update1slot.mat');
pred3 = load('./PredResults/NumSUs8Proact_update1Slot.mat');
thresholdVec = 0.001*(1:1:10);

figure(1)
plotNum(1)=plot(thresholdVec,noPred1.meanSNROfPUs,'-*b');
hold on
plotNum(2)=plot(thresholdVec,pred1.meanSNROfPUs,'-.xr');
plotNum(3)=plot(thresholdVec,noPred2.meanSNROfPUs,'-ob');
plotNum(4)=plot(thresholdVec,pred2.meanSNROfPUs,'-.vr');
plotNum(5)=plot(thresholdVec,noPred3.meanSNROfPUs,'-sb');
plotNum(6)=plot(thresholdVec,pred3.meanSNROfPUs,'-.+r');
legend ([plotNum(1),plotNum(2),plotNum(3),plotNum(4),plotNum(5),plotNum(6)],'W\Prediction 3 Pairs','With prediction 3 Pairs',...
    'W\Prediction 5 Pairs','With prediction 5 Pairs ','W\Prediction 8 Pairs','With prediction 8 Pairs ');
title('Mean SNR of PUs Vs. Interference threshold')
xlabel('Threshold')
ylabel('SNR (dB)')
hold off


figure(2)
plotNum(1)=plot(thresholdVec,noPred1.meanRateOfSus,'-*b');
hold on
plotNum(2)=plot(thresholdVec,pred1.meanRateOfSus,'-.xr');
plotNum(3)=plot(thresholdVec,noPred2.meanRateOfSus,'-ob');
plotNum(4)=plot(thresholdVec,pred2.meanRateOfSus,'-.vr');
plotNum(5)=plot(thresholdVec,noPred3.meanRateOfSus,'-sb');
plotNum(6)=plot(thresholdVec,pred3.meanRateOfSus,'-.+r');
legend ([plotNum(1),plotNum(2),plotNum(3),plotNum(4),plotNum(5),plotNum(6)],'W\Prediction 3 Pairs','With prediction 3 Pairs',...
    'W\Prediction 5 Pairs','With prediction 5 Pairs ','W\Prediction 8 Pairs','With prediction 8 Pairs ');% plotNum(3)=plot(timeSlots/2:timeSlots/2+timeSlots*0.1,meanRateOfSus(timeSlots/2:timeSlots/2+timeSlots*0.1,6),'--k');
% plotNum(4)=plot(timeSlots/2:timeSlots/2+timeSlots*0.1,meanRateOfSus(timeSlots/2:timeSlots/2+timeSlots*0.1,8),'-.b');

title(' Mean data Rate of SUs Vs. Interference threshold');
xlabel('Threshold')
ylabel('Data rate (bps/Hz)')

hold off


figure(3)
plotNum(1)=plot(thresholdVec,noPred1.outageProbability,'-*b');
hold on
plotNum(2)=plot(thresholdVec,pred1.outageProbability,'-.xr');
plotNum(3)=plot(thresholdVec,noPred2.outageProbability,'-ob');
plotNum(4)=plot(thresholdVec,pred2.outageProbability,'-.vr');
plotNum(5)=plot(thresholdVec,noPred3.outageProbability,'-sb');
plotNum(6)=plot(thresholdVec,pred3.outageProbability,'-.+r');
legend ([plotNum(1),plotNum(2),plotNum(3),plotNum(4),plotNum(5),plotNum(6)],'W\Prediction 3 Pairs','With prediction 3 Pairs',...
    'W\Prediction 5 Pairs','With prediction 5 Pairs ','W\Prediction 8 Pairs','With prediction 8 Pairs ');
title('Probability of PU Outage Vs. Threshold');
xlabel('Threshold');
ylabel('Probability of PU Outage');
hold off



figure(4)
plotNum(1)=plot(thresholdVec,noPred1.percentOfThresholdViolation,'-*b');
hold on
plotNum(2)=plot(thresholdVec,pred1.percentOfThresholdViolation,'-.xr');
plotNum(3)=plot(thresholdVec,noPred2.percentOfThresholdViolation,'-ob');
plotNum(4)=plot(thresholdVec,pred2.percentOfThresholdViolation,'-.vr');
plotNum(5)=plot(thresholdVec,noPred3.percentOfThresholdViolation,'-sb');
plotNum(6)=plot(thresholdVec,pred3.percentOfThresholdViolation,'-.+r');
legend ([plotNum(1),plotNum(2),plotNum(3),plotNum(4),plotNum(5),plotNum(6)],'W\Prediction 3 Pairs','With prediction 3 Pairs',...
    'W\Prediction 5 Pairs','With prediction 5 Pairs ','W\Prediction 8 Pairs','With prediction 8 Pairs ');
title('Probability of PU threshold violation Vs. Threshold');
xlabel('Threshold');
ylabel('Probability of PU threshold violation');

hold off



figure(5)
plotNum(1)=scatter(PUPos(:,3),PUPos(:,4),'r^');
hold on
plotNum(2)=scatter(PUPos(1,1),PUPos(1,2),'o');
plotNum(3)=scatter(SUPairPos([1:3,9:11],1),SUPairPos([1:3,9:11],2),'k+');
plotNum(4)=scatter(SUPairPos([4:5,12:13],1),SUPairPos([4:5,12:13],2),'k*');
plotNum(5)=scatter(SUPairPos([6:8,14:16],1),SUPairPos([6:8,14:16],2),'ks');
plotNum(6)=scatter(MSPos(:,1),MSPos(:,2),'bx');
title('Network Topology');
xlabel('X (m)');
ylabel('Y (m)');
legend ([plotNum(1),plotNum(2),plotNum(3),plotNum(4),plotNum(5),plotNum(6)],'PU Receivers','PU Base station','CR pairs 1-3','CR pairs 4-5','CR pairs 6-8', 'MSs');
hold off


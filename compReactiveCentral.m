clear all
clc


PUPos = [50,50,21.77679189106659,58.5173969397831;
         50,50,44.9121421483536,30.2081818024207];
%          50,50,16.8921728989299,7.3033142145689;
%          50,50,59.8929492399754,26.3307248905415;
%          50,50,96.0564997301030,50.1707944368526];
%          50,50,47.1793905369834,99.9780722518256;
%          50,50,84.9294800093629,28.9901957945708;
%          50,50,60.0174238027374,95.8361221865986;];

SUPairPos =  [58.3939406627141,81.8854209754088;
              13.3682967576487,60.8174078727042;
              %63.3448938939210,70.3291107874416; 
              %82.5016308214126,32.8200523927005;
              65.3932779042565,70.7468111874852;
              12.3163300779066,72.9020347185759];
              %56.3537951516341,34.7476745558460;
              %97.8864931620977,50.4527785781561] ;   
SUPairPosCentral = [58.3939406627141,81.8854209754088;
                    65.3932779042565,70.7468111874852;
                    13.3682967576487,60.8174078727042;
                    12.3163300779066,72.9020347185759];
              
MSPos = [25.0   25.0;
         75.0    25.0;
         25.0  75.0;
         75.0   75.0];

central = load('./NumSUs2Central.mat');
distributed = load('./NumSUs21000Slots.mat');

thresholdVec = 0.001*(1:1:10);

figure(1)
plotNum(1)=plot(thresholdVec,central.meanSNROfPUs,'-*b');
hold on
%plotNum(2)=plot(thresholdVec,pred1.meanSNROfPUs,'-.xr');
plotNum(2)=plot(thresholdVec,distributed.meanSNROfPUs,'-or');

legend ([plotNum(1),plotNum(2)],'Central',...
    'Distributed');
title('Mean SINR of PUs Vs. Interference threshold')
xlabel('Threshold')
ylabel('SINR (dB)')
hold off


figure(2)
plotNum(1)=plot(thresholdVec,central.meanRateOfSus,'-*b');
hold on
%plotNum(2)=plot(thresholdVec,pred1.meanRateOfSus,'-.xr');
plotNum(2)=plot(thresholdVec,distributed.meanRateOfSus./2,'-or');

legend ([plotNum(1),plotNum(2)],'Central',...
    'Distributed');


title(' Mean data Rate of SUs Vs. Interference threshold');
xlabel('Threshold')
ylabel('Data rate (bps/Hz)')

hold off


figure(3)
plotNum(1)=plot(thresholdVec,central.outageProbability,'-*b');
hold on
%plotNum(2)=plot(thresholdVec,pred1.outageProbability,'-.xr');
plotNum(2)=plot(thresholdVec,distributed.outageProbability,'-or');
legend ([plotNum(1),plotNum(2)],'Central',...
    'Distributed');
title('Probability of PU Outage Vs. Threshold');
xlabel('Threshold');
ylabel('Probability of PU Outage');
hold off



figure(4)
plotNum(1)=plot(thresholdVec,central.percentOfThresholdViolation,'-*b');
hold on
%plotNum(2)=plot(thresholdVec,pred1.percentOfThresholdViolation,'-.xr');
plotNum(2)=plot(thresholdVec,distributed.percentOfThresholdViolation,'-or');
legend ([plotNum(1),plotNum(2)],'Central',...
    'Distributed');
title('Probability of PU threshold violation Vs. Threshold');
xlabel('Threshold');
ylabel('Probability of PU threshold violation');

hold off



figure(5)
plotNum(1)=scatter(PUPos(:,3),PUPos(:,4),'r^');
hold on
plotNum(2)=scatter(PUPos(1,1),PUPos(1,2),'o');
plotNum(3)=scatter(SUPairPos([1,3],1),SUPairPos([1,3],2),'k+');
plotNum(4)=scatter(SUPairPos([2,4],1),SUPairPos([2,4],2),'k*');
plotNum(6)=scatter(MSPos(:,1),MSPos(:,2),'bx');
title('Network Topology');
xlabel('X (m)');
ylabel('Y (m)');
legend ([plotNum(1),plotNum(2),plotNum(3),plotNum(4),plotNum(6)],'PU Receivers','PU Base station','CR pair 1','CR pairs 2', 'MSs');
hold off


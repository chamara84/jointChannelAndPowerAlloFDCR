%Script to run the optimization program for different settings

%% Parameters
numSUPairs = 4;
numPUs = 5;
gridSize = 200;
timeSlots = 1000;
phi = 0.5; %for unit band width
K = 1;
PMax = 1;

thresholdVec = 0.01*(1:10)*PMax;
% PU and SU position generation
[ SUPairPos,PUPos ] = positionGenSUandPU( gridSize,numSUPairs,numPUs);
index =1;
percentOfQoSSatSU = zeros(1,length(thresholdVec));
percentOfThresholdViolation = zeros(1,length(thresholdVec));

for threshold=thresholdVec
%function to execute optimization
[ rateAchievedbySUs,interferenceOnPU ] = ...
    simMultiSUunderlayFD( numSUPairs,numPUs,gridSize,timeSlots,rxSensitivity,threshold,PMax,SUPairPos,PUPos);


%% evaluation and plotting

percentOfQoSSatSU(1,index) = length(find(rateAchievedbySUs>phi))/(size(rateAchievedbySUs,1)*size(rateAchievedbySUs,2));

percentOfThresholdViolation(1,index)  = length(find(interferenceOnPU>threshold))/(size(interferenceOnPU,1)*size(interferenceOnPU,2));

index=index+1;
end
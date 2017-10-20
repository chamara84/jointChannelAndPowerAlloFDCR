function [ SUPairPos,PUPos ] = positionGenSUandPU( gridSize,numSUPairs,numPUs)
%generates the positions of PUs and SU pairs according to a uniform
%distribution on square grid of size L. The PU network is considered to be
%centralized and downlink is considered




%SUPairPos = [XposTransmitter,YposTransmitter,XposReceiver,YposReceiver]
%PUPos = [XPosPUTx,YPosPUTx,XPosPURx,YPosPURx]

SUPairPos = zeros(numSUPairs*2,2);
PUPos = zeros(numPUs,4);

for SUpoints=1:numSUPairs
    
    if mod(SUpoints/2,2)==0 && mod(SUpoints,4)~=0
        SUPairPos(SUpoints,1)=0.5*gridSize+0.5*gridSize*rand();
        SUPairPos(SUpoints,2)=0.5*gridSize*rand();
        SUPairPos(SUpoints+numSUPairs,1)=0.5*gridSize+0.5*gridSize*rand();
        SUPairPos(SUpoints+numSUPairs,2)=0.5*gridSize*rand();
    elseif mod(SUpoints/2,3)==0
        SUPairPos(SUpoints,1)=0.5*gridSize*rand();
        SUPairPos(SUpoints,2)=0.5*gridSize*rand();
        SUPairPos(SUpoints+numSUPairs,1)=0.5*gridSize*rand();
        SUPairPos(SUpoints+numSUPairs,2)=0.5*gridSize*rand();
    elseif mod(SUpoints/2,4)==0
        SUPairPos(SUpoints,1)=0.5*gridSize*rand();
        SUPairPos(SUpoints,2)=0.5*gridSize+0.5*gridSize*rand();
        SUPairPos(SUpoints+numSUPairs,1)=0.5*gridSize*rand();
        SUPairPos(SUpoints+numSUPairs,2)=0.5*gridSize+0.5*gridSize*rand();
    else 
        SUPairPos(SUpoints,1)=0.5*gridSize+0.5*gridSize*rand();
        SUPairPos(SUpoints,2)=0.5*gridSize+0.5*gridSize*rand();
        SUPairPos(SUpoints+numSUPairs,1)=0.5*gridSize+0.5*gridSize*rand();
        SUPairPos(SUpoints+numSUPairs,2)=0.5*gridSize+0.5*gridSize*rand();
    end
        
end

%% monitoring station position
MonitoringStationx = gridSize/2;
MonitoringStationy = gridSize/2;

%setting PU tx and rx positions
for PUpoints = 1:numPUs
    PUPos(PUpoints,1)= MonitoringStationx;
    PUPos(PUpoints,2)=MonitoringStationy;
    PUPos(PUpoints,3)= gridSize*rand();
    PUPos(PUpoints,4)=gridSize*rand();
end

function [gainMatSUTrnsSURecv,gainMatSUTrnsPURecv,gainMatPUTrnsSURecv,gainMatSUTransMSRecv,gainMatPUTransPURecv] = chanGainCalc(SUPairPos,PUPos,frequency,hs,alpha)
%% channel gain calculation between SU transmitter and SU receiver in gainMatSUTrnsSURecv
%SURecv       1  2  3 ... 2*M        
%SUTrans   1
%          2
%          :
%          2*M
%considering PU system to be centralized and the downlink duration
%assume freespace loss at 10m
gainMatSUTrnsSURecv = zeros(size(SUPairPos,1),size(SUPairPos,1));
gainMatSUTrnsPURecv = zeros(size(SUPairPos,1),size(PUPos,1));
gainMatPUTrnsSURecv = zeros(size(PUPos,1),size(SUPairPos,1));
gainMatSUTransMSRecv = zeros(size(SUPairPos,1),1);
gainMatPUTransPURecv = zeros(size(PUPos,1),1);
n= 2.7;
sigma_fading = 11.8;
FSPL = 20*log10(10)+20*log10(frequency)-147.55;

for SUTransIndex=1:size(SUPairPos,1)
    
    for  SURecvIndex=1:size(SUPairPos,1)
        
        if SUTransIndex~=SURecvIndex
            d = sqrt((SUPairPos(SUTransIndex,1)-SUPairPos(SURecvIndex,1))^2+(SUPairPos(SUTransIndex,2)-SUPairPos(SURecvIndex,2))^2);
            PL = FSPL + 10*n*log10(d/10)+sigma_fading.*randn(1,1);
            gainMatSUTrnsSURecv(SUTransIndex,SURecvIndex)= abs(sqrt(10^(-PL/10)/2)*(randn(1,1)+1i*rand(1,1)));
        else
            gainMatSUTrnsSURecv(SUTransIndex,SURecvIndex)=hs*alpha;
            
        end
    end
end




%% channel gain calculation between SU transmitter and PU receiver in gainMatSUTrnsPURecv
%PURecv       1  2  3 ... N        
%SUTrans   1
%          2
%          :
%          2*M

for SUTransIndex=1:size(SUPairPos,1)
   for  PURecvIndex=1:size(PUPos,1)
       d = sqrt((SUPairPos(SUTransIndex,1)-PUPos(PURecvIndex,3))^2+(SUPairPos(SUTransIndex,2)-PUPos(PURecvIndex,4))^2);
       PL = FSPL + 10*n*log10(d/10)+sigma_fading.*randn(1,1);
       gainMatSUTrnsPURecv(SUTransIndex,PURecvIndex)= abs(sqrt(10^(-PL/10)/2)*(randn(1,1)+1i*rand(1,1)));
       
   end
end
%% channel gain calculation between PU transmitter and SU receiver in gainMatSURecvPUTrns
%SURecv       1  2  3 ... 2*M        
%PUTrans   1
%          2
%          :
%          N

for PUTransIndex=1:size(PUPos,1)
   for  SURecvIndex=1:size(SUPairPos,1)
       d = sqrt((PUPos(PUTransIndex,1)-SUPairPos(SURecvIndex,1))^2+(PUPos(PUTransIndex,2)-SUPairPos(SURecvIndex,2))^2);
       PL = FSPL + 10*n*log10(d/10)+sigma_fading.*randn(1,1);
       gainMatPUTrnsSURecv(PUTransIndex,SURecvIndex)= abs(sqrt(10^(-PL/10)/2)*(randn(1,1)+1i*rand(1,1)));
       
   end
end




%% channel gain calculation between SU transmitter and the monitoring station
%Monitoring station 1
% SUTx  1
%       2
%       3
%       :
%       2*M

for MSIndex=1:1
   for  SUTxIndex=1:size(SUPairPos,1)
       d = sqrt((PUPos(MSIndex,1)-SUPairPos( SUTxIndex,1))^2+(PUPos(MSIndex,2)-SUPairPos( SUTxIndex,2))^2);
       PL = FSPL + 10*n*log10(d/10)+sigma_fading.*randn(1,1);
       gainMatSUTransMSRecv( SUTxIndex,MSIndex)= abs(sqrt(10^(-PL/10)/2)*(randn(1,1)+1i*rand(1,1)));
       
   end
end
%% channel gain calculation between PU transmitter and PU receiver

%PUreceiver         1
%                   2
%                   :
%                   N
for PUTxIndex=1:size(PUPos,1)
   
       d = sqrt((PUPos(PUTxIndex,1)-PUPos( PUTxIndex,3))^2+(PUPos(PUTxIndex,2)-PUPos( PUTxIndex,2))^2);
       PL = FSPL + 10*n*log10(d/10)+sigma_fading.*randn(1,1);
       gainMatPUTransPURecv( PUTxIndex,1)= abs(sqrt(10^(-PL/10)/2)*(randn(1,1)+1i*rand(1,1)));
       
   
end
end


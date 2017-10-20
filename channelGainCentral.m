function [h,g] = channelGainCentral(SUPairPos,MSPos,frequency,hs,alpha)
%h(zeta,theta,j,k) = channel gain between user zeta of pair j and user theta of pair k%
%g(theta,k,m) = channel gain between user theta of pair k and monitoring
%station m

%(theta,k) state name
%SUPairPos  =             x     y
%              (1,1)   
%              (2,1)
%              (1,2)
%              (2,2) 
n= 2.7;
sigma_fading = 11.8;
FSPL = 20*log10(10)+20*log10(frequency)-147.55;
h = zeros(2,2,2,2);
for j=1:2
    for k=1:2
        for zeta=1:2
            for theta=1:2
                if j==k && zeta~=theta
                    d = sqrt((SUPairPos(theta+(k-1)*2,1)-SUPairPos(zeta+(j-1)*2,1))^2+(SUPairPos(theta+(k-1)*2,2)-SUPairPos(zeta+(j-1)*2,2))^2);
                    PL = FSPL + 10*n*log10(d/10)+sigma_fading.*randn(1,1);
                    h(theta,zeta,j,k) = abs(sqrt(10^(-PL/10)/2)*(randn(1,1)+1i*randn(1,1)));
                elseif j==k && zeta==theta
                    h(theta,zeta,j,k)=hs*alpha;
                elseif j~=k 
                    d = sqrt((SUPairPos(theta+(k-1)*2,1)-SUPairPos(zeta+(j-1)*2,1))^2+(SUPairPos(theta+(k-1)*2,2)-SUPairPos(zeta+(j-1)*2,2))^2);
                    PL = FSPL + 10*n*log10(d/10)+sigma_fading.*randn(1,1);
                    h(theta,zeta,j,k) = abs(sqrt(10^(-PL/10)/2)*(randn(1,1)+1i*randn(1,1)));
                end
            end
        end
    end
end
%% channel gain calculation between SU transmitter and the monitoring station
%Monitoring station 1 2 3 4
% SUTx  (1,1)
%       (2,1)
%       (1,2)
%       (2,2)



g = zeros(2,2,2);

for m=1:4
    for k=1:2
        for theta = 1:2
            d = sqrt((MSPos(m,1)-SUPairPos(theta+(k-1)*2,1))^2+(MSPos(m,2)-SUPairPos(theta+(k-1)*2,2))^2);
            PL = FSPL + 10*n*log10(d/10)+sigma_fading.*randn(1,1);
            g(theta,k,m)= abs(sqrt(10^(-PL/10)/2)*(randn(1,1)+1i*randn(1,1)));
        end
    end
end






                 
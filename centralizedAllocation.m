function [SUPower] = centralizedAllocation(Pmax,hs,alpha,gainMatSUTrnsPURecv,minSNR,PUInterference,PUPower,MaxIts)
           numPUs = size(gainMatSUTrnsPURecv,2);
           numSUs = size(gainMatSUTrnsPURecv,1);
           noise = repmat(randn(1,numSUs)*Pmax*1e-6,[numPUs,1]);     
          %noise matrix: %1,2,3,...,2*numSUPairs
                         % 1
                         % :
                         % numPUs
           
           c_i = noise+PUInterference;
           C = minSNR./PUPower';
           lambda = zeros(numPUs);
           lambdaOld = ones(numPUs);
           SUPower = zeros(numPUs,numSUs);
           zeta=1.05;
           for channelNum = 1:numPUs
               for iteration=1:MaxIts
                   for SUIndex=1:numSUs
                       SUPower(channelNum,SUIndex)= max((-c_i(channelNum,SUIndex)*C(channelNum,1)*(lambdaOld(channelNum)-1)...
                           +sqrt(c_i(channelNum,SUIndex)^2*C(channelNum,1)^2*(lambdaOld(channelNum)-1)^2+4*c_i(channelNum,SUIndex)*C(channelNum,1)...
                           *lambdaOld(channelNum)*alpha*hs))/(2*C(channelNum,1)*lambdaOld(channelNum)*hs*alpha),0);
                   end
                   
                   
                   
                   lambda(channelNum)= max(prod(1+C(channelNum,1)* SUPower(channelNum,:))-zeta+lambdaOld(channelNum),0);
                   if abs(lambdaOld(channelNum)- lambda(channelNum))<0.0001
                       disp(SUPower);
                       break;
                   else
                       lambdaOld(channelNum) = lambda(channelNum);
                   end
               end
           end
          
                  
                  
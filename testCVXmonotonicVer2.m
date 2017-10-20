
clear all
%clc

numChannels =2;
n=numChannels; 
SUIndex = 1;
gainVecSUTrnsSURecv = 0.6;
gainVecSUTrnsPURecv = [ 0.2947;0.2307];%rand(2,1)*0.3;
Xinit = ones(numChannels*2,1);
epsilon=0.1;
GammaInit = [0.0581 ,   0.0928
             0.0580 ,   0.0017];
%              0.0121 ,   0.0863;
%              0.0484 ,   0.0845;
%              0.0209 ,   0.0552];%rand(2,numChannels)'*0.1;
Pmax = 1;
Pinit = ones(numChannels,1)/numChannels*Pmax;
maxiterations = 100;
hs = 0.8;
alpha = 0.001;
SUIntAtPU = [ 0.0504;0.0026];%0.0492;0.0290;0.0040];
IntThreshold = ones(numChannels,1)*0.1;
RxSensitivity  = 0.001;
oldSol =  0;
newSol =1e-8 ;


coeff = [gainVecSUTrnsSURecv, -RxSensitivity*hs*alpha];
coeffMat = zeros(n*2);
rhsMat = zeros(n*2,1);

for i=1:2*n
    
    if mod(i,2)~= 0
        coeffMat(i,i:i+1)= coeff;
        rhsMat(i) = GammaInit(floor(i/2)+1,1);
    else
        coeffMat(i,i-1:i)= [-RxSensitivity*hs*alpha, gainVecSUTrnsSURecv];
        rhsMat(i) = GammaInit(i/2,2);
    end
end

pfeasible = coeffMat\rhsMat;
totalPower = zeros(2,n);
powerUpperbound = zeros(2*n,1);
for i=1:n
    totalPower(2,i)=pfeasible((i-1)*2+1); % group B
    totalPower(1,i)=pfeasible((i)*2);% group A
end





if sum(totalPower(1,:))<=Pmax && sum(totalPower(2,:))<=Pmax
    isFeasible = 1;
else
    isFeasible = 0;
end

intAtPU=totalPower(1,:)*gainVecSUTrnsPURecv(1,1)+totalPower(2,:)*gainVecSUTrnsPURecv(2,1);

if isFeasible && all(intAtPU<=(IntThreshold'-SUIntAtPU'))
    isFeasible = 1;
else
    isFeasible = 0;
end


A = zeros(3*n+4,2*n);
AUtopia = zeros(2*n,1);
b = zeros(3*n+4,1);
bUtopia = zeros(2*n,1);
xNew = zeros(2*n,1);
x_k=zeros(2*n,1);
Ired = eye(n);
%% integer linear program
for iterations=1:100
    
    f = zeros(1,2*n);
    for i=1:n
        for theta=1:2
            if theta==2
                f(1,(i-1)*2+1)=Xinit(i*2,1)*(GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(theta,i)*hs*alpha)...
                    /(GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(theta,i)*hs*alpha+totalPower(1,i)*gainVecSUTrnsSURecv)...
                    *totalPower(1,i)*gainVecSUTrnsSURecv/((GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(theta,i)*hs*alpha)^2) ...
                    *totalPower(theta,i)*hs*alpha + log2(1+totalPower(theta,i)*gainVecSUTrnsSURecv/(GammaInit(i,1)+Xinit(i*2,1)*totalPower(1,i)*hs*alpha));
                A(4+(i-1)*2+1,(i-1)*2+1)= -totalPower(theta,i)*gainVecSUTrnsSURecv/(GammaInit(i,1)+totalPower(1,i)*hs*alpha)+RxSensitivity;
                b(4+(i-1)*2+1,1) = 0;
                
                AUtopia((i-1)*2+1,1) = gainVecSUTrnsPURecv(2,1);
                bUtopia((i-1)*2+1,1)= max(IntThreshold(i,1)-SUIntAtPU(i,1),0);
            else
                f(1,i*2)=Xinit((i-1)*2+1,1)*(GammaInit(i,1)+Xinit(i*2,1)*totalPower(theta,i)*hs*alpha)...
                    /(GammaInit(i,1)+Xinit(i*2,1)*totalPower(theta,i)*hs*alpha+totalPower(2,i)*gainVecSUTrnsSURecv)...
                    *totalPower(2,i)*gainVecSUTrnsSURecv/((GammaInit(i,1)+Xinit(i*2,1)*totalPower(theta,i)*hs*alpha)^2)...
                    *totalPower(theta,i)*hs*alpha + log2(1+totalPower(theta,i)*gainVecSUTrnsSURecv/(GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(2,i)*hs*alpha));
                A(4+i*2,i*2)= -totalPower(theta,i)*gainVecSUTrnsSURecv/(GammaInit(i,2)+totalPower(2,i)*hs*alpha)+RxSensitivity;
                b(4+i*2,1) = 0;
                AUtopia(i*2,1) = gainVecSUTrnsPURecv(1,1);
                bUtopia(i*2,1)= max(IntThreshold(i,1)-SUIntAtPU(i,1),0);
            end
        end
        A(2*n+4+i,(i-1)*2+1:(i-1)*2+2) = [totalPower(2,i)*gainVecSUTrnsPURecv(2,1),totalPower(1,i)*gainVecSUTrnsPURecv(1,1)];
        b(2*n+4+i,1)=IntThreshold(i,1)-SUIntAtPU(i,1);
    end
    
    A(1,:)=-kron(ones(1,n),[1,0]);
    A(2,:)=-kron(ones(1,n),[0,1]);
    b(1:2,1) = 0;
    A(3,:)=kron(totalPower(1,:),[0,1]);
    A(4,:)=kron(totalPower(2,:),[1,0]);
    b(3:4,1) = Pmax;
    intcon = 1:2*n;
    lb=0;
    ub=1;
    
    xNew = bintprog(-f,A,b);
    if norm(Xinit-xNew)<0.1
        break;
    end
    Xinit=xNew;
end
%%lower and the upper bound of power

boxesupperbound = bUtopia./AUtopia.*xNew;
GammaInitlinear=[GammaInit(:,1);GammaInit(:,2)];
powerFeasible = zeros(2*n,1);
for i=1:n
    powerUpperbound(n+i,1)=boxesupperbound((i-1)*2+1); % group B
    powerUpperbound(i,1)=boxesupperbound((i)*2);% group A
    powerFeasible(n+i,1)=pfeasible((i-1)*2+1); % group B
    powerFeasible(i,1)=pfeasible((i)*2);% group A
    x_k(n+i,1)=xNew((i-1)*2+1); % group B
    x_k(i,1)=xNew(i*2); % group A
end

if isFeasible
    
   
    cub = prod([(GammaInitlinear(1:n)+powerUpperbound(1:n)*hs*alpha+powerUpperbound(n+1:end)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+powerUpperbound(1:n)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*hs*alpha+powerUpperbound(1:n)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*hs*alpha)]);
    cbv = prod([(GammaInitlinear(1:n)+powerFeasible(1:n)*hs*alpha+powerFeasible(n+1:end)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+powerFeasible(1:n)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+powerFeasible(n+1:end)*hs*alpha+powerFeasible(1:n)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+powerFeasible(n+1:end)*hs*alpha)]);
else
    pfeasible = pfeasible.*xNew;
    
    for i=1:n
       
        powerFeasible(n+i,1)=pfeasible((i-1)*2+1); % group B
        powerFeasible(i,1)=pfeasible((i)*2);% group A
    end

    %cub = prod([(GammaInitlinear(1:n)+powerUpperbound(1:n)*hs*alpha+powerUpperbound(n+1:end)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+powerUpperbound(1:n)*hs*alpha); ...
     %   (GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*hs*alpha+powerUpperbound(1:n)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*hs*alpha)]);
    cbv = prod([(GammaInitlinear(1:n)+powerFeasible(1:n)*hs*alpha+powerFeasible(n+1:end)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+powerFeasible(1:n)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+powerFeasible(n+1:end)*hs*alpha+powerFeasible(1:n)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+powerFeasible(n+1:end)*hs*alpha)]);
end            
cub = prod([(GammaInitlinear(1:n)+powerUpperbound(1:n)*hs*alpha+powerUpperbound(n+1:end)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+powerUpperbound(1:n)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*hs*alpha+powerUpperbound(1:n)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*hs*alpha)]);;
boxeslowerbound = powerFeasible;
boxesupperbound = min(powerUpperbound,Pmax);
p_k = zeros(2*n,maxiterations);

%%polyblock algorithm
initPowerFeasible=powerFeasible;
for k = 1:maxiterations
    %reduce the current boxes
    disp('iteration=')
    disp(k)
    eraseBoxes=[];
    for boxNum=1:size(boxesupperbound,2)
        intAtPU=boxeslowerbound(1:n,boxNum)*gainVecSUTrnsPURecv(1,1)+boxeslowerbound(n+1:end,boxNum)*gainVecSUTrnsPURecv(2,1);
        objFuncUpperbound= prod([(GammaInitlinear(1:n)+boxesupperbound(1:n,boxNum)*hs*alpha+boxesupperbound(n+1:end,boxNum)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+boxesupperbound(1:n,boxNum)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+boxesupperbound(n+1:end,boxNum)*hs*alpha+boxesupperbound(1:n,boxNum)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+boxesupperbound(n+1:end,boxNum)*hs*alpha)]);
        if sum(boxeslowerbound(1:n,boxNum))<=Pmax && sum(boxeslowerbound(n+1:end,boxNum))<=Pmax && all(intAtPU<=(IntThreshold-SUIntAtPU))
            if all(boxesupperbound(n+1:end,boxNum)*gainVecSUTrnsSURecv./(GammaInitlinear(1:n)+boxesupperbound(1:n,boxNum)*hs*alpha)>=RxSensitivity.*x_k(n+1:end,1)) && ...
                    all(boxesupperbound(1:n,boxNum)*gainVecSUTrnsSURecv./(GammaInitlinear(n+1:end)+boxesupperbound(n+1:end,boxNum)*hs*alpha)>=RxSensitivity.*x_k(1:n,1)) && objFuncUpperbound>cbv
                %do the reduction
                alphaVec=[];
                for individual=1:2
                    for row=1:n
                        if individual==1 && x_k(row,1)>0
                            cvx_begin quiet
                            variable alphaOpt nonnegative
                            maximize alphaOpt
                            subject to
                            sum(boxeslowerbound(1:n,boxNum)+alphaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row))<=Pmax;
                            %sum(boxeslowerbound(n+1:end,boxNum)+alphaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))<=Pmax;
                            (boxeslowerbound(row,boxNum)+boxeslowerbound(n+row,boxNum)+alphaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row))*max(gainVecSUTrnsPURecv(1,:))<=IntThreshold(row,1)-SUIntAtPU(row,1);
                            %(boxeslowerbound(n+1:end,boxNum)+alphaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*max(gainVecSUTrnsPURecv(2,:))<=IntThreshold-SUIntAtPU;
                            alphaOpt <= 1;
                            cvx_end
                            if ~isnan(alphaOpt)
                                alphaVec=[alphaVec,alphaOpt];
                            else
                                alphaVec=[alphaVec,1];
                            end
                            
                            boxesupperbound(row,boxNum)=boxeslowerbound(row,boxNum)+alphaVec(row)*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum));
                        elseif individual==1 && x_k(row,1)== 0
                            alphaVec=[alphaVec,1];
                            boxesupperbound(row,boxNum)=boxeslowerbound(row,boxNum)+alphaVec(row)*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum));
                        elseif individual==2 && x_k(n+row,1)>0
                            cvx_begin quiet
                            variable alphaOpt nonnegative
                            maximize alphaOpt
                            subject to
                            %sum(boxeslowerbound(1:n,boxNum)+alphaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))<=Pmax;
                            sum(boxeslowerbound(n+1:end,boxNum)+alphaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))*Ired(:,row))<=Pmax;
                            %(boxeslowerbound(1:n,boxNum)+alphaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))*max(gainVecSUTrnsPURecv(1,:))<=IntThreshold-SUIntAtPU;
                            (boxeslowerbound(row,boxNum)+boxeslowerbound(n+row,boxNum)+alphaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))*Ired(:,row))*max(gainVecSUTrnsPURecv(2,:))<=IntThreshold(row,1)-SUIntAtPU(row,1);
                            alphaOpt <= 1;
                            cvx_end
                            
                            if ~isnan(alphaOpt)
                                alphaVec=[alphaVec,alphaOpt];
                            else
                                alphaVec=[alphaVec,1];
                            end
                            boxesupperbound(n+row,boxNum)=boxeslowerbound(n+row,boxNum)+alphaVec(n+row)*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum));
                        elseif individual==2 && x_k(n+row,1)== 0
                            alphaVec=[alphaVec,1];    
                            boxesupperbound(n+row,boxNum)=boxeslowerbound(n+row,boxNum)+alphaVec(n+row)*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum));
                        end
                        
                    end
                end
                
                
                
                betaVec=[];
                for individual=1:2
                       for row=1:n
                            if individual==1 && x_k(row,1)>0
                                cvx_begin quiet
                                variable betaOpt nonnegative
                                maximize betaOpt
                                subject to
                                sum(boxesupperbound(1:n,boxNum)-betaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row))>=Pmax;
                                %sum(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))>=Pmax;
                                (boxesupperbound(row,boxNum)+boxesupperbound(n+row,boxNum))*max(gainVecSUTrnsPURecv(1,:))-betaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row)*max(gainVecSUTrnsPURecv(1,:))>=IntThreshold(row,1)-SUIntAtPU(row,1);
                                %(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*max(gainVecSUTrnsPURecv(2,:))>=IntThreshold-SUIntAtPU;
                                (boxesupperbound(row,boxNum)-betaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row))*gainVecSUTrnsSURecv./(GammaInitlinear(n+row)+(boxesupperbound(n+row,boxNum))*hs*alpha)>=RxSensitivity;
                                %(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*gainVecSUTrnsSURecv-(GammaInitlinear(1:n)+(boxesupperbound(1:n,boxNum)-betaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))*hs*alpha)*RxSensitivity>=0;
                                betaOpt <= 1;
                                
                                cvx_end
                                if ~isnan(betaOpt)
                                    betaVec=[betaVec,betaOpt];
                                else
                                    betaVec=[betaVec,1];
                                end
                                
                                boxeslowerbound(row,boxNum)=boxesupperbound(row,boxNum)-betaVec(row)*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum));     
                                
                            elseif individual==1 && x_k(row,1)==0
                                betaVec=[betaVec,1];
                                boxeslowerbound(row,boxNum)=boxesupperbound(row,boxNum)-betaVec(row)*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum));
                            elseif individual==2 && x_k(n+row,1)>0
                               
                                    cvx_begin quiet
                                    variable betaOpt nonnegative
                                    maximize betaOpt
                                    subject to
                                    %sum(boxesupperbound(1:n,boxNum)-betaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))>=Pmax;
                                    sum(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))*Ired(:,row))>=Pmax;
                                    (boxesupperbound(row,boxNum)+boxesupperbound(n+row,boxNum)-betaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))'*Ired(:,row))*max(gainVecSUTrnsPURecv(2,:))>=IntThreshold(row,1)-SUIntAtPU(row,1);
                                    %(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*max(gainVecSUTrnsPURecv(2,:))>=IntThreshold-SUIntAtPU;
                                    (boxesupperbound(n+row,boxNum)-betaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))'*Ired(:,row))*gainVecSUTrnsSURecv./(GammaInitlinear(row)+(boxesupperbound(row,boxNum))*hs*alpha)>=RxSensitivity;
                                    %(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*gainVecSUTrnsSURecv-(GammaInitlinear(1:n)+(boxesupperbound(1:n,boxNum)-betaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))*hs*alpha)*RxSensitivity>=0;
                                    betaOpt <= 1;
                                    
                                    cvx_end
                                    
                                    if ~isnan(betaOpt)
                                        betaVec=[betaVec,betaOpt];
                                    else
                                        betaVec=[betaVec,1];
                                    end
                                    boxeslowerbound(n+row,boxNum)=boxesupperbound(n+row,boxNum)-betaVec(n+row)*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum));
                            elseif individual==2 && x_k(n+row,1)==0
                                betaVec=[betaVec,1];
                                boxeslowerbound(n+row,boxNum)=boxesupperbound(n+row,boxNum)-betaVec(n+row)*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum));
                            end
                        end
                end
                
                
            else
                eraseBoxes=[eraseBoxes,boxNum];
                
            end
                
        
            
        end
        
        
        
       
    end
    
    boxeslowerbound(:,eraseBoxes)=[];
    boxesupperbound(:,eraseBoxes)=[];
    cub(eraseBoxes)=[];
    if isempty(boxesupperbound)
        break
    end
    powerUpperbound=boxesupperbound;
    
    %find the maximum feasible
    
    %a_prime = a-norm(boxesupperbound-boxeslowerbound)/4*ones(n*2,1);
    
    
    [maxCUB,maxUpperBoundIndex] =max(cub); 
    
    
    z = powerUpperbound(:,maxUpperBoundIndex);
    a=boxeslowerbound(:,maxUpperBoundIndex);
    a_prime = a-norm(boxesupperbound-boxeslowerbound)/4*ones(n*2,1);
    
%     disp(maxUpperBoundIndex);
%     disp(z);
    %check if z is feasible if so stop
    intAtPU=z(1:n,1)*gainVecSUTrnsPURecv(1,1)+z(n+1:end,1)*gainVecSUTrnsPURecv(2,1);
        objFuncUpperbound= prod([(GammaInitlinear(1:n)+z(1:n,1)*hs*alpha+z(n+1:end,1)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+z(1:n,1)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+z(n+1:end,1)*hs*alpha+z(1:n,1)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+z(n+1:end,1)*hs*alpha)]);
        if sum(z(1:n,1))<=Pmax && sum(z(n+1:end,1))<=Pmax && all(intAtPU<=(IntThreshold-SUIntAtPU))
            if all(z(n+1:end,1)*gainVecSUTrnsSURecv./(GammaInitlinear(1:n)+z(1:n,1)*hs*alpha)>=RxSensitivity.*x_k(n+1:end,1)) && ...
                    all(z(1:n,1)*gainVecSUTrnsSURecv./(GammaInitlinear(n+1:end)+z(n+1:end,1)*hs*alpha)>=RxSensitivity.*x_k(1:n,1)) && objFuncUpperbound>cbv
                powerFeasible=z;
                break;
                
            end
        end
    
    
    cvx_begin quiet
    variable lambda 
    minimize lambda
    subject to
    %z-lambda*(z-a_prime)-powerFeasible<=0
    sum(z(1:n)-lambda*(z(1:n)-a(1:n)))<=Pmax;
    sum(z(n+1:end)-lambda*(z(n+1:end)-a(n+1:end)))<=Pmax;
    (z(1:n)-lambda*(z(1:n)-a(1:n)))*max(gainVecSUTrnsPURecv(1,:))+ (z(n+1:end)-lambda*(z(n+1:end)-a(n+1:end)))*max(gainVecSUTrnsPURecv(2,:))<=IntThreshold-SUIntAtPU;
%     
    %(z_vec(1,:)-lambda*(z_vec(1,:)-a_prime_vec(1,:)))'>=zeros(n,1)
    %(z_vec(2,:)-lambda*(z_vec(2,:)-a_prime_vec(2,:)))'>=zeros(n,1)
    cvx_end
    p_k(:,k) = [(z(1:n)-lambda*(z(1:n)-a(1:n)));(z(n+1:end)-lambda*(z(n+1:end)-a(n+1:end)))];
    
    %powerUpperbound(:,maxUpperBoundIndex) = p_k(:,k);
    powerFeasible = p_k(:,k);
    cbv = prod([(GammaInitlinear(1:n)+powerFeasible(1:n)*hs*alpha+powerFeasible(n+1:end)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+powerFeasible(1:n)*hs*alpha); ...
        (GammaInitlinear(n+1:end)+powerFeasible(n+1:end)*hs*alpha+powerFeasible(1:n)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+powerFeasible(n+1:end)*hs*alpha)]);
    if k>3 && norm(p_k(:,k)-p_k(:,k-1))<=0.00001
        break;
    end
    
    totalPower(1,:)=powerFeasible(1:n);
    totalPower(2,:)=powerFeasible(n+1:end);
    
    
    %%Do the linear integer program%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iterations=1:100
        
        f = zeros(1,2*n);
        for i=1:n
            for theta=1:2
                if theta==2
                    f(1,(i-1)*2+1)=Xinit(i*2,1)*(GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(theta,i)*hs*alpha)...
                        /(GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(theta,i)*hs*alpha+totalPower(1,i)*gainVecSUTrnsSURecv)...
                        *totalPower(1,i)*gainVecSUTrnsSURecv/((GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(theta,i)*hs*alpha)^2) ...
                        *totalPower(theta,i)*hs*alpha + log2(1+totalPower(theta,i)*gainVecSUTrnsSURecv/(GammaInit(i,1)+Xinit(i*2,1)*totalPower(1,i)*hs*alpha));
                    A(4+(i-1)*2+1,(i-1)*2+1)= -totalPower(theta,i)*gainVecSUTrnsSURecv/(GammaInit(i,1)+totalPower(1,i)*hs*alpha)+RxSensitivity;
                    b(4+(i-1)*2+1,1) = 0;
                    
                    AUtopia((i-1)*2+1,1) = gainVecSUTrnsPURecv(2,1);
                    bUtopia((i-1)*2+1,1)= max(IntThreshold(i,1)-SUIntAtPU(i,1),0);
                else
                    f(1,i*2)=Xinit((i-1)*2+1,1)*(GammaInit(i,1)+Xinit(i*2,1)*totalPower(theta,i)*hs*alpha)...
                        /(GammaInit(i,1)+Xinit(i*2,1)*totalPower(theta,i)*hs*alpha+totalPower(2,i)*gainVecSUTrnsSURecv)...
                        *totalPower(2,i)*gainVecSUTrnsSURecv/((GammaInit(i,1)+Xinit(i*2,1)*totalPower(theta,i)*hs*alpha)^2)...
                        *totalPower(theta,i)*hs*alpha + log2(1+totalPower(theta,i)*gainVecSUTrnsSURecv/(GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(2,i)*hs*alpha));
                    A(4+i*2,i*2)= -totalPower(theta,i)*gainVecSUTrnsSURecv/(GammaInit(i,2)+totalPower(2,i)*hs*alpha)+RxSensitivity;
                    b(4+i*2,1) = 0;
                    AUtopia(i*2,1) = gainVecSUTrnsPURecv(1,1);
                    bUtopia(i*2,1)= max(IntThreshold(i,1)-SUIntAtPU(i,1),0);
                end
            end
            A(2*n+4+i,(i-1)*2+1:(i-1)*2+2) = [totalPower(2,i)*gainVecSUTrnsPURecv(2,1),totalPower(1,i)*gainVecSUTrnsPURecv(1,1)];
            b(2*n+4+i,1)=IntThreshold(i,1)-SUIntAtPU(i,1);
        end
        
        A(1,:)=-kron(ones(1,n),[1,0]);
        A(2,:)=-kron(ones(1,n),[0,1]);
        b(1:2,1) = 0;
        A(3,:)=kron(totalPower(1,:),[0,1]);
        A(4,:)=kron(totalPower(2,:),[1,0]);
        b(3:4,1) = Pmax;
        intcon = 1:2*n;
        lb=0;
        ub=1;
        
        xNew = bintprog(-f,A,b);
        if norm(Xinit-xNew)<0.1
            break;
        end
        Xinit=xNew;
    end
    
    for i=1:n
        x_k(n+i,1)=xNew((i-1)*2+1); % group B
        x_k(i,1)=xNew(i*2); % group A
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if abs(maxCUB-cbv)<=epsilon
        break;
    end
    %% new polyblock generation
    colVecz=[];
    colVecy=[];
    rowVector=zeros(2*n,1);
    for col=1:size(powerUpperbound,2)
        for row=1:2*n
            if powerUpperbound(row,col)>p_k(row,k)
                rowVector(row,1)=1;
            elseif powerUpperbound(row,col)==p_k(row,k) && x_k(row,1)==0
                rowVector(row,1)=1;
            end
        end
        if sum(rowVector)==2*n
            colVecz=[colVecz,col];
            colVecy=[colVecy,col];
        elseif all(powerUpperbound(:,col)>=p_k(:,k))
            colVecy=[colVecy,col];
        
        end
    end
    if isempty(colVecz)
        disp(powerFeasible)
        break;
    else
        I = eye(2*n);
        
        for col=fliplr(colVecz)
            rowVec=[];
            z = powerUpperbound(:,col);
            for ycol=fliplr(colVecy)
                comp=find(z>boxesupperbound(:,ycol),2);
                if ~isempty(comp)&& size(comp,1)==1
                    rowVec=[rowVec,comp];
                end
            end
                
            powerUpperbound(:,col)=[];
            boxeslowerbound(:,col)=[];
            cub(col)=[];
            
            
            for row=1:2*n
                if (any(rowVec~=row) || isempty(rowVec))&& x_k(row,1)==1
                    newVertex=z+(p_k(:,k)-z).*I(:,row);
                    upperBoundObj= prod([(GammaInitlinear(1:n)+newVertex(1:n)*hs*alpha+newVertex(n+1:end)*gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+newVertex(1:n)*hs*alpha); ...
                        (GammaInitlinear(n+1:end)+newVertex(n+1:end)*hs*alpha+newVertex(1:n)*gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+newVertex(n+1:end)*hs*alpha)]);
                    if upperBoundObj>cbv
                        powerUpperbound = [powerUpperbound,newVertex];
                        boxeslowerbound = [boxeslowerbound,initPowerFeasible];
                        cub = [cub,upperBoundObj];
                    end
                    
                end
            end
        end
    end
    
    
    removableCols=[];
    for index=1:size(cub,2)
        if all(removableCols~=index)
            for compareIndex=2:size(cub,2)
                if all(powerUpperbound(:,index)>=powerUpperbound(:,compareIndex))
                    removableCols=[removableCols,compareIndex];
                end
            end
        end
    end
    cub(removableCols)=[];
    boxeslowerbound(:,removableCols) =[];
    powerUpperbound(:,removableCols) =[];
    %[~,indices]=sort(cub,'descend');
    if size(cub,2)>=100
        cub=cub(1:100);
        boxeslowerbound =boxeslowerbound(:,1:100);
        powerUpperbound= powerUpperbound(:,1:100);
        boxesupperbound= powerUpperbound;
    else
        boxesupperbound= powerUpperbound;
    end
    
end
% iterations = 1;
% 
% while newSol>oldSol && iterations <1000
%     [gradA,gradB]= gradCalc(GammaInit,hs,alpha,gainVecSUTrnsSURecv,Pinit,Pinit,Xinit');
%     oldSol=newSol;
%     cvx_begin
%     variable p_A(n) nonnegative
%     variable p_B(n) nonnegative
%     minimize -gradA*p_A-gradB*p_B
%     sum(p_A)<=Pmax;
%     sum(p_B)<=Pmax;
%     gainVecSUTrnsSURecv(SUIndex)*p_B - diag(Xinit)*RxSensitivity*hs*alpha*p_A >=(Xinit').*RxSensitivity.*GammaInit(:,1);
%     gainVecSUTrnsSURecv(SUIndex)*p_A - diag(Xinit)*RxSensitivity*hs*alpha*p_B >=(Xinit')*RxSensitivity.*GammaInit(:,2);
%     p_A.*max(gainVecSUTrnsPURecv(1,:))+ p_B.*max(gainVecSUTrnsPURecv(2,:))<=IntThreshold-SUIntAtPU;
%     cvx_end
%     newSol = Xinit*(log(GammaInit(:,2)+alpha*hs*p_B+gainVecSUTrnsSURecv(SUIndex)*p_A))-Xinit*log(GammaInit(:,2) ...
%         +alpha*hs*p_B)+Xinit*(log(GammaInit(:,1)+alpha*hs*p_A+gainVecSUTrnsSURecv(SUIndex)*p_B))-Xinit*log(GammaInit(:,1)+alpha*hs*p_A);
%     
%     [gradA,gradB]= gradCalc(GammaInit,hs,alpha,gainVecSUTrnsSURecv,p_A,p_B,Xinit');
%     iterations=iterations+1;
% end
% disp(p_A)
% disp(p_B)
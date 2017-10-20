classdef SecondaryUserFDJointAlloc < handle
    %This class defines the properties and methods of the secondary user class
    %   OptX function optimizes w.r.t channel allocation
    %   OptP function optimizes w.r.t  power allocation
    %   OptGamma function optimizes w.r.t the interference and noise 
    %   dualMaster function optimizes w.r.t the price of the QoS constraint
    %   PriceOfSUInterRecevd function optimizes w.r.t the interference
    %   received at the given SU
    
    properties
        SUIndex %index of the secondary user pair
        gainVecSUTrnsSURecv  % channel gain from the given transmitting SU to its own receiver
        gainVecSUTrnsMSRecv  % max channel gains from the given transmitting SUs to the monitoring stations
        X                    % channel allocation vec
        P                    % Power allocation Vec
        GammaInit            % interference vec: addition of both SU int and PU int
        alpha                % fraction of residual power after SIC
        maxiterations        % maximum number of iterations
        Pmax                 % Max power usable at an SU
        J                    % interference caused at the PU by SUs in the previous slot
        hs                   % the gain between the SU Tx antenna and Rx antenna
        IntThreshold         %interference threshold at PU
        RxSensitivity        % minimum detectable sinr at detector
        rate                 %data rate of the pair
        
    end
    
    methods
        function SU = SecondaryUserFDJointAlloc(SUIndex,gainVecSUTrnsSURecv,gainVecSUTrnsMSRecv,...
                      GammaInit,Pmax,maxiterations,hs,alpha,SUIntAtPU,IntThreshold,RxSensitivity)
                  SU.SUIndex = SUIndex;
                  SU.gainVecSUTrnsSURecv = gainVecSUTrnsSURecv;
                  SU.gainVecSUTrnsMSRecv = gainVecSUTrnsMSRecv;
                  SU.X = zeros(size(GammaInit,1)*2,1);
                  SU.P = zeros(size(GammaInit,1)*2,1);
                  SU.GammaInit = GammaInit;
                  SU.Pmax = Pmax;
                  SU.maxiterations = maxiterations;
                  SU.hs = hs;
                  SU.alpha = alpha;
                  SU.J = SUIntAtPU;
                  SU.IntThreshold = IntThreshold*ones(size(GammaInit,1),1);
                  SU.RxSensitivity  = RxSensitivity; 
                  SU.rate = 0;
                  
        end
        
        function optXPGamma(SU)
            
            numChannels =size(SU.GammaInit,1);
            n=numChannels;
            Xinit = ones(numChannels*2,1);
            epsilon=0.1;
            
            cbvPrev=-inf;
            
            coeff = [SU.gainVecSUTrnsSURecv(2,SU.SUIndex), -SU.RxSensitivity*SU.hs*SU.alpha];
            coeffMat = zeros(n*2);
            rhsMat = zeros(n*2,1);
            
            for i=1:2*n
                
                if mod(i,2)~= 0
                    coeffMat(i,i:i+1)= coeff;
                    rhsMat(i) = SU.GammaInit(floor(i/2)+1,1);
                else
                    coeffMat(i,i-1:i)= [-SU.RxSensitivity*SU.hs*SU.alpha, SU.gainVecSUTrnsSURecv(2,SU.SUIndex)];
                    rhsMat(i) = SU.GammaInit(i/2,2);
                end
            end
            
            pfeasible = coeffMat\rhsMat;
            totalPower = zeros(2,n);
            powerUpperbound = zeros(2*n,1);
            for i=1:n
                totalPower(2,i)=pfeasible((i-1)*2+1); % group B
                totalPower(1,i)=pfeasible((i)*2);% group A
            end
            
            
            
            
            
            if sum(totalPower(1,:))<=SU.Pmax && sum(totalPower(2,:))<=SU.Pmax
                isFeasible = 1;
            else
                isFeasible = 0;
            end
            
            intAtPU=totalPower(1,:)*SU.gainVecSUTrnsMSRecv(1,1)+totalPower(2,:)*SU.gainVecSUTrnsMSRecv(2,1);
            
            if isFeasible && all(intAtPU<=max(SU.IntThreshold'-SU.J',0))
                isFeasible = 1;
            else
                isFeasible = 0;
            end
            
            
            A = zeros(3*n+2,2*n);
            AUtopia = zeros(2*n,1);
            b = zeros(3*n+2,1);
            bUtopia = zeros(2*n,1);
            xNew = zeros(2*n,1);
            SU.X=zeros(2*n,1);
            Ired = eye(n);
            %% integer linear program
            for iterations=1:100
                
                f = zeros(1,2*n);
                for i=1:n
                    for theta=1:2
                        if theta==2
                            f(1,(i-1)*2+1)= log2(1+totalPower(theta,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)/(SU.GammaInit(i,1)+totalPower(1,i)*SU.hs*SU.alpha));
                            
%                             Xinit(i*2,1)*(SU.GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(theta,i)*SU.hs*SU.alpha)...
%                                 /(SU.GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(theta,i)*SU.hs*SU.alpha+totalPower(1,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))...
%                                 *totalPower(1,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)/((SU.GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(theta,i)*SU.hs*SU.alpha)^2) ...
%                                 *totalPower(theta,i)*SU.hs*SU.alpha + log2(1+totalPower(theta,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)/(SU.GammaInit(i,1)+Xinit(i*2,1)*totalPower(1,i)*SU.hs*SU.alpha));
                             A((i-1)*2+1,(i-1)*2+1)= SU.RxSensitivity;
                             b((i-1)*2+1,1) =totalPower(theta,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)/(SU.GammaInit(i,1)+totalPower(1,i)*SU.hs*SU.alpha);
%                             
                            AUtopia((i-1)*2+1,1) = SU.gainVecSUTrnsMSRecv(2,1);
                            bUtopia((i-1)*2+1,1)= max(SU.IntThreshold(i,1)-SU.J(i,1),0);
                        else
                            f(1,i*2)= log2(1+totalPower(theta,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)/(SU.GammaInit(i,2)+totalPower(2,i)*SU.hs*SU.alpha));
                            
%                             Xinit((i-1)*2+1,1)*(SU.GammaInit(i,1)+Xinit(i*2,1)*totalPower(theta,i)*SU.hs*SU.alpha)...
%                                 /(SU.GammaInit(i,1)+Xinit(i*2,1)*totalPower(theta,i)*SU.hs*SU.alpha+totalPower(2,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))...
%                                 *totalPower(2,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)/((SU.GammaInit(i,1)+Xinit(i*2,1)*totalPower(theta,i)*SU.hs*SU.alpha)^2)...
%                                 *totalPower(theta,i)*SU.hs*SU.alpha + log2(1+totalPower(theta,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)/(SU.GammaInit(i,2)+Xinit((i-1)*2+1,1)*totalPower(2,i)*SU.hs*SU.alpha));
%                             A(i*2+2,i*2)= SU.RxSensitivity;
%                             b(i*2+2,1) = totalPower(theta,i)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)/(SU.GammaInit(i,2)+totalPower(2,i)*SU.hs*SU.alpha);
                            AUtopia(i*2,1) = SU.gainVecSUTrnsMSRecv(1,1);
                            bUtopia(i*2,1)= max(SU.IntThreshold(i,1)-SU.J(i,1),0);
                        end
                    end
                    A(2*n+i+2,(i-1)*2+1:(i-1)*2+2) = [totalPower(2,i)*SU.gainVecSUTrnsMSRecv(2,1),totalPower(1,i)*SU.gainVecSUTrnsMSRecv(1,1)];
                    b(2*n+i+2,1)=max(SU.IntThreshold(i,1)-SU.J(i,1),0);
                end
                
%                 A(1,:)=-kron(ones(1,n),[1,0]);
%                 A(2,:)=-kron(ones(1,n),[0,1]);
%                 b(1:2,1) = 0;
                A(1,:)=kron(totalPower(1,:),[0,1]);
                A(2,:)=kron(totalPower(2,:),[1,0]);
                b(1:2,1) = SU.Pmax;
%                 intcon = 1:2*n;
%                 lb=0;
%                 ub=1;
                
                xNew = bintprog(-f,A,b,[],[],Xinit);
                if norm(Xinit-xNew)<0.1
                    break;
                end
                Xinit=xNew;
            end
            %%lower and the upper bound of power
            
            boxesupperbound = bUtopia./AUtopia.*xNew;
            GammaInitlinear=[SU.GammaInit(:,1);SU.GammaInit(:,2)];
            SU.P = zeros(2*n,1);
            for i=1:n
                powerUpperbound(n+i,1)=max(boxesupperbound((i-1)*2+1),0); % group B
                powerUpperbound(i,1)=max(boxesupperbound((i)*2),0);% group A
                SU.P(n+i,1)=max(pfeasible((i-1)*2+1),0); % group B
                SU.P(i,1)=max(pfeasible((i)*2),0);% group A
                SU.X(n+i,1)=xNew((i-1)*2+1); % group B
                SU.X(i,1)=xNew(i*2); % group A
            end
            
            if isFeasible
                
                
                cub = prod([(GammaInitlinear(1:n)+powerUpperbound(1:n)*SU.hs*SU.alpha+powerUpperbound(n+1:end)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(1:n)+powerUpperbound(1:n)*SU.hs*SU.alpha); ...
                    (GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*SU.hs*SU.alpha+powerUpperbound(1:n)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*SU.hs*SU.alpha)]);
                cbv = prod([(GammaInitlinear(1:n)+SU.P(1:n)*SU.hs*SU.alpha+SU.P(n+1:end)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(1:n)+SU.P(1:n)*SU.hs*SU.alpha); ...
                    (GammaInitlinear(n+1:end)+SU.P(n+1:end)*SU.hs*SU.alpha+SU.P(1:n)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(n+1:end)+SU.P(n+1:end)*SU.hs*SU.alpha)]);
                SU.rate = log2(cbv);
            else
                pfeasible = pfeasible.*xNew;
                
                for i=1:n
                    
                    SU.P(n+i,1)=max(pfeasible((i-1)*2+1),0); % group B
                    SU.P(i,1)=max(pfeasible((i)*2),0);% group A
                end
                
                %cub = prod([(GammaInitlinear(1:n)+powerUpperbound(1:n)*SU.hs*SU.alpha+powerUpperbound(n+1:end)*SU.gainVecSUTrnsSURecv)./(GammaInitlinear(1:n)+powerUpperbound(1:n)*SU.hs*SU.alpha); ...
                %   (GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*SU.hs*SU.alpha+powerUpperbound(1:n)*SU.gainVecSUTrnsSURecv)./(GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*SU.hs*SU.alpha)]);
                cbv = prod([(GammaInitlinear(1:n)+SU.P(1:n)*SU.hs*SU.alpha+SU.P(n+1:end)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(1:n)+SU.P(1:n)*SU.hs*SU.alpha); ...
                    (GammaInitlinear(n+1:end)+SU.P(n+1:end)*SU.hs*SU.alpha+SU.P(1:n)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(n+1:end)+SU.P(n+1:end)*SU.hs*SU.alpha)]);
                SU.rate = log2(cbv);
            end
            cub = prod([(GammaInitlinear(1:n)+powerUpperbound(1:n)*SU.hs*SU.alpha+powerUpperbound(n+1:end)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(1:n)+powerUpperbound(1:n)*SU.hs*SU.alpha); ...
                (GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*SU.hs*SU.alpha+powerUpperbound(1:n)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(n+1:end)+powerUpperbound(n+1:end)*SU.hs*SU.alpha)]);;
            boxeslowerbound = SU.P;
            boxesupperbound = max(min(powerUpperbound,SU.Pmax),0);
            p_k = zeros(2*n,SU.maxiterations);
            
            %%polyblock algorithm
            initPowerFeasible=SU.P;
            for k = 1:SU.maxiterations
                %reduce the current boxes
                
                eraseBoxes=[];
                for boxNum=1:size(boxesupperbound,2)
                    intAtPU=boxeslowerbound(1:n,boxNum)*SU.gainVecSUTrnsMSRecv(1,1)+boxeslowerbound(n+1:end,boxNum)*SU.gainVecSUTrnsMSRecv(2,1);
                    objFuncUpperbound= prod([(GammaInitlinear(1:n)+boxesupperbound(1:n,boxNum)*SU.hs*SU.alpha+boxesupperbound(n+1:end,boxNum)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(1:n)+boxesupperbound(1:n,boxNum)*SU.hs*SU.alpha); ...
                        (GammaInitlinear(n+1:end)+boxesupperbound(n+1:end,boxNum)*SU.hs*SU.alpha+boxesupperbound(1:n,boxNum)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(n+1:end)+boxesupperbound(n+1:end,boxNum)*SU.hs*SU.alpha)]);
                    if sum(boxeslowerbound(1:n,boxNum))<=SU.Pmax && sum(boxeslowerbound(n+1:end,boxNum))<=SU.Pmax && all(intAtPU<=max(SU.IntThreshold-SU.J,0))
                        if all(boxesupperbound(n+1:end,boxNum)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)./(GammaInitlinear(1:n)+boxesupperbound(1:n,boxNum)*SU.hs*SU.alpha)>=SU.RxSensitivity.*SU.X(n+1:end,1)) && ...
                                all(boxesupperbound(1:n,boxNum)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)./(GammaInitlinear(n+1:end)+boxesupperbound(n+1:end,boxNum)*SU.hs*SU.alpha)>=SU.RxSensitivity.*SU.X(1:n,1)) && objFuncUpperbound>cbv
                            %do the reduction
                            alphaVec=[];
                            for individual=1:2
                                for row=1:n
                                    if individual==1 && SU.X(row,1)>0
                                        cvx_begin quiet
                                        variable alphaOpt nonnegative
                                        maximize alphaOpt
                                        subject to
                                        sum(boxeslowerbound(1:n,boxNum)+alphaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row))<=SU.Pmax;
                                        %sum(boxeslowerbound(n+1:end,boxNum)+alphaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))<=SU.Pmax;
                                        (boxeslowerbound(row,boxNum)+boxeslowerbound(n+row,boxNum)+alphaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row))*max(SU.gainVecSUTrnsMSRecv(1,:))<=max(SU.IntThreshold(row,1)-SU.J(row,1),0);
                                        %(boxeslowerbound(n+1:end,boxNum)+alphaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*max(SU.gainVecSUTrnsMSRecv(2,:))<=SU.IntThreshold-SU.J;
                                        alphaOpt <= 1;
                                        cvx_end
                                        if ~isnan(alphaOpt)
                                            alphaVec=[alphaVec,alphaOpt];
                                        else
                                            alphaVec=[alphaVec,1];
                                        end
                                        
                                        boxesupperbound(row,boxNum)=boxeslowerbound(row,boxNum)+alphaVec(row)*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum));
                                    elseif individual==1 && SU.X(row,1)== 0
                                        alphaVec=[alphaVec,1];
                                        boxesupperbound(row,boxNum)=boxeslowerbound(row,boxNum)+alphaVec(row)*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum));
                                    elseif individual==2 && SU.X(n+row,1)>0
                                        cvx_begin quiet
                                        variable alphaOpt nonnegative
                                        maximize alphaOpt
                                        subject to
                                        %sum(boxeslowerbound(1:n,boxNum)+alphaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))<=SU.Pmax;
                                        sum(boxeslowerbound(n+1:end,boxNum)+alphaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))*Ired(:,row))<=SU.Pmax;
                                        %(boxeslowerbound(1:n,boxNum)+alphaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))*max(SU.gainVecSUTrnsMSRecv(1,:))<=SU.IntThreshold-SU.J;
                                        (boxeslowerbound(row,boxNum)+boxeslowerbound(n+row,boxNum)+alphaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))*Ired(:,row))*max(SU.gainVecSUTrnsMSRecv(2,:))<=max(SU.IntThreshold(row,1)-SU.J(row,1),0);
                                        alphaOpt <= 1;
                                        cvx_end
                                        
                                        if ~isnan(alphaOpt)
                                            alphaVec=[alphaVec,alphaOpt];
                                        else
                                            alphaVec=[alphaVec,1];
                                        end
                                        boxesupperbound(n+row,boxNum)=boxeslowerbound(n+row,boxNum)+alphaVec(n+row)*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum));
                                    elseif individual==2 && SU.X(n+row,1)== 0
                                        alphaVec=[alphaVec,1];
                                        boxesupperbound(n+row,boxNum)=boxeslowerbound(n+row,boxNum)+alphaVec(n+row)*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum));
                                    end
                                    
                                end
                            end
                            
                            
                            
                            betaVec=[];
                            for individual=1:2
                                for row=1:n
                                    if individual==1 && SU.X(row,1)>0
                                        cvx_begin quiet
                                        variable betaOpt nonnegative
                                        maximize betaOpt
                                        subject to
                                        sum(boxesupperbound(1:n,boxNum)-betaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row))>=SU.Pmax;
                                        %sum(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))>=SU.Pmax;
                                        (boxesupperbound(row,boxNum)+boxesupperbound(n+row,boxNum))*max(SU.gainVecSUTrnsMSRecv(1,:))-betaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row)*max(SU.gainVecSUTrnsMSRecv(1,:))>=max(SU.IntThreshold(row,1)-SU.J(row,1),0);
                                        %(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*max(SU.gainVecSUTrnsMSRecv(2,:))>=SU.IntThreshold-SU.J;
                                        (boxesupperbound(row,boxNum)-betaOpt*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum))*Ired(:,row))*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)./(GammaInitlinear(n+row)+(boxesupperbound(n+row,boxNum))*SU.hs*SU.alpha)>=SU.RxSensitivity;
                                        %(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*SU.gainVecSUTrnsSURecv-(GammaInitlinear(1:n)+(boxesupperbound(1:n,boxNum)-betaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))*SU.hs*SU.alpha)*SU.RxSensitivity>=0;
                                        betaOpt <= 1;
                                        
                                        cvx_end
                                        if ~isnan(betaOpt)
                                            betaVec=[betaVec,betaOpt];
                                        else
                                            betaVec=[betaVec,1];
                                        end
                                        
                                        boxeslowerbound(row,boxNum)=boxesupperbound(row,boxNum)-betaVec(row)*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum));
                                        
                                    elseif individual==1 && SU.X(row,1)==0
                                        betaVec=[betaVec,1];
                                        boxeslowerbound(row,boxNum)=boxesupperbound(row,boxNum)-betaVec(row)*(boxesupperbound(row,boxNum)-boxeslowerbound(row,boxNum));
                                    elseif individual==2 && SU.X(n+row,1)>0
                                        
                                        cvx_begin quiet
                                        variable betaOpt nonnegative
                                        maximize betaOpt
                                        subject to
                                        %sum(boxesupperbound(1:n,boxNum)-betaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))>=SU.Pmax;
                                        sum(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))*Ired(:,row))>=SU.Pmax;
                                        (boxesupperbound(row,boxNum)+boxesupperbound(n+row,boxNum)-betaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))'*Ired(:,row))*max(SU.gainVecSUTrnsMSRecv(2,:))>=max(SU.IntThreshold(row,1)-SU.J(row,1),0);
                                        %(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*max(SU.gainVecSUTrnsMSRecv(2,:))>=SU.IntThreshold-SU.J;
                                        (boxesupperbound(n+row,boxNum)-betaOpt*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum))'*Ired(:,row))*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)./(GammaInitlinear(row)+(boxesupperbound(row,boxNum))*SU.hs*SU.alpha)>=SU.RxSensitivity;
                                        %(boxesupperbound(n+1:end,boxNum)-betaOpt*(boxesupperbound(n+1:end,boxNum)-boxeslowerbound(n+1:end,boxNum))'*Ired(:,row))*SU.gainVecSUTrnsSURecv-(GammaInitlinear(1:n)+(boxesupperbound(1:n,boxNum)-betaOpt*(boxesupperbound(1:n,boxNum)-boxeslowerbound(1:n,boxNum))'*Ired(:,row))*SU.hs*SU.alpha)*SU.RxSensitivity>=0;
                                        betaOpt <= 1;
                                        
                                        cvx_end
                                        
                                        if ~isnan(betaOpt)
                                            betaVec=[betaVec,betaOpt];
                                        else
                                            betaVec=[betaVec,1];
                                        end
                                        boxeslowerbound(n+row,boxNum)=boxesupperbound(n+row,boxNum)-betaVec(n+row)*(boxesupperbound(n+row,boxNum)-boxeslowerbound(n+row,boxNum));
                                    elseif individual==2 && SU.X(n+row,1)==0
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
                %a_prime = a-norm(boxesupperbound-boxeslowerbound)/4*ones(n*2,1);
                
                %     disp(maxUpperBoundIndex);
                %     disp(z);
                %check if z is feasible if so stop
                intAtPU=z(1:n,1)*SU.gainVecSUTrnsMSRecv(1,1)+z(n+1:end,1)*SU.gainVecSUTrnsMSRecv(2,1);
                objFuncUpperbound= prod([(GammaInitlinear(1:n)+z(1:n,1)*SU.hs*SU.alpha+z(n+1:end,1)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(1:n)+z(1:n,1)*SU.hs*SU.alpha); ...
                    (GammaInitlinear(n+1:end)+z(n+1:end,1)*SU.hs*SU.alpha+z(1:n,1)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(n+1:end)+z(n+1:end,1)*SU.hs*SU.alpha)]);
                if sum(z(1:n,1))<=SU.Pmax && sum(z(n+1:end,1))<=SU.Pmax && all(intAtPU<=max(SU.IntThreshold-SU.J,0))
                    if all(z(n+1:end,1)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)./(GammaInitlinear(1:n)+z(1:n,1)*SU.hs*SU.alpha)>=SU.RxSensitivity.*SU.X(n+1:end,1)) && ...
                            all(z(1:n,1)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex)./(GammaInitlinear(n+1:end)+z(n+1:end,1)*SU.hs*SU.alpha)>=SU.RxSensitivity.*SU.X(1:n,1)) && objFuncUpperbound>cbv
                        SU.P=z;
                        break;
                        
                    end
                end
                
                
                cvx_begin quiet
                variable lambda
                minimize lambda
                subject to
                %z-lambda*(z-a_prime)-SU.P<=0
                sum(z(1:n)-lambda*(z(1:n)-a(1:n)))<=SU.Pmax;
                sum(z(n+1:end)-lambda*(z(n+1:end)-a(n+1:end)))<=SU.Pmax;
                (z(1:n)-lambda*(z(1:n)-a(1:n)))*max(SU.gainVecSUTrnsMSRecv(1,:))+ (z(n+1:end)-lambda*(z(n+1:end)-a(n+1:end)))*max(SU.gainVecSUTrnsMSRecv(2,:))<=max(SU.IntThreshold-SU.J,0);
                %
                %(z_vec(1,:)-lambda*(z_vec(1,:)-a_prime_vec(1,:)))'>=zeros(n,1)
                %(z_vec(2,:)-lambda*(z_vec(2,:)-a_prime_vec(2,:)))'>=zeros(n,1)
                cvx_end
                p_k(:,k) = max([(z(1:n)-lambda*(z(1:n)-a(1:n)));(z(n+1:end)-lambda*(z(n+1:end)-a(n+1:end)))],0);
                
                %powerUpperbound(:,maxUpperBoundIndex) = p_k(:,k);
                SU.P = max(p_k(:,k),0);
                cbvPrev=cbv;
                cbv = prod([(GammaInitlinear(1:n)+SU.P(1:n)*SU.hs*SU.alpha+SU.P(n+1:end)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(1:n)+SU.P(1:n)*SU.hs*SU.alpha); ...
                    (GammaInitlinear(n+1:end)+SU.P(n+1:end)*SU.hs*SU.alpha+SU.P(1:n)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(n+1:end)+SU.P(n+1:end)*SU.hs*SU.alpha)]);
                
                SU.rate = log2(cbv);
%                 if k>3 && abs(cbvPrev-cbv)<=0.00001
%                     break;
%                 end
                powerUpperbound(:,maxUpperBoundIndex)=p_k(:,k);
                totalPower(1,:)=SU.P(1:n);
                totalPower(2,:)=SU.P(n+1:end);
                
                
                
                
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
                        elseif powerUpperbound(row,col)==p_k(row,k) && SU.X(row,1)==0
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
                    disp(p_k(:,k))
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
                            if (any(rowVec~=row) || isempty(rowVec))&& SU.X(row,1)==1
                                newVertex=z+(p_k(:,k)-z).*I(:,row);
                                upperBoundObj= prod([(GammaInitlinear(1:n)+newVertex(1:n)*SU.hs*SU.alpha+newVertex(n+1:end)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(1:n)+newVertex(1:n)*SU.hs*SU.alpha); ...
                                    (GammaInitlinear(n+1:end)+newVertex(n+1:end)*SU.hs*SU.alpha+newVertex(1:n)*SU.gainVecSUTrnsSURecv(2,SU.SUIndex))./(GammaInitlinear(n+1:end)+newVertex(n+1:end)*SU.hs*SU.alpha)]);
                                if upperBoundObj>cbv
                                    powerUpperbound = [powerUpperbound,newVertex];
                                    boxeslowerbound = [boxeslowerbound,initPowerFeasible.*SU.X];
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
        end
    end
    
end


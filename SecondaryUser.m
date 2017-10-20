classdef SecondaryUser < handle
    %This class defines the properties and methods of the secondary user class
    %   OptX function optimizes w.r.t channel allocation
    %   OptP function optimizes w.r.t  power allocation
    %   OptGamma function optimizes w.r.t the interference and noise 
    %   dualMaster function optimizes w.r.t the price of the QoS constraint
    %   PriceOfSUInterRecevd function optimizes w.r.t the interference
    %   received at the given SU
    
    properties
        SUIndex %index of the secondary user pair
        gainVecSUTrnsSURecv  % channel gain from the given transmitting SU to other SU receivers including its own receiver
        gainVecSUTrnsPURecv  % channel gain from the given transmitting SU to other PU receivers including its own receiver
        gainVecSURecvPUTrns  % channel gain from the given transmitting PU to the given SU receiver last time slot
        X                    % channel allocation vec
        P                    % Power allocation Vec
        Gamma                % interference vec (auxiliary variable)
        alpha                % residual power after SIC
        maxiterations        % maximum number of iterations
        step                 % step size used in subgradient method
        Pmax                 % Max power usable at an SU
        J                    % interference cased at the PU by SUs
        hs                   % the gain between the SU Tx antenna and Rx antenna
        Threshold            %interference threshold at PU
        RxSensitivity        % minimum detectable sinr at detector
        
    end
    
    methods
        function SU = SecondaryUser(SUIndex,gainVecSUTrnsSURecv,gainVecSUTrnsPURecv,Xinit,...
                      Pinit,GammaInit,Pmax,maxiterations,hs,alpha,SUIntAtPU,IntThreshold,RxSensitivity)
                  SU.SUIndex = SUIndex;
                  SU.gainVecSUTrnsSURecv = gainVecSUTrnsSURecv;
                  SU.gainVecSUTrnsPURecv = gainVecSUTrnsPURecv;
                  SU.X = Xinit;
                  SU.P = Pinit;
                  SU.Gamma = GammaInit;
                  SU.Pmax = Pmax;
                  SU.maxiterations = maxiterations;
                  SU.hs = hs;
                  SU.alpha = alpha;
                  SU.J = SUIntAtPU;
                  SU.Threshold = IntThreshold;
                  SU.RxSensitivity  = RxSensitivity; 
                  
        end
        
        function optXPGamma(SU)
            
            %   Detailed explanation goes here
            %maxiterations = 5000;
            Xold = zeros(1,length(SU.P));
            Pold = zeros(1,length(SU.X));
            Gammaold = zeros(1,length(SU.X));
            iteration = 0;
            oldObjective = 0; 
                       
            while (oldObjective < objective || iteration <= SU.maxiterations)
                
                %% optimization w.r.t. channel allocation vector
                SU.betaX = SU.betaX+SU.step*(sum(SU.X)-SU.K);
                SU.X = ones(length(SU.X),1);
                Xold = SU.X;
                
                
                %% optimizing w.r.t Power allocation
                cvx_begin
                variable p_A(n) nonnegative
                variable p_B(n) nonnegative
                maximize SU.X*(log(SU.Gamma(2,:)+SU.alpha*SU.hs*p_B+SU.gainVecSUTrnsSURecv(SU.SUIndex)*p_A))-SU.X*log(SU.Gamma(2,:) ...
                        +SU.alpha*SU.hs*p_B)+SU.X*(log(SU.Gamma(1,:)+SU.alpha*SU.hs*p_A+SU.gainVecSUTrnsSURecv(SU.SUIndex)*p_B))-SU.X*log(SU.Gamma(1,:)+SU.alpha*SU.hs*p_A)
                sum(p_A)<=SU.Pmax;
                sum(p_B)<=SU.Pmax;
                SU.gainVecSUTrnsSURecv(SU.SUIndex)*p_B - diag(SU.X)*SU.RxSensitivity*SU.hs*SU.alpha*p_A >=SU.X*SU.RxSensitivity*SU.Gamma(1,:);
                SU.gainVecSUTrnsSURecv(SU.SUIndex)*p_A - diag(SU.X)*SU.RxSensitivity*SU.hs*SU.alpha*p_B >=SU.X*SU.RxSensitivity*SU.Gamma(2,:);
                p_A.*SU.gainVecSUTrnsPURecv(1,:)'+ p_B.*SU.gainVecSUTrnsPURecv(2,:)'<=SU.Threshold-SU.J
                cvx_end
                disp(p_A,p_B)
                iteration =iteration+ 1;
            end
            
        end
    end
    
end


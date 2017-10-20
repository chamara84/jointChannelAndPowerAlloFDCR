function [xhatK,PK] = DkalmanUpdate(A,Q,H,R,zK,PKneg1,xhatKneg1)

%xK = A xKneg1 + B uKneg1 + wKneg1
%zK = H xK + vK
%wK~N(0,Q)
%vK~N(0,R)
%PK = a posteriori estimate error covariance
%xhatKneg1 = a posteriory estimate of state

%% Time update equations
xhatMinusK = A*xhatKneg1;
PMinusK = A*PKneg1*A'+ Q;


%% Measurement update equations
if rcond(H*PMinusK*H'+R)>1e-6 || ~isnan(rcond(H*PMinusK*H'+R))
    KK = PMinusK*H'*inv(H*PMinusK*H'+R);
else
    disp(H);
    %disp(R);
    %disp(Q);
    KK = 0*PMinusK*H'*eye(size(R,1));
end

xhatK = xhatMinusK + KK*(zK-H*xhatMinusK);
xhatK = max(xhatK,0);
PK = (eye(size(KK*H,1))-KK*H)*PMinusK;


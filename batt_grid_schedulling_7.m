%% Social welfare schedulling program (c) 2023
% Prof. Paulo M. De Oliveira pdeoliv@gmail.com
%%
clear all
clc
close all
global T epsilon2 epsilon3 lambda1 lambda2_max lambda3_max Smax...
       G B PF QD2 QD3 PG1 QG1 S12 S13 S23 P12 P13 P23 P21 P31 P32...
       Q12 Q13 Q23 Q21 Q31 Q32 alpha beta gamma ...


       reply = input('Reactive power optimization? Y/N [Y]:','s');
       if isempty(reply)
          reply = 'Y';
       end
%        reply2 = input('Congestion considered? Y/N [Y]:','s');
%        if isempty(reply2)
           reply2 = 'N';
%        end
Sbase=10;%MVA
Vbase=69;%kV
Ebase=10;%MWh
%% SIMULATION DATA
T=24;%Hours
%% WHOLESALE SPOT MARKET
% lambda1 =[1004.14339092153,1004.14339092153,1000,1000,1000,...
%     1004.14339092153,1004.14339092153,1004.14339092153,...
%     1004.14339092153,1004.14339092153,1040.65702341756,...
%     1077.68857977877,1077.68857977877,1040.65702341756,...
%     1040.65702341756,1014.90880840516,1004.14339092153,...
%     1004.14339092153,1040.65702341756,1077.68857977877,...
%     1077.68857977877,1014.90880840516,1004.14339092153,1000];%pu XM data
 lambda1=[1042.90 957.10 957.10 952.81 940.59 947.19 1000.00...
 1196.04 1319.14 1324.09 1273.60 1179.87 1221.12 1217.16 1162.38...
 1133.00 1101.32 1244.88 1333.66 1353.47 1359.74 1358.42 1314.52 1148.51];% OMIE data
%% LOAD MODEL
lambda2_max  = 10000;%$/pu.h
lambda3_max  = 10000;%$/pu.h
% loadcurve=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
% for k=1:T
% epsilon2(k)=(lambda2_max-lambda1(k))/0.5*loadcurve(k);%$/pu^2.h
% epsilon3(k)=(lambda3_max-lambda1(k))/3.5*loadcurve(k);%$/pu^2.h
% PD2max(k)=lambda2_max/epsilon2(k);
% PD3max(k)=lambda3_max/epsilon3(k);
% PD2_0(k)=(lambda2_max-lambda1(k))/epsilon2(k); 
% PD3_0(k)=(lambda3_max-lambda1(k))/epsilon3(k); 
% elasticity2_min(k)=-inv(epsilon2(k)*PD2_0(k)/lambda1(k));%inelastic
% elasticity3_min(k)=-inv(epsilon3(k)*PD3_0(k)/lambda1(k));%inelastica
% end
%% LOAD MODEL
lambda2_max  = 10000;%$/pu.h
lambda3_max  = 10000;%$/pu.h
epsilon2=(lambda2_max-1000)/0.5;%$/pu^2.h
epsilon3=(lambda3_max-1000)/3.5;%$/pu^2.h
PD2max=lambda2_max/epsilon2;
PD3max=lambda3_max/epsilon3;
PD2_0min=(lambda2_max-min(lambda1))/epsilon2; 
PD3_0min=(lambda3_max-min(lambda1))/epsilon3; 
elasticity2_min=-inv(epsilon2*PD2_0min/min(lambda1));%inelastic
elasticity3_min=-inv(epsilon3*PD3_0min/min(lambda1));%inelastica
PD2_0max=(lambda2_max-max(lambda1))/epsilon2;
PD3_0max=(lambda3_max-max(lambda1))/epsilon3;
elasticity2_max=-inv(epsilon2*PD2_0max/max(lambda1));
elasticity3_max=-inv(epsilon3*PD3_0max/max(lambda1));
PF = 0.8;%power factor (lag)
%% WIND GENERATOR
% PG2forecast = [0.63, 0.58, 0.63, 0.68, 0.79, 0.79, 0.84, 0.84, 1 , 0.79, 0.68, ...
%        0.58, 0.32, 0.74, 0.95, 0.89, 0.79, 0.53, 0.26, 0.32, 0.53, 0.47, ...
%        0.58, 0.63]*50/Sbase;%pu wind producer
%% SOLAR PV Generator   
PG2forecast = [ 0 0 0 0 0 0 0.811697096 1.82361281 2.884230349...
 3.804153725 4.480513858 5.043290512 5.254331757 4.848483209 4.588740138 3.809565039...
 2.873407721 1.829024124 0.768406584 0 0 0 0 0];  
PG2forecast=PG2forecast';
lambda1=lambda1';
beta=0;alpha=0;% OPEX
gamma=0; %CAPEX
%% 69 kV Network model
% Phase: 266.8 MCM 65/7 ACSR 460A 
% Shield wire/neutral: Copperweld 3/8"
% ATP Analisis: 69kV_7.acp
%
%   G *<--3m--> --
%     |         
%     |________ 6m 
%     |     A ! --
%     |________ 3m   
%     |     B ! --
%     |________ 3m   
%     |     C ! --
%     |
%     |         15.2m
%     |
%     |___________
%      \\\\\\\\\\\
%
% Sbase=10MVA, Vbase=69kV
zline=0.240544+0.481652i;%ohms/km
Zbase=Vbase^2/Sbase;%ohms
 zpu=zline/Zbase;%pu/km
 L12=30;%km
 L13=10;%km
 L23=15;%km
 y12=inv(L12*zpu);%siemens
 y13=inv(L13*zpu);%siemens
 y23=inv(L23*zpu);%siemens
Ybus=[y12+y13 -y12 -y13;-y12 y12+y23 -y23;-y13 -y23 y13+y23];%siemens
G=real(Ybus);
B=imag(Ybus);
if reply2=='Y'
    disp('Congestion enabled')
 Smax=2;   
else
Smax=0.8*sqrt(3)*.460*Vbase/Sbase;%pu
end
%% BATTERY MODEL
C_3 = 10;%pu
rho1=0.2;% RHO pct charge2
rho2=1;% RHO pct charge2
PG3max=3;%pu
lfactor=C_3/(24*PG3max);
HSE=lfactor*24;
%% OPTIMIZATION BEGINS HERE
time000=cputime;
%Initial conditions
     %1*T  2*T  3*T   4*T  5*T    6*T  7*T     8*T
   %24Pg3 24Qg3 24Qg2 24V2 24V3   24T2 24T3    25E3  /    
x0 = [zeros(1,T*3)    ones(1,T*2) zeros(1,T*2) ones(1,T+1)...
      2*ones(1, T*2) 100*ones(1, T*2)  zeros(1,T) ]; 
%     24PD2 24PD3    24L2   24L3     24PG2
%     9*T+1 10*T+1   11*T+1 12*T+1   13*T+1
%Bounds
ub = Inf(1,12*T+1);
lb = -Inf(1,12*T+1);
%% Reactive power capability
if reply=='Y'
    disp('Reactive control enabled')
else
lb(1*T+1:2*T) = 0; % QG2min reactive power
ub(1*T+1:2*T) = 0; % QG2max reactive power
lb(2*T+1:3*T) = 0; % QS3min reactive power
ub(2*T+1:3*T) = 0; % QS3max reactive power
end
%%
ub(7*T+1:8*T+1) = C_3*rho2; % Batt capacity
lb(7*T+1:8*T+1) = C_3*rho1; %min batt usage
lb(0*T+1:1*T) = -PG3max; % Batt max power
ub(0*T+1:1*T) = PG3max; % Batt max power
lb(3*T+1:5*T) = 0.8; % v  lower bound
ub(3*T+1:5*T) = 1.2; % v  lower bound
lb(8*T+2:10*T+1) = 0; % no negative demand
ub(8*T+2:9*T+1) = lambda2_max./epsilon2;% PD2max
ub(9*T+2:10*T+1) = lambda3_max./epsilon3;% PD3max
lb(10*T+2:11*T+1) = 0;%lambda min2=0 
lb(11*T+2:12*T+1) = 0;%lambda min3=0
ub(10*T+2:11*T+1) = lambda2_max;%lambda max2 
ub(11*T+2:12*T+1) = lambda3_max;%lambda max3
lb(12*T+2:13*T+1) = 0;%min wind production in node 2 
ub(12*T+2:13*T+1) = PG2forecast;%max wind production in node 2 
%Equality linear constraints - Storage model
Aeq = zeros(T+1, 13*T+1);
beq = zeros(1, T+1);
Aeq(1,7*T+1)=1;
Aeq(1,8*T+1)=-1;% t=1 E1-E25=0 
for t = 2:T+1
   Aeq(t,7*T+(t-1))=1; %Et
   Aeq(t,(t-1))=-1; %-Pg3 
   Aeq(t,7*T+(t))=-1;%-Et+1
end
%FMINCON calculation
options = optimoptions('fmincon');
options.MaxFunctionEvaluations = 5000000;
options.ConstraintTolerance = 1.0000e-12;
options.MaxIterations = 100000;
options.OptimalityTolerance = 1.0000e-12;
options.StepTolerance = 1.0000e-20;
options.Display='iter';
options.Algorithm='interior-point';
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@objective_func, x0, [], [], Aeq, beq, lb, ub, @network_model, options);
elapsedtime000=cputime-time000 % Set simulation time
exitflag;
%Lagrange multipliers
for k=1:length(lambda.eqnonlin)/2
lambdap(k,1)=lambda.eqnonlin(2*(k-1)+1);%
end
%% Results
R(1,:) = x(0*T+1:1*T)*Sbase;%PS3
R(2,:) = x(1*T+1:2*T)*Sbase;%QG2
R(3,:) = x(2*T+1:3*T)*Sbase;%QS3
R(4,:) = x(3*T+1:4*T);%V2 
R(5,:) = x(4*T+1:5*T);%V3 
R(6,:) = x(5*T+1:6*T);% Th2
R(7,:) = x(6*T+1:7*T);% Th3 
R(8,:) = x(8*T+2:9*T+1)*Sbase;% PD2
R(9,:) = x(9*T+2:10*T+1)*Sbase;% PD3
R(10,:)= x(10*T+2:11*T+1)/Sbase;%lambda2 
R(11,:) = x(11*T+2:12*T+1)/Sbase;%lambda3
R(12,:) = PG1*Sbase;% PG1*Sbase
R(13,:) = QD2*Sbase;% QD2*Sbase
R(14,:) = QD3*Sbase;% QD3
R(15,:)=x(7*T+1:8*T)*Ebase;% SOC
R(16,:) = QG1*Sbase;% QG1
R(17,:) = PG2forecast*Sbase;% PG2 forecast
R(18,:) = x(12*T+2:13*T+1)*Sbase;% PG2 dispached
R(19,:) = S12*Sbase;% S12
R(20,:) = S13*Sbase;% S13 
R(21,:) = S23*Sbase;% S23
R(22,:) = lambda1/Sbase;% lambda1
R(23,:) = (P12+P21+P31+P13+P23+P32)*Sbase;% LossP
R(24,:) = (Q12+Q21+Q31+Q13+Q23+Q32)*Sbase;% LossQ
PS3 = x(0*T+1:1*T);%PS3
PG2=x(12*T+2:13*T+1);
PD3= x(9*T+2:10*T+1);% PD3
PD2=x(8*T+2:9*T+1);% PD2
lambda2=x(10*T+2:11*T+1);
lambda3=x(11*T+2:12*T+1);
E(1,:)=x(7*T+1:8*T+1);% SOC
figure
plot(E)
hold
plot(PS3)

%% Analysis
PSP1=sum(lambda1.*PG1); %Market profit
SellPG1=0;BuyPG1=0;ESellPG1=0;EBuyPG1=0;
for t=1:T
    if PG1(t) > 0
        SellPG1=SellPG1+PG1(t)*lambda1(t); %Eur/day
        ESellPG1=ESellPG1+PG1(t); %puh/day      
    else
        BuyPG1=BuyPG1+PG1(t)*lambda1(t); %Eur/day
        EBuyPG1=EBuyPG1+PG1(t); %puh/day
    end
end

SellPG2=sum(lambda2.*PG2);%Eur/day
EG2=sum(PG2);%puh/day
ED2=sum(PD2);%puh/day
ED3=sum(PD3);%puh/day
CostPG2=sum(0.5*beta.*PG2.^2+alpha.*PG2); 
PSP2=SellPG2-CostPG2;%Wind farm utility
PSP3=sum(lambda3.*PS3); %Battery Uyility
SellPG3=0;BuyPG3=0;ESellPG3=0;EBuyPG3=0;
for t=1:T
    if PS3(t) > 0
        SellPG3=SellPG3+PS3(t)*lambda3(t);%Eur/day
        ESellPG3=ESellPG3+PS3(t);%puh/day
    else
        BuyPG3=BuyPG3+PS3(t)*lambda3(t);%Eur/day
        EBuyPG3=EBuyPG3+PS3(t);%puh/day
    end
end

DemandSurplus=  sum(lambda2_max*PD2-(epsilon2.*(PD2.^2)/2))+...
    sum(lambda3_max*PD3-(epsilon3.*(PD3.^2)/2))-...
sum(lambda2.*PD2+lambda3.*PD3); 
DemandUtility=sum(lambda2_max*PD2-(epsilon2.*(PD2.^2)/2))+...
    sum(lambda3_max*PD3-(epsilon3.*(PD3.^2)/2));
DemandCost=sum(lambda2.*PD2+lambda3.*PD3);

NetworkSurplus=sum(lambda1'.*(-PG1)+lambda2.*(PD2-PG2)+lambda3.*(PD3-PS3));
TPSP=sum(SellPG1)+PSP2+sum(SellPG3+BuyPG3);
SWelfare=DemandUtility-sum(SellPG1+BuyPG1)-CostPG2;

R2(1,:)=-fval;
R2(2,:)=SellPG2;
R2(3,:)=CostPG2;
R2(4,:)=PSP2;
R2(5,:)=sum(SellPG3);
R2(6,:)=sum(-BuyPG3);
R2(7,:)=sum(SellPG3+BuyPG3);
R2(8,:)=DemandUtility;
R2(9,:)=DemandCost;
R2(10,:)=DemandSurplus;
R2(11,:)=NetworkSurplus;
R2(12,:)=sum(SellPG1);
R2(13,:)=sum(-BuyPG1);
R2(14,:)=EG2*Sbase;
R2(15,:)=ESellPG3*Sbase;
R2(16,:)=ESellPG1*Sbase;
R2(17,:)=EBuyPG3*Sbase;
R2(18,:)=-(ED2+ED3)*Sbase;
R2(19,:)=EBuyPG1*Sbase;
R2(20,:)=Sbase*(-sum(P12+P21)-sum(P13+P31)-sum(P23+P32));

 
 
 disp('*******************************************************')
fprintf('Optimization results:\n') 
fprintf('Robust Energy Community Social Welfare %6.2f Eur/day\n',-fval)
fprintf(' \n')
fprintf('-----------Economic Balance---------------------\n')
fprintf('Solar PV Energy Sold                Eur/day %6.2f \n',SellPG2)
fprintf('Solar PV Cost                       Eur/day %6.2f \n',CostPG2)
fprintf('Solar PV Surplus                    Eur/day %6.2f \n',PSP2)
fprintf('Energy sold by the storage          Eur/day %6.2f \n',sum(SellPG3))
fprintf('Energy bought by the storage        Eur/day %6.2f \n',sum(-BuyPG3))
fprintf('Storage Surplus                     Eur/day %6.2f \n',sum(SellPG3+BuyPG3))
fprintf('Demand Benefit                      Eur/day %6.2f \n',DemandUtility)
fprintf('Demand Energy bougth                Eur/day %6.2f \n',DemandCost)
fprintf('Consumer Surplus                    Eur/day %6.2f \n',DemandSurplus)
fprintf('Network Surplus                     Eur/day %6.2f \n',NetworkSurplus)
fprintf('Energy sold by the spot market      Eur/day %6.2f \n',sum(SellPG1))
fprintf('Energy bought by the spot market    Eur/day %6.2f \n',sum(-BuyPG1))
fprintf(' \n')
fprintf('-----------Energy Balance-----------------------\n')
fprintf('Energy Injected by the PV           MWh/day %6.2f \n',EG2*Sbase)
fprintf('Energy Injected by the storage      MWh/day %6.2f \n',ESellPG3*Sbase)
fprintf('Energy Injected by the spot market  MWh/day %6.2f \n',ESellPG1*Sbase)
fprintf('                                            -------\n')
fprintf('Total Energy Injected               MWh/day %6.2f \n',ESellPG1*Sbase+ESellPG3*Sbase+EG2*Sbase)
fprintf('  \n')
fprintf('Energy Consumed by the storage      MWh/day %6.2f \n',EBuyPG3*Sbase)
fprintf('Energy Consumed by the demand       MWh/day %6.2f \n',-(ED2+ED3)*Sbase)
fprintf('Energy Consumed by the spot market  MWh/day %6.2f \n',EBuyPG1*Sbase)
fprintf('Energy Losses                       MWh/day %6.2f \n',Sbase*(-sum(P12+P21)-sum(P13+P31)-sum(P23+P32)))
fprintf('                                            -------\n')
fprintf('Total Energy Consumed               MWh/day %6.2f \n',-(ED2+ED3)*Sbase+...
EBuyPG3*Sbase+EBuyPG1*Sbase+Sbase*(-sum(P12+P21)-sum(P13+P31)-sum(P23+P32)));


function [f] = objective_func(x)
global T epsilon2 epsilon3 lambda1 lambda2_max lambda3_max G B...
    PG1 QG1 PD2 PD3 P12 P13 P23 P21 P31 P32 alpha beta gamma
PD2 = x(8*T+2:9*T+1);
PD3 = x(9*T+2:10*T+1);
PG2 = x(12*T+2:13*T+1);
lambda2=x(10*T+2:11*T+1);
lambda3=x(11*T+2:12*T+1);
PS3 = x(0*T+1:1*T);%PS3
V2 = x(3*T+1:4*T);
V3 = x(4*T+1:5*T);
V1 = ones(1,length(V2));%slack bus
Th2 = x(5*T+1:6*T);
Th3 = x(6*T+1:7*T);
Th1 = zeros(1,length(Th2));
V = [V1;V2;V3];
Th = [Th1;Th2;Th3];
for i = 1:T
    flag = 0;
    for k = 1:3
        flag=flag+V(k,i)*(G(1,k)*cos(Th(1,i)-Th(k,i))+B(1,k)*sin(Th(1,i)-Th(k,i)));
    end
    PG1(i) = V(1,i)*flag;
        flag2 = 0;
    for k = 1:3
        flag2=flag2+V(k,i)*(G(1,k)*sin(Th(1,i)-Th(k,i))-B(1,k)*cos(Th(1,i)-Th(k,i)));
    end
    QG1(i) = V(1,i)*flag2;
end
SW = 0;
for t = 1:T
    U_2 = lambda2_max*PD2(t)-(epsilon2*(PD2(t)^2)/2);
    U_3 = lambda3_max*PD3(t)-(epsilon3*(PD3(t)^2)/2);
    CD_2= lambda2(t)*PD2(t);
    CD_3= lambda3(t)*PD3(t);
    CSP_2=U_2-CD_2;%Consumer 2 surplus
    CSP_3=U_3-CD_3;%Consumer 3 surplus
    PSP_1 = lambda1(t)*PG1(t);%Wholesale market surplus 
    I_2=lambda2(t)*PG2(t); 
    C_2=((1/2)*beta*PG2(t)^2+alpha*PG2(t)+gamma);
    PSP_2 = I_2-C_2;%Generator 2 surplus
    PSP_3= lambda3(t)*PS3(t);%Battery surplus 
    NetworkSurplus=lambda1(t)*(-PG1(t))+lambda2(t)*(PD2(t)-PG2(t))...
        +lambda3(t)*(PD3(t)-PS3(t)); %Network surplus  
   %
SW = SW+CSP_2+CSP_3+PSP_2+PSP_3+NetworkSurplus;
end
f=-SW;%minimization of the social cost
for t = 1:T
    P12(t) = V(1, t) * V(2, t) * (G(1,2) * cos( Th(1, t)-Th(2, t) ) + B(1,2) * sin( Th(1, t)-Th(2, t) )) - G(1,2) * V(1, t)^2;
    P21(t) = V(1, t) * V(2, t) * (G(1,2) * cos( Th(2, t)-Th(1, t) ) + B(1,2) * sin( Th(2, t)-Th(1, t) )) - G(1,2) * V(2, t)^2;
    P13(t) = V(1, t) * V(3, t) * (G(1,3) * cos( Th(1, t)-Th(3, t) ) + B(1,3) * sin( Th(1, t)-Th(3, t) )) - G(1,3) * V(1, t)^2;
    P31(t) = V(1, t) * V(3, t) * (G(1,3) * cos( Th(3, t)-Th(1, t) ) + B(1,3) * sin( Th(3, t)-Th(1, t) )) - G(1,3) * V(3, t)^2;  
    P23(t) = V(2, t) * V(3, t) * (G(2,3) * cos( Th(2, t)-Th(3, t) ) + B(2,3) * sin( Th(2, t)-Th(3, t) )) - G(2,3) * V(2, t)^2;
    P32(t) = V(2, t) * V(3, t) * (G(2,3) * cos( Th(3, t)-Th(2, t) ) + B(2,3) * sin( Th(3, t)-Th(2, t) )) - G(2,3) * V(3, t)^2;
end

%losses=sum(P12+P21+P31+P13+P23+P32)%minimal active losses
end
function [c,ceq] = network_model(x)
global T epsilon2 epsilon3 lambda2_max lambda3_max Smax...
       G B PF QD2 QD3 S12 S13 S23 P12 P13 P23 P21 P31 P32...
       Q12 Q13 Q23 Q21 Q31 Q32 PS3
V2 = x(3*T+1:4*T);
V3 = x(4*T+1:5*T);
V1 = ones(1,length(V2));
Th2 = x(5*T+1:6*T);
Th3 = x(6*T+1:7*T);
Th1 = zeros(1,length(Th2));
V=[V1;V2;V3];
Th = [Th1;Th2;Th3];
PD2 = x(8*T+2:9*T+1);
PD3 = x(9*T+2:10*T+1);
PG2 = x(12*T+2:13*T+1);
QD2 = PD2*tan(acos(PF));
QD3 = PD3*tan(acos(PF));
lambda2 = x(10*T+2:11*T+1);
lambda3 = x(11*T+2:12*T+1);
ceq1 = zeros(1,4*T);
for t = 1:T
    %%
    term1 = 0;
    for k = 1:3
        term1 = term1 + V(k, t)*( G(2,k)*cos(Th(2, t)-Th(k, t)) + B(2,k)*sin(Th(2, t)-Th(k, t)) );
    end 
     P2_sp = -PG2(t) + PD2(t) + V(2, t)*term1;  
    %%
    PS3 = x(1:T);
    flag2 = 0;
    for k = 1:3
    flag2=flag2+V(k,t)*(G(3,k)*cos(Th(3,t)-Th(k,t))+B(3,k)*sin(Th(3,t)-Th(k,t)));
    end
    P3_sp = -PS3(t) + PD3(t) + V(3,t)*flag2;
    %%
    QG2 = x(2*T+1:3*T);
    flag3 = 0;
    for k = 1:3
    flag3=flag3+V(k,t)*(G(2,k)*sin(Th(2,t)-Th(k,t))-B(2,k)*cos(Th(2,t)-Th(k,t)));
    end
    Q2_sp = -QG2(t) + QD2(t) + V(2, t)*flag3;
%%
    QS3 = x(T+1:2*T);
    flag4 = 0;
    for k = 1:3
    flag4=flag4+V(k,t)*(G(3,k)*sin(Th(3,t)-Th(k,t))-B(3,k)*cos(Th(3,t)-Th(k,t)));
    end
    Q3_sp = -QS3(t) + QD3(t) + V(3, t)*flag4;
    ceq1(4*t - 3) = P2_sp; 
    ceq1(4*t - 2) = Q2_sp;
    ceq1(4*t - 1) = P3_sp;
    ceq1(4*t) = Q3_sp;
end
ceq2 = zeros(1,2*T);
for t = 1:T
    ceq2(2*t - 1) = -lambda2(t) + lambda2_max - (epsilon2*PD2(t));
    ceq2(2*t) = -lambda3(t) + lambda3_max - (epsilon3*PD3(t)); 
end
ceq = [ceq1 ceq2];

%% CAPACITY CONSTRAINTS
nq = zeros(1,3*T);
for t = 1:T
    P12(t) = V(1, t) * V(2, t) * (G(1,2) * cos( Th(1, t)-Th(2, t) ) + B(1,2) * sin( Th(1, t)-Th(2, t) )) - G(1,2) * V(1, t)^2;
    P21(t) = V(1, t) * V(2, t) * (G(1,2) * cos( Th(2, t)-Th(1, t) ) + B(1,2) * sin( Th(2, t)-Th(1, t) )) - G(1,2) * V(2, t)^2;
    P13(t) = V(1, t) * V(3, t) * (G(1,3) * cos( Th(1, t)-Th(3, t) ) + B(1,3) * sin( Th(1, t)-Th(3, t) )) - G(1,3) * V(1, t)^2;
    P31(t) = V(1, t) * V(3, t) * (G(1,3) * cos( Th(3, t)-Th(1, t) ) + B(1,3) * sin( Th(3, t)-Th(1, t) )) - G(1,3) * V(3, t)^2;  
    P23(t) = V(2, t) * V(3, t) * (G(2,3) * cos( Th(2, t)-Th(3, t) ) + B(2,3) * sin( Th(2, t)-Th(3, t) )) - G(2,3) * V(2, t)^2;
    P32(t) = V(2, t) * V(3, t) * (G(2,3) * cos( Th(3, t)-Th(2, t) ) + B(2,3) * sin( Th(3, t)-Th(2, t) )) - G(2,3) * V(3, t)^2;
    Q12(t) = V(1, t) * V(2, t) * (G(1,2) * sin( Th(1, t)-Th(2, t) ) - B(1,2) * cos( Th(1, t)-Th(2, t) )) + B(1,2) * V(1, t)^2;
    Q21(t) = V(1, t) * V(2, t) * (G(1,2) * sin( Th(2, t)-Th(1, t) ) - B(1,2) * cos( Th(2, t)-Th(1, t) )) + B(1,2) * V(2, t)^2;
    Q13(t) = V(1, t) * V(3, t) * (G(1,3) * sin( Th(1, t)-Th(3, t) ) - B(1,3) * cos( Th(1, t)-Th(3, t) )) + B(1,3) * V(1, t)^2;
    Q31(t) = V(1, t) * V(3, t) * (G(1,3) * sin( Th(3, t)-Th(1, t) ) - B(1,3) * cos( Th(3, t)-Th(1, t) )) + B(1,3) * V(3, t)^2;
    Q23(t) = V(2, t) * V(3, t) * (G(2,3) * sin( Th(2, t)-Th(3, t) ) - B(2,3) * cos( Th(2, t)-Th(3, t) )) + B(2,3) * V(2, t)^2;
    Q32(t) = V(2, t) * V(3, t) * (G(2,3) * sin( Th(3, t)-Th(2, t) ) - B(2,3) * cos( Th(3, t)-Th(2, t) )) + B(2,3) * V(3, t)^2;
    pL12(t) = abs(P12(t)-P21(t))/2;
    pL13(t) = abs(P13(t)-P31(t))/2;
    pL23(t) = abs(P23(t)-P32(t))/2;
    qL12(t) = abs(Q12(t)-Q21(t))/2;
    qL13(t) = abs(Q13(t)-Q31(t))/2;
    qL23(t) = abs(Q23(t)-Q32(t))/2;
    S12(t) = sqrt(pL12(t)^2 + qL12(t)^2);
    S13(t) = sqrt(pL13(t)^2 + qL13(t)^2);
    S23(t) = sqrt(pL23(t)^2 + qL23(t)^2);
    
     nq(3*t - 2) = S12(t) - Smax;
     nq(3*t - 1) = S13(t) - Smax;
     nq(3*t) = S23(t) - Smax;
    
end
%% END CAPACITY CONSTRAINTS %%%

c = [nq];
end








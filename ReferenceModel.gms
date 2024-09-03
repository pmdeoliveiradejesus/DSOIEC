$title ENERGIES - FIORITTI CASE

Set
   t 'hours'         / t1*t24 /;

Table data(t,*)
    lambda  Pd2u    Pd3u    Ppvu
t1  0.1043  0.9921  0.9914  0.000
t2  0.0957  0.9985  0.9963  0.000
t3  0.0957  0.9985  0.9963  0.000
t4  0.0953  0.9987  0.9963  0.000
t5  0.0941  0.9991  0.9963  0.000
t6  0.0947  0.9989  0.9963  0.000
t7  0.1000  0.9979  0.9963  0.148
t8  0.1196  0.9781  0.9753  0.332
t9  0.1319  0.9689  0.966   0.524
t10 0.1324  0.97    0.9661  0.692
t11 0.1274  0.9758  0.971   0.815
t12 0.1180  0.9835  0.9765  0.917
t13 0.1221  0.9822  0.9765  0.955
t14 0.1217  0.9819  0.9765  0.882
t15 0.1162  0.9851  0.9788  0.834
t16 0.1133  0.9863  0.9805  0.693
t17 0.1101  0.9882  0.9834  0.522
t18 0.1245  0.9727  0.9698  0.333
t19 0.1334  0.9608  0.959   0.140
t20 0.1353  0.9589  0.959   0.000
t21 0.1360  0.9586  0.959   0.000
t22 0.1358  0.9587  0.959   0.000
t23 0.1315  0.9612  0.9604  0.000
t24 0.1149  0.9801  0.9793  0.000;    
* -----------------------------------------------------

Variable
   CW           'Community Welfare $/period'
   
Positive Variables
   Pb2(t)       'Power bought from the public grid 1 kW'
   Ps2(t)       'Power sold to the public grid 1 kW'
   SOC3(t)      'State of Charge 2 kWh'
   Pd3(t)       'Power discharge 2 kW'
   Pch3(t)      'Power Charge 2 kW' 
   Pb3(t)       'Power bought from the public grid 2 kW'
   Ps3(t)       'Power sold to the public grid 2 kW'
   P23CDS(t)    'Power sent from 1 to 2 through a CDS kW'
   P32CDS(t)    'Power sent from 2 to 1 through a CDS kW'
   epsilon(t)   'Internal Local Market Price $/kWh'
   SOCf3        'Final SOC MWh'
   AP2          'Benefit Load 1 + PV'
   AP3          'Benefit Load 2 + BESS'
   APA          'Benefit Aggregator';
         
Binary Variable
   w1(t)     
   w2(t)
   u1(t)     
   u2(t)   
   v12(t) 

Scalar
   Ppvmax2     /55000/
   Pdmax2      /5000/
   Pdmax3      /35000/ 
   Pmax21       /30000/
   Pmax31       /30000/ 
   C3          /100000/
   eff_c3      /1.0/
   eff_d3      /1.0/
   CR_c3       /0.33/
   CR_d3       /0.30/
   psi2        /0.0/
   psi3        /0.0/
   PmaxCDS     /30000/
   varphib /0.0075/
   varphis /0.0025/
   peso1   /1/
   peso2   /1/
   Share3  /0.2/;

SOC3.up(t)     = 1*C3;
SOC3.lo(t)     = 0.2*C3;
Pch3.lo(t) = 0;
Pd3.lo(t) = 0;
Pb2.lo(t) = 0;
Ps2.lo(t) = 0;
Pb3.lo(t) = 0;
Ps3.lo(t) = 0;

Equation CWcalc, balance1, balance2,r1, r2, r3, r4, r5, r6, r7, rcds1, rcds2, rcds3, rcds4, rcds5, rcds6;
CWcalc.. CW =e= peso1*365*(sum((t), (data(t,'lambda')*ps2(t)-(data(t,'lambda')+psi2)*pb2(t)))+
sum((t), (data(t,'lambda')*ps3(t)-(data(t,'lambda')+psi3)*pb3(t)))+
sum((t), (epsilon(t)-varphis)*P23CDS(t))-
sum((t), (epsilon(t)+varphib)*P32CDS(t))+
sum((t), (epsilon(t)-varphis)*P32CDS(t))-
sum((t), (epsilon(t)+varphib)*P23CDS(t)))+
peso2*365*(sum((t), (varphis+varphib)*(P23CDS(t)+P32CDS(t))));
balance1(t)..  Pb2(t) + Ppvmax2*data(t,'Ppvu')+P32CDS(t) =e= Ps2(t) + Pdmax2*data(t,'Pd2u')+P23CDS(t);
balance2(t)..  Pd3(t) + Pb3(t) +P23CDS(t) =e= Pch3(t) +Ps3(t) + Pdmax3*data(t,'Pd3u')+P32CDS(t);
r1(t)..  Pch3(t)  =l= CR_c3*C3*u1(t);
r2(t)..  Pd3(t) =l= CR_d3*C3*(1-u1(t));
r3(t)..  SOC3(t) =e= SOCf3$(ord(t)=1) + SOC3(t-1)$(ord(t)>1) + Pch3(t)*eff_c3  - Pd3(t)/eff_d3;
r4(t)..  Pb2(t) =l= Pmax21*w2(t);
r5(t)..  Ps2(t) =l= Pmax21*(1-w2(t));
r6(t)..  Pb3(t) =l= Pmax31*u2(t);
r7(t)..  Ps3(t) =l= Pmax31*(1-u2(t));
rcds1(t).. P23CDS(t) =l= PmaxCDS*v12(t);
rcds2(t).. P32CDS(t) =l= PmaxCDS*(1-v12(t));
rcds3(t).. epsilon(t) =l= data(t,'lambda')+psi2;
rcds4(t).. epsilon(t) =l= data(t,'lambda')+psi3;
rcds5(t).. epsilon(t) =g= data(t,'lambda');
rcds6..   SOCf3=e=Share3*C3;

Model modelCW / all /;
solve modelCW using MINLP maximizing CW;

execute_unload 'results.gdx';

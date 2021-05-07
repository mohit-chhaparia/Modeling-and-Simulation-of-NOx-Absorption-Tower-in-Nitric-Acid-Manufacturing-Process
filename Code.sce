clc;
clear;
P=input('Enter the absorption column total pressure (in atm) = ');
T=input('Enter the temperature inside the absorption column (in K) = ');
nNO=input('Enter the number of kmoles of nitrogen monoxide in the inlet gas stream = ');
nNO2=input('Enter the number of kmoles of nitrogen dioxide in the inlet gas stream = ');
nN2O3=input('Enter the number of kmoles of dinitrogen trioxide in the inlet gas stream = ');
nN2O4=input('Enter the number of kmoles of dinitrogen tetroxide in the inlet gas stream = ');
nN2=input('Enter the number of kmoles of nitrogen gas in the inlet gas stream = ');
nO2=input('Enter the number of kmoles of oxygen gas in the inlet gas stream = ');
nINERTS=input('Enter the number of kmoles of inert gases in the inlet gas stream (excluding nitrogen gas) = ');
nHNO2=input('Eneter the number of kmoles of nitrous acid gas in the inlet gas stream = ');
nHNO3=input('Enter the number of kmoles of nitric acid gas in the inlet gas stream = ');
nH2O=input('Enter the number of kmoles of H2O vapor in the inlet gas stream = ');
nT=nNO+nNO2+nN2O4+nN2O3+nN2+nO2+nINERTS+nHNO2+nHNO3+nH2O;    //Total number of gaseous kmoles in the inlet gas stream
nTINERTS=nN2+nINERTS;    //Total number of gaseous kmoles of inerts in the inlet gas stream
MONO(1)=nNO;

WHNO2(1)=0;    //Weight fraction of HNO2 in the inlet liquid stream of tray 1
WHNO3(1)=0.6;    //Weight fraction of HNO3 in the inlet liquid stream of tray 1
WH2O(1)=0.4;    //Weight fraction of H2O in the inlet liquid stream of tray 1

X(1)=1;    //Initial guess to find Epsilone 1
Z(1)=1;    //Initial guess to find Epsilone 2
Y(1)=1;    //Initial guess to find Epsilone 3
W(1)=1;    //Initial guess to find Epsilone 4
Q(1)=0.60;
TP(1)=350;
AP=[0:1:50];

for i=1:50
yNO=nNO/nT;    //Mole fraction of NO in the inlet gas stream
yNO2=nNO2/nT;    //Mole fraction of NO2 in the inlet gas stream
yN2O3=nN2O3/nT;    //Mole fraction of N2O3 in the inlet gas stream
yN2O4=nN2O4/nT;    //Mole fraction of N2O4 in the inlet gas stream
yN2=nN2/nT;    //Mole fraction of N2 in the inlet gas stream
yO2=nO2/nT;    //Mole fraction of O2 in the inlet gas stream
yH2O=nH2O/nT;    //Mole fraction of H2O in the inlet gas stream
yINERTS=nINERTS/nT;    //Mole fraction of inerts in the inlet gas stream (excluding nitrogen)
yHNO2=nHNO2/nT;    //Mole fraction of HNO2 in the inlet gas stream
yHNO3=nHNO3/nT;    //Mole fraction of HNO3 in the inlet gas stream
yTINERTS=nTINERTS/nT;    //Mole fraction of total gaseous inerts in the inlet gas stream
pNO=yNO*P;    //Partial pressure of NO in the inlet gas stream (in atm)
pH2O=yH2O*P;    //Partial pressure of H2O vapor in the inlet gas stream (in atm)
pNO2=yNO2*P;    //Partial pressure of NO2 in the inlet gas stream (in atm)
pN2O3=yN2O3*P;    //Partial pressure of N2O3 in the inlet gas stream (in atm)(in atm)
pN2O4=yN2O4*P;    //Partial pressure of N2O4 in the inlet gas stream (in atm)
pN2=yN2*P;    //Partial pressure of N2 in the inlet gas stream (in atm)
pO2=yO2*P;    //Partial pressure of O2 in the inlet gas stream (in atm)
pINERTS=yINERTS*P;    //Partial pressure of INERT gases excluding nitrogen gas in the inlet gas stream (in atm)
pHNO2=yHNO2*P;    //Partial pressure of HNO2 gas in the inlet gas stream (in atm)
pHNO3=yHNO3*P;    //Partial pressure of HNO3 gas in the inlet gas stream (in atm)
pTINERTS=yTINERTS*P;    //Partial pressure of Total INERT gases in the inlet gas stream (in atm)
K(1)=0;    //Garbage
K(2)=10^(2993/T-9.226);    //Equilibrium constant for NO2 to N2O4 (in atm inverse)
K(3)=10^(2072/T-7.234);    //Equilibrium constant for (NO & NO2) to N2O3 (in atm inverse)
K(4)=10^(2051.17/T-6.7328);    //Equilibrium constant for (NO, NO2 & H2O) to HNO2 (in atm inverse)
K(5)=10^(2003.8/T-8.757);    //Equilibrium constant for (NO2 & H2O) to (HNO3 & NO) (in atm inverse)
function ep1=EP1(X)
    ep1 = K(2)*(((nNO2-2*X(1))^2)*P)-(nN2O4+X(1))*(nN2O4+nNO2-X(1));    //Equation to find extent of reaction for NO2 to N2O4
endfunction
function ep2=EP2(Z)
    ep2 = K(3)*((nNO-Z(1))*(nNO2-Z(1))*P)-(nN2O3+Z(1))*(nN2O3+nNO+nNO2-Z(1));    //Equation to find extent of reaction for (NO & NO2) to N2O3
endfunction
function ep3=EP3(Y)
    ep3 = K(4)*((nNO-Y(1))*(nNO2-Y(1))*(nH2O-Y(1))*P)-(((nH2O+Y(1))^2)*(nHNO2+nNO+nNO2+nH2O-2*Y(1)));    //Equation to find extent of reaction for (NO, NO2 & H2O) to HNO2
endfunction
function ep4=EP4(W)
    ep4 = K(5)*(((nNO2-3*W(1))^3)*(nH2O-W(1))*P)-(((nHNO3+2*W(1))^2)*(nNO+W(1))*(nHNO3+nNO+nNO2+nH2O-W(1)));    //Equation to find extent of reaction for (NO2 & H2O) to (HNO3 & NO)
endfunction

ep(1)=fsolve(X(1),EP1);    //Epsilone 1 is extent of reaction for NO2 to N2O4
ep(2)=fsolve(Z(1),EP2);    //Epsilone 2 is extent of reaction for (NO & NO2) to N2O3
ep(3)=fsolve(Y(1),EP3);    //Epsilone 3 is extent of reaction for (NO, NO2 & H2O) to HNO2
ep(4)=fsolve(W(1),EP4);    //Epsilone 4 is extent of reaction for (NO2 & H2O) to (HNO3 & NO)
X(1)=ep(1);
Z(1)=ep(2);
Y(1)=ep(3);
W(1)=ep(4);
YNO=(nNO-ep(2)-ep(3)+ep(4))/(nT-ep(2)-2*ep(3)-ep(4));    //mole fraction of NO gas in the outlet gas stream of tray 1
YNO2=(nNO2-ep(1)-ep(2)-2*ep(3)-ep(4))/(nT-ep(1)-ep(2)-2*ep(3)-ep(4));    //mole fraction of NO2 in the outlet gas stream of tray 1
YN2O3=(nN2O3+ep(2))/(nT-ep(2));    //mole fraction of N2O3 in the outlet gas stream of tray 1
YN2O4=(nN2O4+ep(1))/(nT-ep(1));    //mole fraction of N2O4 in the outlet gas stream of tray 1
YH2O=(nH2O-ep(3)-ep(4))/(nT-ep(3)-ep(4));    //mole fraction of H2O in the outlet gas stream of tray 1
YHNO3=(nHNO3+2*ep(4))/(nT-ep(4));    //mole fraction of HNO3 in the outlet gas stream of tray 1
YO2=(nO2)/(nT);    //mole fraction of O2 in the outlet gas stream of tray 1
YHNO2=(nHNO2+ep(3))/(nT-2*ep(3));    //mole fraction of HNO2 in the outlet gas stream of tray 1
YTINERTS=1-(YNO+YNO2+YN2O3+YN2O4+YH2O+YHNO3+YHNO2+YO2);    //Mole fraction of Total INERTS in the outlet gas stream of tray 1
NINERTS=nINERTS;    //Number of kmoles of INERTS (Excluding N2) in the outlet gas stream of tray 1
NT=nTINERTS/YTINERTS;    //Total number of kmoles in the outlet gas stream of tray 1
YINERTS=NINERTS/NT;    //Mole fraction of INERTS (Excluding N2) in the outlet gas stream of tray 1
NN2=nN2;    //Number of kmoles of N2 in the outlet gas stream of tray 1
NNO=YNO*NT;    //Number of kmoles of NO in the outlet gas stream of tray 1
NNO2=YNO2*NT;    //Number of kmoles of NO2 in the outlet gas stream of tray 1
NN2O3=YN2O3*NT;    //Number of kmoles of N2O3 in the outlet gas stream of tray 1
NN2O4=YN2O4*NT;    //Number of kmoles of N2O4 in the outlet gas stream of tray 1
NO2=YO2*NT;    //Number of kmoles of O2 in the outlet gas stream of tray 1
NH2O=YH2O*NT;    //Number of kmoles of H2O in the outlet gas stream of tray 1
NHNO2=YHNO2*NT;    //Number of kmoles of HNO2 in the outlet gas stream of tray 1
NHNO3=YHNO3*NT;    //Number of kmoles of HNO3 in the outlet gas stream of tray 1
YN2=nN2/NT;    //Mole fraction of N2 in the outlet gas stream of tray 1
NTINERTS=nTINERTS;    //Number of kmoles of total inerts in the outlet gas stream of tray 1
PNO=YNO*P;    //Partial pressure of NO in the outlet gas stream of tray 1 (in atm)
PNO2=YNO2*P;    //Partial pressure of NO2 in the outlet gas stream of tray 1 (in atm)
PN2O3=YN2O3*P;    //Partial pressure of N2O3 in the outlet gas stream of tray 1 (in atm)
PN2O4=YN2O4*P;    //Partial pressure of N2O4 in the outlet gas stream of tray 1 (in atm)
PO2=YO2*P;    //Partial pressure of O2 in the outlet gas stream of tray 1 (in atm)
PH2O=YH2O*P;    //Partial pressure of H2O in the outlet gas stream of tray 1 (in atm)
PHNO2=YHNO2*P;    //Partial pressure of HNO2 in the outlet gas stream of tray 1 (in atm)
PHNO3=YHNO3*P;    //Partial pressure of HNO3 in the outlet gas stream of tray 1 (in atm)
PTINERTS=YTINERTS*P;    //Partial pressure of Total INERTS in the outlet gas stream of tray 1 (in atm)
PN2=YN2*P;    //Partial pressure of N2 in the outlet gas stream of tray 1 (in atm)
PINERTS=YINERTS*P;    //Partial pressure of INERTS (Excluding N2) in the outlet gas stream of tray 1 (in atm)
G3=nT;    //Gaseous flowrate (in kmol/s)
printf('Tray - %d\n',i);
printf('Component                            Inlet stream                                                                 Outlet stream\n');
printf('                Number of moles      Mole fraction       Partial pressure                      Number of moles    Mole fraction    Partial pressure\n');
printf('    NO             %10.8f         %10.8f          %10.8f                            %10.8f         %10.8f        %10.8f\n',nNO,yNO,pNO,NNO,YNO,PNO);
printf('    NO2            %10.8f        %10.8f          %10.8f                            %10.8f        %10.8f        %10.8f\n',nNO2,yNO2,pNO2,NNO2,YNO2,PNO2);
printf('    N2O3           %10.8f         %10.8f          %10.8f                            %10.8f         %10.8f        %10.8f\n',nN2O3,yN2O3,pN2O3,NN2O3,YN2O3,PN2O3);
printf('    N2O4           %10.8f        %10.8f          %10.8f                            %10.8f        %10.8f        %10.8f\n',nN2O4,yN2O4,pN2O4,NN2O4,YN2O4,PN2O4);
printf('    N2             %10.8f      %10.8f          %10.8f                            %10.8f      %10.8f        %10.8f\n',nN2,yN2,pN2,NN2,YN2,PN2);
printf('    O2             %10.8f        %10.8f          %10.8f                            %10.8f        %10.8f        %10.8f\n',nO2,yO2,pO2,NO2,YO2,PO2);
printf('    H2O            %10.8f         %10.8f          %10.8f                            %10.8f         %10.8f        %10.8f\n',nH2O,yH2O,pH2O,NH2O,YH2O,PH2O);
printf('    HNO2           %10.8f         %10.8f          %10.8f                            %10.8f         %10.8f        %10.8f\n',nHNO2,yHNO2,pHNO2,NHNO2,YHNO2,PHNO2);
printf('    HNO3           %10.8f         %10.8f          %10.8f                            %10.8f         %10.8f        %10.8f\n',nHNO3,yHNO3,pHNO3,NHNO3,YHNO3,PHNO3);
printf('    INERTS         %10.8f        %10.8f          %10.8f                            %10.8f        %10.8f        %10.8f\n',nINERTS,yINERTS,pINERTS,NINERTS,YINERTS,PINERTS);
G2=G3*yTINERTS/YTINERTS;
MNO=30; //molecular mass
MNO2=46;
MN2O3=76;
MN2O4=92;
MN2=28;
MO2=32;
MH2O=18;
MHNO2=48;
MHNO3=64;
Mavgi=yNO*MNO + yNO2*MNO2 + yN2O3*MN2O3 + yN2O4*MN2O4 + yN2*MN2 + yO2*MO2 + yH2O*MH2O + yHNO2*MHNO2 + yHNO3*MHNO3;
L2=Mavgi*G3;    //Liquid flowrate (in kg/s)

L1=L2-Mavgi*(abs(G2-G3));
a(1,1)=2*L1/18;    //Matrix calculations to find weight fractions in the inlet liquid stream of tray 1 (Start)
a(1,2)=L1/47;
a(1,3)=L1/63;
a(2,1)=L1/18;
a(2,2)=2*L1/47;
a(2,3)=3*L1/63;
a(3,1)=0;
a(3,2)=L1/47;
a(3,3)=L1/63;
b(1)=L2*(2*WH2O(i)/18+WHNO2(i)/47+WHNO3(i)/63)+G2*(YHNO2+YHNO3+2*YH2O)-G3*(yHNO2+yHNO3+2*yH2O);
b(2)=L2*(WH2O(i)/18+2*WHNO2(i)/47+3*WHNO3(i)/63)+G2*(YNO+2*YNO2+3*YN2O3+4*YN2O4+2*YHNO2+3*YHNO3+YH2O+2*YO2)-G3*(yNO+2*yNO2+3*yN2O3+4*yN2O4+2*yHNO2+3*yHNO3+yH2O+2*yO2);
b(3)=L2*(WHNO2(i)/47+WHNO3(i)/63)+G2*(YNO+YNO2+2*YN2O3+2*YN2O4+YHNO2+YHNO3)-G3*(yNO+yNO2+2*yN2O3+2*yN2O4+yHNO2+yHNO3);
a(3,1)=a(3,1)-a(1,1)*a(3,3)/a(1,3);
a(3,2)=a(3,2)-a(1,2)*a(3,3)/a(1,3);
b(3)=b(3)-b(1)*a(3,3)/a(1,3);
a(2,1)=a(2,1)-a(1,1)*a(2,3)/a(1,3);
a(2,2)=a(2,2)-a(1,2)*a(2,3)/a(1,3);
b(2)=b(2)-b(1)*a(2,3)/a(1,3);
a(3,1)=a(3,1)-a(2,1)*a(3,2)/a(2,2);
b(3)=b(3)-b(2)*a(3,2)/a(2,2);
x1=b(3)/a(3,1);
x2=(b(2)-a(2,1)*x1)/a(2,2);
x3=(b(1)-a(1,1)*x1-a(1,2)*x2)/a(1,3);    //Matrix calculations to find weight fractions in the inlet liquid stream of tray 1 (End)
WHNO2(i+1)=x2;
WHNO3(i+1)=x3;
WH2O(i+1)=x1;
Q(i+1)=x3;
printf('Weight fraction of HNO3 in the inlet liquid stream of tray %d = %10.6f\n',i,x3);
printf('Weight fraction of HNO2 in the inlet liquid stream of tray %d = %10.6f\n',i,x2);
printf('Weight fraction of H2O in the inlet liquid stream of tray %d = %10.6f\n',i,x1);    //Display of weight fractions found using the above Matrix calculations

printf('Liquid mass flow rate (inlet) for tray %d (in g/s) = %10.6f\n',i,L1);
printf('Liquid mass flow rate (Outlet) for tray %d (in g/s) = %10.6f\n',i,L2);
printf('Gas molar flow rate (inlet) for tray %d (in mol/s) = %10.6f\n',i,G3);
printf('Gas molar flow rate (outlet) for tray %d (in mol/s) = %10.6f\n',i,G2);

R=0.08206; //Universal Gas Constant (in (m3 atm)/(kmol K))
V= 16.631; //volume of one tray spacing (oxidiser in m3)
MNO=30; //molecular mass
MNO2=46;
MN2O3=76;
MN2O4=92;
MN2=28;
MO2=32;
MH2O=18;
MHNO2=48;
MHNO3=64;
Mavg=YNO*MNO + YNO2*MNO2 + YN2O3*MN2O3 + YN2O4*MN2O4 + YN2*MN2 + YO2*MO2 + YH2O*MH2O + YHNO2*MHNO2 + YHNO3*MHNO3; //Average molecular weight (kg/kmol) (excluding inerts))
Davg=3.4; //Average density of outlet gas stream from tray 1 (kg/m3)
MFR=G2*Mavg; //Mass Flow Rate (kg/s)
VFR=MFR/Davg; //Volumetric flow rate (m3/s)
t=V/VFR; //Space time (tau) in seconds
k(1)=10^((652.1/T)-0.7356); //Rate constant (atm^(-2)*s^(-1))
CpNO=32.82-(0.02*T)+(5.304*10^(-5)*T^(2))-(3.883*10^(-8)*T^(3))+(9.602*10^(-12)*T^(4));
CpO2=28.31-(2.535*10^(-3)*T)+(2.683*10^(-5)*T^(2))-(2.446*10^(-8)*T^(3))+(6.764*10^(-12)*T^(4));
CpNO2=24.26+(0.04*T)-(1.628*10^(-5)*T^(2))-(6.878*10^(-9)*T^(3))+(3.97*10^(-12)*T^(4));
CpN2O4=33.054+(18.661*10^(-2)*T)-(11.339*10^(-5)*T^(2));
CpN2O3=32.301+(17.961*10^(-2)*T)-(16.41*10^(-5)*T^(2));
CpHNO2=24.89974+(91.37563*10^(-3)*T)-(64.84614*10^(-6)*T^(2))+(17.92007*10^(-9)*T^(3));
CpHNO3=19.63229+(153.9599*10^(-3)*T)-(115.8378*10^(-6)*T^(2))+(32.87955*10^(-9)*T^(3));
CpH2O=34.51-(0.01*T)+(4.686*10^(-5)*T^(2))-(3.781*10^(-8)*T^(3))+(1.173*10^(-11)*T^(4));
CpAvgi=(CpNO*yNO)+(CpO2*yO2)+(CpNO2*yNO2)+(CpN2O4*yN2O4)+(CpN2O3*yN2O3)+(CpHNO2*yHNO2)+(CpHNO3*yHNO3)+(CpH2O*yH2O);
CpAvgo=(CpNO*YNO)+(CpO2*YO2)+(CpNO2*YNO2)+(CpN2O4*YN2O4)+(CpN2O3*YN2O3)+(CpHNO2*YHNO2)+(CpHNO3*YHNO3)+(CpH2O*YH2O);
CNOi=PNO/(R*T); //Concentration of NO entering the oxidizer (kmol/m3)
CNOo=CNOi/2; //Concentration of NO leaving the oxidizer (initial guess for fsolve) (in kmol/m3)
function abcd=ABCD(CNOo)
    abcd = k*t*(CNOo^2)+CNOo-CNOi;
endfunction
CNOo=fsolve(CNOo,ABCD); //Concentration of O2 leaving the oxidizer (in kmol/m3)
PNOo=CNOo*R*T; //Partial pressure of NO in the outlet stream from oxidizer (in atm)
YNOo=PNOo/P; //Mole fraction of NO in the outlet stream from oxidizer
NNOo=YNOo*NT; //Number of moles of NO in the outlet stream from the oxidizer
CO2i=PO2/(R*T); //Concentration of O2 entering the oxidizer (kmol/m3)
CO2o=CO2i-0.5*(CNOi-CNOo); // Concentraion of O2 leaving the oxidizer (in kmol/m3)
PO2o=CO2o*R*T; //Partial pressure of O2 in the outlet stream from oxidizer (in atm)
YO2o=PO2o/P; //Mole fraction of O2 in the outlet stream from oxidizer
NO2o=YO2o*NT; //Number of moles of O2 in the outlet stream from the oxidizer
CNO2i=PNO2/(R*T); //Concentration of NO2 entering the oxidizer (kmol/m3)
CNO2o=(CNOi-CNOo)+CNO2i; //Concentraion of NO2 leaving the oxidizer (in kmol/m3)
PNO2o=CNO2o*R*T; //Partial pressure of NO2 in the outlet stream from oxidizer (in atm)
YNO2o=PNO2o/P; //Mole fraction of NO2 in the outlet stream from oxidizer
NNO2o=YNO2o*NT; //Number of moles of NO2 in the outlet stream from the oxidizer
NTo=NNOo+NNO2o+NN2O3+NN2O4+NN2+NO2o+NH2O+NHNO2+NHNO3+NINERTS; //Total number of moles in the outlet stream from oxidizer
N2O3f=G2*(YN2O3-yN2O3);
dH1=-1573.353 - 24.779*T + 119.805*10^-3*T^2 - 66.953*10^-6*T^3 + 11.427*10^-9*T^4;
Q1=N2O3f*dH1;
N2O4f=G2*(YN2O4-yN2O4);
dH2=1561.237 - 15.466*T + 53.505*10^-3*T^2 - 64.74*10^-6*T^3 + 3.439*10^-9*T^4;
Q2=N2O4f*dH2;
HNO3f=G2*(YHNO3-yHNO3);
dH3=17606.292 - 74.25*T + 0.045*T^2 + 18.34*10^-6*T^3 + 4.904*10^-9*T^4;
Q3=HNO3f*dH3;
HNO2f=G2*(YHNO2-yHNO2);
dH4=6469.15 - 41.812*T + 0.08637*T^2 - 71.104*10^-6*T^3 + 29.839*10^-9*T^4;
Q4=HNO2f*dH4;
NOf=G2*(YNO-YNOo-YN2O3+yN2O3-0.5*YHNO2+0.5*yHNO2+0.5*YHNO3-0.5*yHNO3);
dH5=12924.19 - 45.43*T + 21.267*10^-3*T^2 - 55.156*10^-6*T^3 + 22.091*10^-9*T^4;
Q5=NOf*dH5;
QT=Q1+Q2+Q3+Q4+Q5;
Mavgi=yNO*MNO + yNO2*MNO2 + yN2O3*MN2O3 + yN2O4*MN2O4 + yN2*MN2 + yO2*MO2 + yH2O*MH2O + yHNO2*MHNO2 + yHNO3*MHNO3;
Mavgo=YNOo*MNO + YNO2o*MNO2 + YN2O3*MN2O3 + YN2O4*MN2O4 + YN2*MN2 + YO2o*MO2 + YH2O*MH2O + YHNO2*MHNO2 + YHNO3*MHNO3;
mi=G2*Mavgi;
me=NTo*Mavgo;
Te=(QT+(mi*CpAvgi*T))/(me*CpAvgo);
printf('Total heat evolved (in KJ/(Kmol).(K)) on tray %d = %10.6f\n',i,QT);
printf('Temperature of the gas stream (in K) before reaching the cooling coil after tray %d = %10.6f\n',i,Te);
TP(i+1)=Te;

nNO=NNOo;
MONO(i+1)=nNO;
nNO2=NNO2o;
nN2O3=NN2O3;
nN2O4=NN2O4;
nN2=NN2;
nO2=NO2o;
nINERTS=NINERTS;
nHNO2=0;
nHNO3=0;
nH2O=NH2O;
nT=nNO+nNO2+nN2O4+nN2O3+nN2+nO2+nINERTS+nHNO2+nHNO3+nH2O;    //Total number of gaseous kmoles in the inlet gas stream
nTINERTS=nN2+nINERTS;    //Total number of gaseous kmoles of inerts in the inlet gas stream
end
figure(1);
plot(AP,Q);
figure(2);
plot(AP,MONO);
figure(3);
plot(AP,TP);

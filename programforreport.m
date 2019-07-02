%initial parameters
SF=0.8;
ETI=0.85;
ETS=0.72;
B1s=-60;
Nu=0.45;
al2=65;
r1sr2=0.65;
g=1.4;
prompt = 'What is the pressure ratio? ';
PR = input(prompt);
prompt = 'what is the blade backsweep angle? ';
B2i = input(prompt);
prompt = 'what is the initial Pressure? ';
P01 = input(prompt);
prompt = 'what is the initial Temperature? ';
T01 = input(prompt);
D01=P01/(287*T01);
X = sprintf('Stage1\nPressure- %d\nTemperature- %d\nDensity- %d',P01,T01,D01);
disp(X);
B1=torad(B1s);
B2=torad(B2i);
a2=torad(al2);
AL=SF/(1-tan(B2)/tan(a2));                      %work factor
U2A01=((PR^((g-1)/g)-1)/(0.4*ETS*AL))^0.5;      %Mach no from pressure ratio for compressible flow
CT2A01=AL*U2A01;                                %work factor x Mu
T02T01=1+(g-1)*CT2A01*U2A01;
P02P01=(ETI*(T02T01-1)+1)^(g/(g-1));
P02=P02P01*P01;
T02=T02T01*T01;
D02=P02/(287*T02);
X = sprintf('Stage2\nPressure- %d\nTemperature- %d\nDensity- %d',P02,T02,D02); %Props at 2 
disp(X);
%% 

%outlet velocities
C2A01=CT2A01/sin(a2);
C2A02=C2A01/T02T01^0.5;
T2T02=1-(g-1)*C2A02^2/2;
AM2=C2A02/T2T02^0.5;                            %Discharge mach no (can be slightly more than 1)
U2A02=U2A01/T02T01^0.5;
U2A2=U2A02/T2T02^0.5;
CT2A2=CT2A01/(T02T01*T2T02)^0.5;
WT2A2=CT2A2-U2A2;
CM2A2=AM2*cos(a2);
AM2R=(WT2A2^2+CM2A2^2)^0.5;
B2=atand(WT2A2/CM2A2);
CM2A01=CT2A01/tan(a2);
%% 

%inlet velocities
U1A01=U2A01*r1sr2;
C1A01=-U1A01/tan(B1);
AM1=C1A01*(1-((g-1)*C1A01^2)/2)^-0.5;
AM1R=AM1/cos(B1);
MR=AM2R/AM1R;                                     
W1SA01=(U1A01^2+C1A01^2)^0.5;
W2A01=AM2R*(T02T01*T2T02)^0.5;
DR=W2A01/W1SA01;                                %reaches max 1.4 for subsonic inducer
                                                %DR is a measure of diffusion in impeller
%% 
                                               
%PERFORMANCE PARAMETERS
D1D01=(1+((g-1)*AM1^2)/2)^(-1/(g-1));
P2P02=T2T02^(g/(g-1));
D2D01=P2P02*P02P01/(T2T02*T02T01);
PSI=2*AL*ETS;
PHI=D1D01*r1sr2^2*(1-Nu^2)*C1A01/U2A01;
THETA=PHI*U2A01;                                %ratio of actual mass flow to that through orifice at a01
B2R2=THETA/(2*D2D01*CM2A01);                   
SS=(pi*PHI)^0.5/(PSI/2)^0.75;                   
SSG=(pi*PHI/D1D01)^0.5/(PSI/2)^0.75;
WND=PSI*THETA*(U2A01)^2;
%% 

%GEOMETRY VARIABLES
prompt = 'what is the mass flow rate? ';        
MF = input(prompt);
A01=(g*287*T01)^0.5;
A02=(g*287*T02)^0.5;
LD2=(0.28*(AM1R+0.8)*(1-0.7754*r1sr2)*r1sr2*(1-Nu))^0.5;              %Nu=0.45
A2=MF/(THETA*D01*A01);
r2=(A2/pi)^0.5;
r1s=r1sr2*r2;
b2=B2R2*r2;
r1h=Nu*r1s;
L=LD2*2*r2;

%% inlet velocity triangle
figure(1);
o = [0 0];  %# Origin
va2 = [U1A01 0];  %# Vector 1
vb2 = [0 C1A01];  %# Vector 2
vc2 = va2+vb2;      %# Resultant
mag=norm(vc2);      %magnitude of resultant
arrowStarts = [o; va2; o];        %# Starting points for arrows
arrowEnds = [va2; vc2; vc2];          %# Ending points for arrows
arrow(arrowStarts,arrowEnds);   %# Plot arrows
axis equal;
%% outlet velocity triangle
figure;
o = [0 0];  %# Origin
va2 = [CT2A2 0];  %# Vector 1
vb2 = [0 CM2A2];  %# Vector 2 
vc2 = va2+vb2;      %# Resultant
vd2 = [U2A2 0];
arrowStarts = [o; va2; o; vd2; vd2];        %# Starting points for arrows
arrowEnds = [va2; vc2; vc2; va2; vc2];          %# Ending points for arrows
arrow(arrowStarts,arrowEnds);   %# Plot arrows
axis equal;
%% iteration with new formulae
AM1R=(AM1^2+U2A01^2*THETA*(1+(g-1)*AM1^2/2)^((3*g-1)/2*(g-1))/((1-Nu^2)*AM1))^0.5;
N=U2A01*60*A01/(2*pi*r2);



%close all;
%clear; clc;

P0 = 50000;
TPool_R = 20;
T_In = TPool_R;

% - Geometric data --------------------------------------------------------
N_Fuel = 83;                                                     % - Number of fuel elements
N_Channel = 91;                                               % - Number of channels in the core
V_Fuel = 3.5825e-4;                                         % - Fuel volume [m^3]
R_FuelO = 0.01791;                                       % - Fuel outer radius [m]
RCladO = 0.0187;                                           % - Cladding outer radius [m]
R_CladI = 0.01794;                                          % - Cladding inner radius [m]
H_Core = 0.7224;                                             % - Core height [m]
H_Active = 0.381;                                            % - Core active height [m]
F_D = 0.03;                                                       % - Darcy friction coefficient
AFlow = 0.06352;                                           % - Coolant core flow area [m^2]
W = 0.5;

R_GapO = 0.5*(R_FuelO+R_CladI);                % - Gap radius [m]
PWet = 2*pi*RCladO*91;                             % - Wetted Perimeter [m]
DEq = 4*AFlow/PWet;                               % - Core equivalent diameter [m]

% - Neutronic data --------------------------------------------------------
BetaTOT = 730*10^-5;                                      % - Overall delayed neutron fraction
Life = 60*10^-6;                                                % - Mean neutron generation time [s]

Lambda = [3.01 1.14 0.301 0.111 0.0305 0.0124];
Fraction = [0.042 0.115 0.395 0.196 0.219 0.033];
Beta = BetaTOT*Fraction;

% - Poison data -----------------------------------------------------------
Y_I = 0.0639;                                                      % - Iodine fission yield
Y_Xe = 0.00237;                                                   % - Xenon fission yield
Y_Pm = 0.01071;                                                   % - Promethium fission yield
LambdaI = 2.87*10^-5;                                      % - Iodine decay constant [1/s]
LambdaXe = 2.09*10^-5;                                  % - Xenon decay constant [1/s]
LambdaPm = 3.63*10^-6;                                  % - Promethium decay constant [1/s]
SigmaXe = 2.65*10^-18;                                   % - Xenon absorption cross section [cm^2]
SigmaSm = 4.1*10^-20;                                     % - Samarium absorption cross section [cm^2]
SigmaU = 680.8*10^-24;                                   % - Uranium absorption cross section [cm^2]
SigmaFiss = 580.2*10^-24;                               % - Fuel fission cross section [cm^2]
U_0 = 9.813*10^19;                                     % - Fuel atom volumetric density [1/cm^3]
N_Mean = 2.43;                                                 % - Average neutrons per fission of uranium

EFiss = (200*10^6)*(1.6*10^-19);
ESigmaF = U_0*SigmaFiss;
Phi_0 = P0/(EFiss*ESigmaF*N_Fuel*V_Fuel*10^6);
I_0 = (Y_I*SigmaFiss*U_0*Phi_0)/LambdaI;
Xe_0 = (Y_Xe*SigmaFiss*U_0*Phi_0 + LambdaI*I_0)/(SigmaXe*Phi_0 + LambdaXe);
Pm_0 = (Y_Pm*SigmaFiss*U_0*Phi_0)/LambdaPm;
Sm_0 = (LambdaPm*Pm_0)/(SigmaSm*Phi_0);

% - Thermal-hydraulic data ------------------------------------------------
f = 1;                                                                  % - Fraction of heat generated within the fuel
K_Fuel = 17.58;                                                 % - Fuel thermal conductivity [W/m?C]
H_Gap = 2550;                                                  % - Gap heat transfer coefficient [W/m^2?C]
K_Clad = 228;                                                   % - Cladding thermal conductivity [W/m?C]
Pressure = 1.5;                                                   % - Reactor pressure [bar]
CCool = XSteam('Cp_pT', Pressure, T_In);
CCool = CCool*1000;
Rho_In = XSteam('rho_pT', Pressure, T_In);
Nu_Visc = XSteam('my_pT',Pressure, T_In)/Rho_In;

NuExp = 245*10^-6;                                          % - Coolant thermal expansion coefficient [1/?C]
MCool = 22.57;                                               % - Coolant total mass [kg]

% - Reactor pool data -----------------------------------------------------
MPool = 18963;                                              % - Mass of water in the pool [kg]
M_Res = 800;                                                    % - Mass of water in the reservoir [kg]
M1 = 7;                                                            % - Primary flow rate in the secondary pool [kg/s]
M3 = 7.6;                                                         % - Mass flow rate of the acqueduct water [kg/s]
T_Aq = 15;                                                      % - Average aqueduct water temperature [?C]

% - Evaluation of the stationary mass flow rate ---------------------------
A_Frict = (H_Core*F_D)./(2*DEq*(AFlow^2)*Rho_In);
GCore_0 = (H_Core*9.81*Rho_In*NuExp.*P0./(A_Frict.*CCool)).^(1/3);
U_Cool0 = GCore_0./(AFlow*Rho_In);

% - Reactivity Feedback ---------------------------------------------------
Alpha_F = -8.5*10^-5;                                      % - Fuel temperature feedback coefficient [1/?C]
Alpha_C = 4*10^-5;                                          % - Coolant temperature feedback coefficient [1/?C]
Alpha_V = -2.5e-06;                                         % - Coolant void feedback coefficient [1/cm^3]

% - Stationary quantities (coolant) ---------------------------------------------
TPool_0 = TPool_R+(P0./(MPool*CCool));
TCool_0 = TPool_0+(P0./(20*GCore_0.*CCool));

MuI = XSteam('my_pT', Pressure, T_In);
Rho_0 = XSteam('rho_pT', Pressure, TCool_0);
Mu = XSteam('my_pT', Pressure, TCool_0);
C_Cool0 = XSteam('Cp_pT', Pressure, TCool_0);
C_Cool0 = C_Cool0*1000;
KCool = XSteam('tc_pT',Pressure,TCool_0);

KPool = M1*C_Cool0;

% - Evaluation of the coolant properties at mean temperature --------------
Re_Cool0 = GCore_0.*DEq./(Mu*AFlow);
Pr = (Mu.*CCool)./KCool;
Nu_Cool = 0.023*(Re_Cool0.^0.8).*(Pr.^0.4);
H_Cool0 = (KCool.*Nu_Cool./DEq);

% - Evaluation of the overall heat transfer coefficient -------------------
RT_Fuel = 1/(8*pi*K_Fuel);                                
RT_Gap = 1/(2*pi*R_GapO*H_Gap);                              
RT_Clad = (log(RCladO/R_CladI))/(2*pi*K_Clad);   
RT_Cool = 1./(2*pi*RCladO*H_Cool0);           
RT_Tot = RT_Fuel+RT_Gap+RT_Clad+RT_Cool;
K = (N_Fuel*H_Active)./(RT_Tot)/2;    

% - Evaluation of the reactor time constants ------------------------------
TFuel_0 = (f*P0./K)+TCool_0;
C_Fuel = (750+1.55*(TFuel_0-25));
TauF = C_Fuel./K;
TauC = MCool.*CCool./K;

% - Input Reactivity Block ------------------------------------------------
R_Rod = 0.00089;
Step = 1;

% A = [R_Rod-BetaTOT/Life Beta/Life Alpha_F/Life Alpha_C/Life 0 0 0 -SigmaXe/(ESigmaF*N_Mean*Life) 0 -SigmaSm/(ESigmaF*N_Mean*Life) 0 0
%     Lambda' -diag(Lambda) zeros(6,10)
%     P0/(TauF*K) zeros(1,6) -1/TauF 1/TauF 0 0 0 0 0 0 0 0
%     0 zeros(1,6) 1/TauC -(1/TauC+2*CCool*GCore_0/(TauC*K)) 2*CCool*(TCool_0-TPool_0)/(TauC*K) 0 0 0 0 0 2*CCool*GCore_0/(TauC*K) 0
%     0 zeros(1,6) 0 -2*Rho_In*Nu_Visc*9.81*AFlow -2*GCore_0*AFlow*A_Frict/H_Core 0 0 0 0 0 2*Rho_In*Nu_Visc*9.81*AFlow 0
%     -SigmaU*U_0 zeros(1,6) 0 0 0 -SigmaU*Phi_0 0 0 0 0 0 0
%      Y_I*SigmaFiss*Phi_0 zeros(1,6) 0 0 0 Y_I*SigmaFiss*Phi_0 -LambdaI 0 0 0 0 0
%      (Y_Xe*SigmaFiss*Phi_0-SigmaXe*Phi_0) zeros(1,6) 0 0 0 Y_Xe*SigmaFiss*Phi_0 LambdaI  -LambdaXe-SigmaXe*Phi_0 0 0 0 0
%      Y_Pm*SigmaFiss*Phi_0 zeros(1,6) 0 0 0 Y_Pm*SigmaFiss*Phi_0  0 0 -LambdaPm 0 0 0
%      -SigmaSm*Phi_0 zeros(1,6) 0 0 0 0 0 0 LambdaPm -SigmaSm*Phi_0 0 0
%      0 zeros(1,6) 0 GCore_0/MPool (TCool_0-TPool_0)/MPool 0 0 0 0 0 (GCore_0+M1)/MPool M1/MPool
%      0 zeros(1,6) 0 0 0 0 0 0 0 0 M1/M_Res (M1-M3)/M_Res];

disp('Starting Simulink simulation');tic;
sim('TRIGA_model.slx');toc;
disp('Simulink run ended');
[A,B,C,D]=linmod('TRIGA_model');

% save('ASimulink2.mat','A');
% load('ASimulink2.mat');

% A_New = A;
% A_New(1,1) = -0.833;
% 

Y0 = [1,1,1,1,1,1,1,TFuel_0,TCool_0,GCore_0,U_0,I_0,Xe_0,Pm_0,Sm_0,TPool_0,TPool_R];
T_Span = [0,1800];

Y = ode15s(@(t,Y2) myode(t,Y2,A),T_Span,Y0);

plot(Y.x,Y.y(1,:)*P0)

% plot(Y.x,Y.y(1,:),'b-'); hold on;
% plot(Y2.x,Y2.y(1,:),'r-');
% plot(CoolTData(:,1),CoolTData(:,2)+273.15,'g-');
% legend('TCool-Analytic','TCool-Simu');
% 

%Y_S = ode15s(@(t,Y) myode(t,Y,A_AnS),T_Span,Y0);


%sim('TRIGA_model_basic.slx');

% Rho_Cool = Alpha_C*(CoolTData(1:end,2)-TCool_0);
% Rho_Fuel = Alpha_F*(FuelTData(1:end,2)-TFuel_0);
% Rho_Xe = (SigmaXe+SigmaSm)*(ScopePoison(1:end,2)*10^18)/(2.43*ESigmaF);
% Rho_Rod = ReactivityData(:,2) - Rho_Cool - Rho_Fuel - Rho_Xe;


% reduce_plot(Rho_Cool);hold on;
% PointFinder = findobj(gca,'Type','line');
% Stability_Red = (transpose([get(PointFinder,'Xdata');get(PointFinder,'Ydata')]));
% Stability = table(Stability_Red(:,1),Stability_Red(:,2));
% writetable(Stability,'RhoCool50.txt','Delimiter',',');
% type 'RhoCool50.txt';
% 
% reduce_plot(Rho_Fuel);hold on;
% PointFinder = findobj(gca,'Type','line');
% Stability_Red = (transpose([get(PointFinder,'Xdata');get(PointFinder,'Ydata')]));
% Stability = table(Stability_Red(:,1),Stability_Red(:,2));
% writetable(Stability,'RhoFuel50.txt','Delimiter',',');
% type 'RhoFuel50.txt';
% 
% reduce_plot(Rho_Xe);hold on;
% PointFinder = findobj(gca,'Type','line');
% Stability_Red = (transpose([get(PointFinder,'Xdata');get(PointFinder,'Ydata')]));
% Stability = table(Stability_Red(:,1),Stability_Red(:,2));
% writetable(Stability,'RhoXe50.txt','Delimiter',',');
% type 'RhoXe50.txt';
% 
% reduce_plot(Rho_Rod);hold on;
% PointFinder = findobj(gca,'Type','line');
% Stability_Red = (transpose([get(PointFinder,'Xdata');get(PointFinder,'Ydata')]));
% Stability = table(Stability_Red(:,1),Stability_Red(:,2));
% writetable(Stability,'RhoRod50.txt','Delimiter',',');
% type 'RhoRod50.txt';
% 
% reduce_plot(ReactivityData(:,2));hold on;
% PointFinder = findobj(gca,'Type','line');
% Stability_Red = (transpose([get(PointFinder,'Xdata');get(PointFinder,'Ydata')]));
% Stability = table(Stability_Red(:,1),Stability_Red(:,2));
% writetable(Stability,'RhoTOT50.txt','Delimiter',',');
% type 'RhoTOT50.txt';plot(Rho_Rod);hold on;

% Topic: FEMM AND ANALYTICAL MODELS OF ELECROMAGNETIC CHARACTERISATION
% Name: Omigbodun Tirenioluwa Lois
% EENG20002
%%
clear all;
addpath(genpath('C:\femm42')); % calls the path of the actuator file
openfemm() % opens the femm file in the window 
%% Loading Coordinates and Opening a New Document 
opendocument('femm_template.FEM'); % opens a femm temeplate document ensuring that it is a nagnetic file 
mi_saveas('actuator.fem') ;% changes the name of the femmm templlate to chosen name 
load('corep.mat'),load('moverp.mat'),load('coil1p.mat'),load('coil2p.mat'),load('coil3p.mat'),load('coil4p.mat')
list_of_variable={corep;moverp;coil1p;coil2p;coil3p;coil4p};

%% Block and Segment Addition, Labelling and Properties 
for y = 1: length(list_of_variable)
   var = list_of_variable{y};
   mi_addnode(var(:,1),var(:,2))

   for i = 1:length(var) % length( array) shows the length of the item, is array a x b , output will be a 
     % Drawing and selecting nodes
     mi_selectnode(var(i,1),var(i,2))
     mi_setnodeprop('<None>',y)
     mi_clearselected()

     if i < length(var) 
        mi_addsegment(var(i,1),var(i,2),var((i+1),1),var((i+1),2));   
        mi_selectsegment((1/2*(var(i,1)+ var(i+1,1))), (1/2 *(var(i,2)+ var(i+1,2))));
        mi_setsegmentprop('<none>',1,0,0,y)
        mi_clearselected()
     end
     if i == length(var) 
       mi_addsegment(var(length(var) ,1),var(length(var) ,2),var((1),1),var((1),2));    
       mi_selectsegment((1/2*(var(end,1)+ var(1,1))), (1/2 *(var(1,2)+ var(end,2))));
       mi_setsegmentprop('<None>',1,0,0,y)
       mi_clearselected()
     end                 
   end
    % Block Label
    if y == 1 
       mi_addblocklabel(10,9);
      mi_selectlabel(10,9);
      mi_setblockprop('core_linear',0,0,'<None>',0,1,0)
      mi_clearselected()
    elseif y == 2 
        mi_addblocklabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
        mi_selectlabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
        mi_setblockprop('core_linear',0,0,'<None>',0,2,0)
        mi_clearselected()
    elseif y == 3 
      mi_addblocklabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
      mi_selectlabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
      mi_setblockprop('copper',0,0,'winding_1',0,3,100)
      mi_clearselected()        
        
    elseif y == 4 
        mi_addblocklabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
        mi_selectlabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
        mi_setblockprop('copper',0,0,'winding_1',0,4,-100)
        mi_clearselected()
    elseif y == 5 
       mi_addblocklabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
       mi_selectlabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
       mi_setblockprop('copper',0,0,'winding_2',0,5,100)
       mi_clearselected()        
    elseif y == 6
       mi_addblocklabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
       mi_selectlabel((1/2*(var(1,1)+ var(3,1))), (1/2 *(var(1,2)+ var(3,2))));
       mi_setblockprop('copper',0,0,'winding_2',0,6,-100)
       mi_clearselected()
    end 
end  

 %% Mesh refinement for the Air-gap
   % line segment 1
   mi_selectsegment((1/2*(corep(10,1)+ corep(11,1))), (1/2 *(corep(10,2)+ corep(11,2)))); 
   mi_setsegmentprop('<none>',0.5,0,0,1)
   mi_clearselected()
   % line segment 2
   mi_selectsegment((1/2*(moverp(4,1)+ moverp(1,1))), (1/2 *(moverp(4,2)+ moverp(1,2)))); 
   mi_setsegmentprop('<none>',0.5,0,0,2)
   mi_clearselected()

%% Circular shells for Geometrically boundary
% Using the coordinates shown in the femm interface 
mi_makeABC()
% Problem Definition and Excitiation
mi_probdef(0,'millimeters','planar', 1E-8, 20,30,0) 
% Block of outer shell
 mi_addblocklabel(-6,9);
 mi_selectlabel(-6,9);
 mi_setblockprop('air',0,0,'<None>',0,7,0)
 mi_clearselected()
 
%% VALUES 
%---------MESHING THE MODEL--------------
% Smart Mesh Control
% smartmesh(1)
% mi_createmesh()
% mi_showmesh()
%-------Analyse the Model----------------- 
% mi_analyse()
% mi_loadsolution()

%% ---Circuit properties Outputs ---------------
linear_test = output(1);
nonlinear_test = output(2);

%% Measurement of Resistance from the 2-D FOR ONLY A WINDING    
for num = 1:4 
mo_groupselectblock((num+2));
Volu{num}= mo_blockintegral(10);
mo_clearblock()
end
% Active Length of winding (L) and Resitance of the windings (R)
Volume_2D ={(Volu{1}+Volu{2}),(Volu{3}+Volu{4})};
for num1 = 1:2 
mo_groupselectblock((num1+2));
Area{num1} = mo_blockintegral(5);
mo_clearblock();
L_length{num1} = division(Area{num1},Volume_2D{num1});
Resistance_2D{num1} = resitance(100,L_length{num1},58e06,0.6,Area{num1});% resistance for 2-d model 
end
%%  Measurement of Resistance from the 3-D FOR ONLY A WINDING
Volume_3D = 81153e-9; % Volume of a Winding
Area_3D = 49e-3 * 13.2e-3;
Length_3D = division(Area_3D,Volume_3D);
Resistance_3D  = resitance(100,Length_3D,58e06,0.6,Area_3D);% resistance for 3-d model 
%% POWER VALUES
Current = 0:1:10;
Power_Femm = linear_test{5};% Power values for FEMM Linear 
Power_3D_model = Resistance_3D * Current.^2;% Power values for 3d model 
Power_Femm_ana = Resistance_2D{1} * Current.^2;% Power values for FEMM Analytical Linear 
%% Inductance Values 
air_gap = [0.005, 0.004, 0.003, 0.002, 0.001, 0.0001,];
L_Inductance_core_fringing= core_fringing(1) ;% Indcutance values for MEC with fringing 
L_Inductance_core_no_fringing= core_fringing(2) ;% Indcutance values for MEC without fringing 
inductance_linear = linear_test{4};% Indcutance values for linear Femm 
induct_non_linear = nonlinear_test{4};% Indcutance values for non-linear Femm 
%%%FOR circuit Model
N=100;
I=10;
BHcurve= xlsread('B-Hcurve.xlsx');
plot(BHcurve(:,2),BHcurve(:,1));

g_v = [0.005, 0.004, 0.003, 0.002, 0.001,0.0001];
Refffringe = core_fringing_reluc(1);
Reffnofringe = core_fringing_reluc(2) ;
Ac = 400e-6;
uo = 4*pi*1e-7;
ur = 1000;

 for m =1:length(Refffringe) 
     set_param('BHMODEL/Reffective','R','Refffringe(m)')
     simOut = sim('BHMODEL');
     flux_fringe(m)= simOut.flux(1);
     L_fringe(m)= (N*flux_fringe(m))/I;
     B_in (m)=  flux_fringe(m)/Ac ;
     H_in(m) =  B_in (m)/(uo*ur);
end
 for m =1:length(Reffnofringe) 
     set_param('BHMODEL/Reffective','R','Reffnofringe(m)')
    simOut = sim('BHMODEL');
     flux_nofringe(m)= simOut.flux(1);
     L2_nofringe(m)= (N*flux_nofringe(m))/I;
      B2_in(m) =  flux_nofringe(m)/Ac ;
      H2_in(m) =  B2_in(m)/(uo*ur);
 end
figure
plot(BHcurve(:,2),BHcurve(:,1),'k-');
hold on
plot(H2_in,B2_in,'rd')
plot(H_in,B_in,'bs')
hold off


%% PSI Values
Psi_fringing = Psi(L_Inductance_core_fringing);%  for MEC with fringing 
Psi_nofringing =Psi(L_Inductance_core_no_fringing);% for MEC without fringing  
non_linear_psi = nonlinear_test{1}; %  the nonlinear_test output list 
linear_psi = linear_test{1};% the linear_test output list 
%% Co energy Values
Co_energy_linear = 2*Integral(linear_psi);% the linear_test output list
Co_energy_nonlinear = 2*Integral(non_linear_psi);% the nonlinear_test PSI-I METHOD output list 
Co_energy_fring = 2*Integralside(Psi_fringing);% for MEC with fringing 
Co_energy_nofring = 2*Integralside(Psi_nofringing);% for MEC without fringing 
Nonlinear_Femm_coenergy = Coenergy(nonlinear_test{6}); % the FEA MODEL nonlinear_test output list
Armature_distance = [0,0.001,0.002,0.003,0.004,0.0049];
%% Force Values 
Force_Fring = Gradient(Co_energy_fring);%  for MEC with fringing 
Force_No_Fring = Gradient(Co_energy_nofring);%  for MEC without fringing 
Force_linear= Gradient(Co_energy_linear);% the linear_test output list
Force_nonlinear= Gradient(Co_energy_nonlinear);% the nonlinear_test output list
mid_points =[0,0.0010,0.0015,0.0025,0.0035,0.00445];

%% Plots
% Plot for the Inductance
h1=figure;
plot(air_gap,L_Inductance_core_fringing,'k-s','LineWidth',1)
hold on 
plot(air_gap,L_Inductance_core_no_fringing,'k-o','LineWidth',1)
plot(air_gap,induct_non_linear,'k--d','LineWidth',1)
plot(air_gap,inductance_linear,'k--x','LineWidth',1)
plot(air_gap,L_fringe,'rd','LineWidth',1)
plot(air_gap,L2_nofringe,'ro','LineWidth',1)
legend('MEC with Air-gap Fringing','MEC without Air-gap Fringing','FEMM:Linear ','FEMM:Non Linear','SIMSCAPE MEC with Air-gap Fringing','SIMSCAPE MEC without Air-gap Fringing')
%set(gca,'XGrid','on','YGrid','on')
xlabel('Air Gap Distance (m)');ylabel('Inductance (H)')
title('Inductance to Air Gap Distance')
hold off 
ax = gca;
ax.XAxis.Exponent =0;
ax.YAxis.Exponent =0;
print(h1,'-djpeg','-r200','Inductance1')
%%
h10 = figure;
plot(BHcurve(:,2),BHcurve(:,1),'k-','LineWidth',1.5);
hold on
plot(H2_in,B2_in,'rd','LineWidth',1.5)
plot(H_in,B_in,'bs','LineWidth',1.5)
hold off
legend('FEMM','SIMSCAPE MEC with Air-gap Fringing','SIMSCAPE MEC without Air-gap Fringing','Location','east')
xlabel('Magnetic Intensity H(A/m)');ylabel('Magnetic Density B (T)')
title('B-H Curve')
hold off 
print(h10,'-djpeg','-r300','B-Hcurve')
%%
% Plot for the Power Loss 
h7= figure;
plot(Current,Power_Femm,'k-x','LineWidth',1)
hold on
plot(Current,Power_3D_model,'k-d','LineWidth',1)
plot(Current,Power_Femm_ana,'k-o','LineWidth',1)
xlabel('Current (A)');ylabel('Power (W)')
title('Power Loss to Current')
legend('FEA Model','Analytical 3-d model','Analytical FEMM model','Location','west')
hold off 
ax = gca;
ax.XAxis.Exponent =0;
ax.YAxis.Exponent =0;
print(h7,'-djpeg','-r200','Power')


% Plots for the PSI i CURVE
% Analytical
analyz = figure; 
ana = tiledlayout(1,2,"TileSpacing","tight");
nexttile
%PSI i CURVE No fringing
Psi_plot(2,Psi_fringing);
title('MEC with Air-gap Fringing');
xlabel('Current (A)');ylabel('Flux Linkage(Wb)')
nexttile
%PSI i CURVE No fringing
Psi_plot(2,Psi_nofringing);
title('MEC with Air-gap Fringing ');
xlabel('Current (A)');ylabel('Flux Linkage(Wb)')
title(ana,'Analytical Solution PSI-I Curves')
leg1 =legend('g_v=0.05 open','g_v=0.004','g_v=0.003','g_v=0.002','g_v=0.001','g_v=0.0001 closed');
title(leg1,'Air Gap Length (m)')
leg1.NumColumns = 3;
leg1.Layout.Tile = 'north';
print(analyz,'-djpeg','-r200','Analytical')

% FEMM
femmz = figure;
fem = tiledlayout(1,2,"TileSpacing","tight");
nexttile
%PSI i CURVE Non_Linear
Psi_plot(1,non_linear_psi);
title('FEMM:Non-Linear')
xlabel('Current (A)');ylabel('Flux Linkage(Wb)')
nexttile
%PSI i CURVE Linear 
Psi_plot(1,linear_psi);
title('FEMM:Linear')
title(fem,'FEA Numerical Solution PSI-I Curves')
xlabel('Current (A)');ylabel('Flux Linkage(Wb)')
leg =legend('g_v=0.05 open','g_v=0.004','g_v=0.003','g_v=0.002','g_v=0.001','g_v=0.0001 closed');
title(leg,'Air Gap Length (m)')
leg.NumColumns = 3;
leg.Layout.Tile = 'north';
print(femmz,'-djpeg','-r200','Femm Psi')


% Co energy Plots 
h2 = figure;
subplot(1,2,1)
plot(Armature_distance,Co_energy_linear,'-x','LineWidth',1.2)
title( 'FEMM:Linear')
xlabel('Armature Distance(m)');ylabel('Co-Energy(J)')
strValues = strtrim(cellstr(num2str([Co_energy_linear(:) ],'%.2d')));
text(Armature_distance,Co_energy_linear,strValues,'VerticalAlignment','bottom','FontSize',8);
ax = gca;
ax.XAxis.Exponent =0;
ax.YAxis.Exponent =0;
subplot(1,2,2)
plot(Armature_distance,Co_energy_nonlinear,'-x','LineWidth',1.2)
title('FEMM:Non Linear')
xlabel('Armature Distance(m)');ylabel('Co-Energy(J)')
ax = gca;
ax.XAxis.Exponent =0;
ax.YAxis.Exponent =0;
strValues2 = strtrim(cellstr(num2str([Co_energy_nonlinear(:) ],'%.2d')));
text(Armature_distance,Co_energy_nonlinear,strValues2,'HorizontalAlignment','right','FontSize',8)
print(h2,'-djpeg','-r300','Co_energy')

h4 = figure;
subplot(1,2,1)
plot(Armature_distance,Co_energy_fring,'-x','LineWidth',1.2)
title('MEC with Air-gap Fringing')
xlabel('Armature Distance(m)');ylabel('Co-Energy(J)')
ax = gca;
ax.XAxis.Exponent =0;
ax.YAxis.Exponent =0;
strValues2 = strtrim(cellstr(num2str([Co_energy_fring(:) ],'%.2d')));
text(Armature_distance,Co_energy_fring,strValues2,'HorizontalAlignment','right','FontSize',8)
subplot(1,2,2)
plot(Armature_distance,Co_energy_nofring,'-x','LineWidth',1.2)
title('MEC without Air-gap Fringing')
xlabel('Armature Distance(m)');ylabel('Co-Energy(J)')
ax = gca;
ax.XAxis.Exponent =0;
ax.YAxis.Exponent =0;
strValues2 = strtrim(cellstr(num2str([Co_energy_nofring(:)],'%.2d')));
text(Armature_distance,Co_energy_nofring,strValues2,'HorizontalAlignment','right','FontSize',8)
ylim([-0.01 0.45])
print(h4,'-djpeg','-r300','Co_energy_Ana')

% Comparison of FEA MODEL AND PSI-I CURVE METHOD
h5 = figure;
plot(Armature_distance,Co_energy_nonlinear,'k-x','LineWidth',1.2)
hold on 
plot(Armature_distance,Nonlinear_Femm_coenergy,'ro','MarkerSize',10,'LineWidth',1.2)
title('Change in Co-Energy Curve')
xlabel('Armature Distance(m)');ylabel('Co-Energy(J)')
legend({'PSI-I curve','FEA Model'},'Location','northwest')
print(h5,'-djpeg','-r300','Co_energy_compare')

%%
% Plot for Force Displacement 
h6=figure
title('Force Displacement Curves')
plot(mid_points,Force_Fring,'k-x','LineWidth',1.2)
hold on 
plot(mid_points,Force_No_Fring,'k-o','LineWidth',1.2)
plot(mid_points,Force_linear,'k--s','LineWidth',1.2)
plot(mid_points,Force_nonlinear,'k--d','LineWidth',1.2)
hold off
legend('MEC with Air-gap Fringing','MEC without Air-gap Fringing','FEMM:Linear','FEMM:Non Linear','Location','northwest')
ax = gca;
xlabel('Armature Distance(m)');ylabel('Force(N)')
ax.XAxis.Exponent =0;
ax.YAxis.Exponent =0;
xlim([0 0.005])
print(h6,'-djpeg','-r300','Forces')
%% Functions 
% FEA EQUATION
function select = output(prop_num)
% Changes properties and gets output values(Flux_Linkage;Resitance_1;Resitance_2;Inductance;Power )
    armature_motion = [0,-1,-2,-3,-4,-4.9];
    current = 0:1:10;
for gap = 1:length(armature_motion) 
   mi_selectgroup(2);
   % moving the armature
   mi_movetranslate(armature_motion(gap),0); 
   mi_clearselected;
   if prop_num == 1
           mi_selectlabel(10,9);
           mi_setblockprop('core_linear',0,0,'<None>',0,1,0)
           mi_clearselected()
           mi_selectlabel(70.1,0);
           mi_setblockprop('core_linear',0,0,'<None>',0,2,0)
           mi_clearselected()
   else
            mi_selectlabel(10,9);
            mi_setblockprop('core_nonlinear',0,0,'<None>',0,1,0)
            mi_clearselected()
            mi_selectlabel(70.1,0);
            mi_setblockprop('core_nonlinear',0,0,'<None>',0,2,0)
            mi_clearselected()
   end
   for I_value = 1:length(current)
    %-------- Setting Winidng Cuurent------------- 
    %winding 1
    mi_setcurrent('winding_1',current(I_value));
    %winding  2
    mi_setcurrent('winding_2',current(I_value));
    %---------MESHING THE MODEL--------------
     % Smart Mesh Control
%     mi_createmesh()
%     mi_showmesh()
    %-------Analyse the Model----------------- 
    mi_analyse()
    mi_loadsolution()
    % extract the circuit properties for the numerical method 
    CP = mo_getcircuitproperties('winding_1');
    CP2 = mo_getcircuitproperties('winding_2');

     % Measured Output from Femm  
    Flux_Linkage(I_value,gap) = CP(3);  
    Resitance_1(I_value) = division(CP(1),CP(2));
    Resitance_2(I_value)= division(CP2(1),CP2(2));
    Inductance(gap)= division(CP2(1),CP2(3));
    Power(I_value) =  CP(1)*CP(2);
    mo_groupselectblock()
    Co_energy(gap)= mo_blockintegral(17);
    mo_clearblock()
   end
    mi_selectgroup(2);
   % moving the armature
   mi_movetranslate((-1*armature_motion(gap)),0); 
   mi_clearselected;
end 
   select = {Flux_Linkage;Resitance_1;Resitance_2;Inductance;Power;Co_energy};
end 

% RESISTANCE
function Cw = division(Aw,Bw)
% Aw is the divisor  and Bw is the dividend
Cw = (Bw/Aw);
end
function Rw = resitance(n,l,c,k,a)
% n is the number of windings , l is the length of the windings,  c is the conductivity  = 58MS/m ,k is the packing factor 
Rw = ((n*l)/((c*(k*a))/n)); 
end

% Inductance 
function L_Inductance = core_fringing(prop_num)
%Calculates the Inductance for the Analytical Solution.
air_gap = [0.005, 0.004, 0.003, 0.002, 0.001,0.0001];
if prop_num == 1
    for y = 1: length(air_gap)
         x = air_gap(y);
            %L_Inductance(y) = 16/(((4000*x)/(pi*(x^2 + (x/25) + (1/2500)))) + (6831/(pi)));
           L_Inductance(y)= 10000/(((2500000 *x)/(pi*x^2 + ((pi* x)/25) + ((pi)/2500))) + (4268842.12018141/(pi)));
     end 
else 
   for y = 1: length(air_gap)
       x = air_gap(y);
       L_Inductance(y) = (16*pi)/(10000000*x + 7295);
   end 
end 
end

function Reluctance = core_fringing_reluc(prop_num)
%Calculates the Relucatance for the Analytical Simscape Solution.
air_gap = [0.005, 0.004, 0.003, 0.002, 0.001,0.0001];
if prop_num == 1
    for y = 1: length(air_gap)
         g_v = air_gap(y);
            
           Reluctance(y)= ((2500000 *g_v)/(pi*g_v^2 + ((pi*g_v)/25) + ((pi)/2500))) + (4268842.12018141/(pi));
     end 
else 
   for y = 1: length(air_gap)
       g_v = air_gap(y);
       Reluctance(y) = ((6250000000*g_v)/(pi)) + (4559375/(pi));
   end 
end 
end

% PSI - I 
function Psi_value = Psi(listname) 
   % Used to find PSI values 
   Current = 0:1:10;
   for i = 1: length(listname)
     for j = 1:length(Current)
       y = listname(i);
       Psi_value(i,j) = y.*Current(j);   
     end   
   end
end 

function plot_psi_curves  = Psi_plot(type,listname)
% Plots PSI-I CURVES
  Current = 0:1:10;
  if type == 1  
   for psi_number = 1:6
         plot(Current,listname(:,psi_number),'LineWidth',1.2)   
         hold on 
   end
   hold off
  else
    for psi_number = 1:6
         plot(Current,listname(psi_number,:),'LineWidth',1.2)   
         hold on 
    end
    hold off
  end

end 

% Co-energy
 function Coenergy = Integral(listname) 
% Used to find the change in Coenergy
Current = 0:1:10;
for inter = 1:min(size(listname))
    if inter == 1
    Co(inter)= trapz(Current,(listname(:,(inter))));    
    else 
    Co(inter)= trapz(Current,(listname(:,(inter)))) - trapz(Current,(listname(:,(inter-1))));
    end
    Coenergy = Co;
end 

 end
function Coenergy = Integralside(listname) 
% Used to find the change in Coenergy
Current = 0:1:10;
for inter = 1:min(size(listname))
    if inter == 1 
     Co(inter)=  trapz(Current,(listname(inter,:)));
    else 
     Co(inter)= trapz(Current,(listname((inter),:))) - trapz(Current,(listname((inter-1),:)));
    end
    Coenergy = Co;
end 

end

function Non_Linear_Coenergy = Coenergy(listname)
% Used to find the change in Coenergy
for inter = 1:6
    if inter == 1 
        Co(inter)= listname(inter);
    else 
        Co(inter)= listname(inter) - listname((inter-1));
    end
end
 Non_Linear_Coenergy = Co;
end 

% Force 
function Force = Gradient(listname)
% Used to find the Force
Armature_distance = [0,0.001,0.002,0.003,0.004,0.0049];
for i = 1: length(listname)
  if i > 1
        dW = listname(i); 
        dx = (Armature_distance(i) - Armature_distance(i-1));
        Force(i) = dW/dx;
  elseif i == 1
        dW = listname(i); 
        dx = Armature_distance(i);
        Force(i) = dW/dx;      
  end 
end
end 


  
 

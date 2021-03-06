%% -- MAIN -- %% 

%%%% Traitement de données de glider 
%%%% Calcul de vitesses verticales du courant à partir du modèle théorique
%%%% de Merchelback 
%%%% Auteur : Yann Darchen (yann.darchen@gmail.com) Juin 2020

%% ----- Initialisation ----- %%
% ------------------------------------- %
close all;clc;clear;clear global;
%Data path
addpath('Data_Stage');
addpath('Fonction');
addpath('SeaWater');


%data_name ='BioswotSeaExplorer_NavCTD_26Juin2020';
%data_name ='BioswotSeaExplorer_NavCTD_26Juin2020(2)';
%data_name ='BioswotSeaExplorer_Nav&CTD.mat';
%data_name ='Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat';
data_name='Fumseck_SeaExplorer_Nav&CTD_03Juillet2020.mat';
data=load(data_name);
%data.tableau=data.data; %% be careful to the variable of the load 
%Select Year
Year = [2019];
i_ign= data.tableau(:,2) ~= Year ;
ign=i_ign(:,1);
if length(Year)>1
   for i=1:length(Year)-1
       
       ign=logical(ign.*i_ign(:,i+1));
   end
end
data.tableau(ign,:)=[];  

%Delete False data  
to_igno = data.tableau(:,21) > 98;
data.tableau(to_igno,10)=NaN;%lat
data.tableau(to_igno,11)=NaN;%lon
data.tableau(to_igno,8)=NaN;%time
data.tableau(to_igno,9)=NaN;%depth
data.tableau(to_igno,14)=NaN;%oil
data.tableau(to_igno,18)=NaN;%pitch
data.tableau(to_igno,20)=NaN;%pressure
data.tableau(to_igno,21)=NaN;%temp
data.tableau(to_igno,1)=NaN;%dive
data.tableau(to_igno,22)=NaN;%cond

data.tableau(:,10) =fillmissing(data.tableau(:,10),'previous');
data.tableau(:,11)=fillmissing(data.tableau(:,11),'previous');
data.tableau(:,8)=fillmissing(data.tableau(:,8),'previous');
data.tableau(:,9)=fillmissing(data.tableau(:,9),'previous');
data.tableau(:,14)=fillmissing(data.tableau(:,14),'previous');
data.tableau(:,18)=fillmissing(data.tableau(:,18),'previous');
data.tableau(:,20)=fillmissing(data.tableau(:,20),'previous');
data.tableau(:,21)=fillmissing(data.tableau(:,21),'previous');
data.tableau(:,1)=fillmissing(data.tableau(:,1),'previous');
data.tableau(:,22)=fillmissing(data.tableau(:,22),'previous');


%%% Choose what dives you want
Quintuplet_desc_1 =[];
Quintuplet_desc_2 =[];
Quintuplet_mont_1 =[];
Quintuplet_mont_2 =[];
Precision_d1=[];
Precision_d2=[];
Precision_m1=[];
Precision_m2=[];
Fb_desc_1=[];
Ww_d1_final=[];
Ww_d2_final=[];
Ww_m1_final=[];
Ww_m2_final=[];
Triplet=[];
Triplet_diving=[];
Triplet_climbing=[];

for fn=0:12:24
  
explorer.first_dive =6+fn; %default is 1 
%explorer.last_dive = data.tableau(end,1); %default is last 
explorer.last_dive =17+fn;
explorer.first=explorer.first_dive; % in case you need to change this value 
explorer.last=explorer.last_dive;

%%% Choose what you want to display
yes=1;
no=0;

Display_Web_Map = no;
Display_Map_Figure = no;
Lat_Lon = no;
Depth_oil_pitch = no;
Colored_Temperature_Salinity = no;
Display_3D_trip = no; Animation_3D = no; %be careful animation can be very long 5 dives is recommended
Temperature_Salinity_Profiles = no;

%For vertical velocity 
vertical_velocities = yes; % always yes in order to calculate vertical velocities

%%% les optimisations 
opt_descente_1 = yes;
opt_descente_2 = yes;
opt_montee_1 = yes;
opt_montee_2 = yes;
opt_all = no;
opt_bydive = no;
opt_diving = no;
opt_climbing=no;

valeurs_opt=yes;
valeurs_moyennes=no;

Descente_1 =yes;
Descente_2 =yes;
Montee_1 =yes;
Montee_2 =yes;
Together = yes;
All = no;
Diving=no;
Climbing=no; 

W_glider_W_Model = no;  attack_angle_bydive = no; attack_angle_all = no;
Parameters_evolution = no; Parameters_evolution_bydive =no;
water_velocity_descent_bydive = no; hist_desc = no;
water_velocity_ascent_by_dive=no; 
water_velocity_descent_all = no;
water_velocity_ascent_all = no; hist_asc=no; sbplot = no;
sub_asc_desc = no;
water_velocity_all_bydive = no;
water_velocity_diving = no;


%Constant
Const.d2s = 86400;
Const.g=9.81; % gravité
Const.R = 6371;
Const.V0 = 0.05944;
Const.eps_d= 3.9230e-10;
Const.eps_m=4.0410e-10;
Const.Cd_d = 0.1039;
Const.Cd_m = 0.995;
Const.aw=2.6;
Const.Cdw=2.6;

%load(nom_de_fichier);
load('blue_red_cmap.mat')
%% ----- Read Data ----- %%
% ------------------------------------- %

%dives wanted
i_dive= data.tableau(:,1) >= explorer.first_dive & data.tableau(:,1) <= explorer.last_dive ;
tableau = data.tableau(i_dive,:);   

%Read Data
explorer_all = read_EXPLORER(tableau,explorer,Const);
explorer_bydive=by_dive2(tableau,explorer,Const);

%% ----- Display Map Web Browser ----- %%
% ------------------------------------- %

if Display_Web_Map == 1
    webmap('Ocean Basemap')
    indices = linspace(1,explorer_all.size,200);% nombre de points à afficher
    indices=fix(indices);

    %wmmarker(explorer.lat(indices),explorer.lon(indices));
    wmmarker(explorer_all.lat(1),explorer_all.lon(1),'Color','green');
    wmmarker(explorer_all.lat(end),explorer_all.lon(end),'Color','red');

    wmline(explorer_all.lat(indices),explorer_all.lon(indices));
end

%% ----- Display Map Figure + gif -------------------------- %%
% ----------------------------------------------------------- %
if Display_Map_Figure == 1
    
    pause(1.5)
    %h = figure;
    %axis off
    %filename = 'trajectoire.gif';
    figure('Name','Map','NumberTitle','off','Units','normalized','Position',[0,0.5,0.33,0.42]);
    load coastlines
    axesm('ortho','origin',[45 0]);
    Mean_lat = mean(explorer_all.lat);
    Mean_lon = mean(explorer_all.lon);
    axesm('mercator','MapLatLimit',[Mean_lat-2 Mean_lat+2],'MapLonLimit',[Mean_lon-2 Mean_lon+2])
    axis off;
    gridm off;
    framem on;
    %mlabel('equator')
    %plabel('fontweight','bold')

    hold on 
    plotm(coastlat,coastlon)
    plotm(explorer_all.lat(1),explorer_all.lon(1),'+g')
     for i=2:200:explorer_all.size-100
         plotm(explorer_all.lat(i),explorer_all.lon(i),'.r')
         %title(['Trajectoire du glider jusquau ' num2str(position(i,3)) '/' num2str(position(i,4)) '/' num2str(position(i,5)) '  ' num2str(position(i,6)) ':' num2str(position(i,7))])
         drawnow
         % Capture the plot as an image 
          %frame = getframe(h); 
          %im = frame2im(frame); 
          %[imind,cm] = rgb2ind(im,256); 
          % Write to the GIF File 
          %if i == 2 
            %  imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05); 
          %else 
           %   imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05); 
          %end 
    end

    plotm(explorer_all.lat(end),explorer_all.lon(end),'+b')
    hold off

end
%% ----- Display latitude and longitude ----- %%
% ----------------------------------------------- %


if Lat_Lon == 1
    
    pause(1.5)
    figure('Name','Latitude and Longitude','NumberTitle','off','Units','normalized','Position',[0.33,0.5,0.33,0.42]);
    subplot(2,1,1)
    plot(explorer_all.time,explorer_all.lat,'LineWidth',2)
    %plot(explorer_all.lat,'LineWidth',2)
    title('Evolution de la latitude')
    ylabel('latitude en °')
    datetick('x',2,'keepticks')

    subplot(2,1,2)
    plot(explorer_all.time,explorer_all.lon,'LineWidth',2)
    title('Evolution de la longitude')
    ylabel('longitude en °')
    datetick('x',2,'keepticks')

end

%% %% ----- Display depth, oil volume, pitch ----- %%
% ------------------------------------- %

if Depth_oil_pitch == 1
    

    pause(1.5)
    figure('Name','depth oil volume and pitch','NumberTitle','off','Units','normalized','Position',[0.66,0.5,0.33,0.42]);
    yyaxis left
    plot(explorer_all.time,-explorer_all.depth,'k','LineWidth',1)
    hold on 
    plot(explorer_all.time,explorer_all.oil,'-+','LineWidth',1)
    ylabel('profondeur(m), huile (ml)')
    yyaxis right
    %plot(explorer_all.time,explorer_all.pitch_filter,'-+','LineWidth',1)
    plot(explorer_all.time,explorer_all.pitch,'-+','LineWidth',1)
    ylabel('pitch (rad)')
    legend('profondeur','huile','pitch')
    title('Profondeur, volume huile et pitch')
    datetick('x',0,'keepticks')
    
    figure('Name','depth oil volume and pitch','NumberTitle','off','Units','normalized','Position',[0.66,0.5,0.33,0.42]);
    yyaxis left
    plot(explorer_all.dive,-explorer_all.depth,'LineWidth',1)
    hold on 
    plot(explorer_all.dive,explorer_all.oil,'LineWidth',1)
    ylabel('profondeur(m), huile (ml)')
    %yyaxis right
    %plot(explorer_all.dive,explorer_all.pitch_filter*(pi/180),'-+','LineWidth',1)
    
    hold off
    title('Profondeur, huile, pitch')
    ylabel('pitch (°)')
    %datetick('x',0,'keepticks')
    legend('profondeur','volume huile','pitch')

    figure()
    yyaxis left
    plot(tableau(:,8),tableau(:,19))
    ylabel('Roll (°)')
    hold on 
    yyaxis right
    plot(explorer_all.time,explorer_all.pitch_filter,'-+','LineWidth',1)
    title('Roll and pitch')
    ylabel('Pitch (°)')
    datetick('x',0,'keepticks')
    legend('Roll','pitch')
    
end

%% %% ----- Display temperature and salinity ----- %%
% ------------------------------------- %

% %partition in dive

if Colored_Temperature_Salinity == 1 % only descents 
    pause(1.5)
    
    tab_ti=[];
    tab_si=[];

    for k=1:explorer.last_dive-explorer.first_dive+1 
      
        explorer2_bydive.pressure=explorer_bydive.pressure(:,k);% struct editable
        explorer2_bydive.temp=explorer_bydive.temp(:,k);
        explorer2_bydive.s=explorer_bydive.s(:,k);
        
%         to_ign=[];
%         max_ind = find(explorer2_bydive.pressure == max(explorer2_bydive.pressure));
%         i=1;
%         while explorer2_bydive.pressure(i) ~= max(explorer_bydive.pressure)
%           if  not(explorer2_bydive.pressure(i)<explorer2_bydive.pressure(i+1)) % delete variations
%             to_ign = [to_ign i+1];
%           end
%           i=i+1;
%         end
%         to_ign = [to_ign max_ind:length(explorer2_bydive.pressure)]; % taking only descent 
          explorer2_bydive.p_sorted = explorer2_bydive.pressure;
          C=NaN(length(explorer2_bydive.pressure),3);
          C(:,1)=explorer2_bydive.p_sorted;
          C(:,2)=explorer2_bydive.temp;
          C(:,3)=explorer2_bydive.s;
          D=sortrows(C,1);
%         explorer2_bydive.p_sorted(to_ign)=[];
%         explorer2_bydive.temp(to_ign)=[];
%         explorer2_bydive.s(to_ign)=[];
        to_igno = [];
        for i=1:length(D)-1
         if  D(i,1) == D(i+1,1) % delete doublon
           to_igno = [to_igno i+1];
         end
        end

        D(to_igno,:)=[];
        ign = isnan(D(:,1));
        D(ign,:)=[];
        D=real(D);
%         explo(to_igno)=[];
%         explorer2_bydive.s(to_igno)=[];
%         explorer2_bydive.temp=fillmissing(explorer2_bydive.temp,'previous');
%         explorer2_bydive.s=fillmissing(explorer2_bydive.s,'previous');
        pi=[0:1:max(explorer_all.pressure)];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
        ti=interp1(D(:,1),D(:,2),pi); 
        si=interp1(D(:,1),D(:,3),pi);
        
        ind = abs(ti) > mean(explorer_all.temp)+4 | abs(ti) < mean(explorer_all.temp)-4;
        ti(ind)=NaN;
        ind = abs(si) > real(mean(explorer_all.s))+4 | abs(si) < real(mean(explorer_all.s))-4;
        si(ind)=NaN;
        
    
        tab_ti=[tab_ti ti'];
        tab_si=[tab_si si'];

    end

    figure('Name','Temperature and Salinity','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38]);
    subplot(2,1,1)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tab_ti);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
   % H.Limits=[mean(explorer_all.temp)-4,mean(explorer_all.temp)];
    ylabel(H,'       (°C)','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Dive number')
    ylabel('- Pressure (dbar)')
    title('Temperature')

    subplot(2,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tab_si);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'      psu','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('dive number')
    ylabel('- Pressure (dbar)')
    title('Salinity')

end


  
%% %% ----- Display 3D ----- %%
% ------------------------------------- %

if Display_3D_trip == 1
    pause(0.1)
    ind = [];
    for j = 1:explorer_all.size-1
        if explorer_all.dive(j) ~= explorer_all.dive(j+1) && mod(explorer_all.dive(j),2)==0 %collect first points of dives 
           ind  = [ind j+1]; 
        end
    end

    explorer.lat_i=explorer_all.lat(ind);
    explorer.time_i=explorer_all.time(ind);
    time_interp = linspace(explorer.time_i(1),explorer.time_i(end),(explorer.last_dive-explorer.first_dive)*150);
    %figure()
    explorer.lat_interp = interp1(explorer.time_i,explorer.lat_i,time_interp,'spline');
    %plot(explorer.time_i,explorer.lat_i,'o',time_interp,explorer.lat_interp,':.');

    explorer.lon_i=explorer_all.lon(ind);
    %figure()
    explorer.lon_interp = interp1(explorer.time_i,explorer.lon_i,time_interp,'spline');
    %plot(explorer.time_i,explorer.lon_i,'o',time_interp,explorer.lon_interp,':.');
    
    %Filter on explorer.time
    to_ign=[];
    for i=1:length(explorer_all.time)-1
        if explorer_all.time(i)> explorer_all.time(i+1) || explorer_all.time(i)==explorer_all.time(i+1)
            to_ign=[to_ign i];
        end
    end
    explorer2=explorer_all;
    
    explorer2.time(to_ign)=[];
    explorer2.depth(to_ign)=[];
    explorer2.temp(to_ign)=[];

    explorer.depth_interp = interp1(explorer2.time,explorer2.depth,time_interp,'spline');
    explorer.temp_interp = interp1(explorer2.time,explorer2.temp,time_interp,'pchip');
    %plot(explorer.time,explorer.temp,'o',time_interp,explorer.temp_interp,':.')
    
    X = explorer.lat_interp;
    Y = explorer.lon_interp;
    Z = -explorer.depth_interp;
    C = explorer.temp_interp;
    Z(end)=NaN;
    C(end)=NaN;
%     h = figure;
%     %axis off
%     filename = 'trajectoire.gif';
    figure('Name','3D_Trip_Temperature','NumberTitle','off','Units','normalized','Position',[0.33,0.05,0.33,0.38]);

    fill3(X,Y,Z,C,'EdgeColor','w','LineWidth',1)
    c = colorbar;
    c.Label.String = 'Température °C';
    title('Trajectoire du glider')
    xlabel('latitude')
    ylabel('longitude')
    zlabel('profondeur (m)')
    hold on 
    for j= 1:length(explorer.lat_interp)

        fill3(X(j),Y(j),Z(j),C(j),'EdgeColor','interp','Marker','o','MarkerFaceColor','flat')
        if Animation_3D == 1 
            drawnow
%             % Capture the plot as an image 
%           frame = getframe(h); 
%           im = frame2im(frame); 
%           [imind,cm] = rgb2ind(im,256); 
%           % Write to the GIF File 
%           if i == 2 
%               imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05); 
%           else 
%               imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05); 
%           end 
        end
    end
    hold off
    v = [1 -3 3];
    [caz,cel] = view(v);
end

%% %% ----- Display temperature and salinity profile ----- %%
% ------------------------------------- %

if Temperature_Salinity_Profiles == 1
    pause(1.5)
    
    
    
    figure('Name','Temperature_Salinity_Density_Profiles','NumberTitle','off','Units','normalized','Position',[0.66,0.05,0.33,0.38]);
    subplot(1,3,1)
    plot(explorer_all.temp,-explorer_all.pressure,'b')
    title('Température')
    ylabel('Pression (dbar)')
    xlabel('°C')
    subplot(1,3,2)
    plot(explorer_all.s,-explorer_all.pressure,'r')
    title('Salinité')
    xlabel('psu')
    subplot(1,3,3)
    plot(explorer_all.dens,-explorer_all.pressure,'k')
    title('Masse volumique')
    xlabel('kg/m^3')

end

%% %% ----- Display vertical velocities and attack angle ----- %%
% ------------------------------------- %

% W : vertical water velocity
% W_model : velocity from the model
% W_glider : Pressure/time
% W = W_glider - W_model
if vertical_velocities == 1
    pause(1.5)
  
%%%---------------------------------------------------------------------%
%%%%%%%%%%%%%%%% Méthode fminsearch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---------------------------------------------------------------------%

 %%%%% Optimisation sur couple descente-montée %%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%---------------------------------------------------------------------%
 if opt_bydive == 1 
     
     W_model_all=[];
     W_glider_all=[];
     time_all=[];
     Tri=[];
     tab_att=[];
     
      explorer3_bydive=explorer_bydive;
       explorer3_bydive.W_glider_acc = explorer3_bydive.W_glider_acc*explorer_all.M;
       ign_acc = explorer3_bydive.W_glider_acc > 0.007 | explorer3_bydive.W_glider_acc < -0.007;
       explorer3_bydive.pressure(ign_acc)=NaN;
       explorer3_bydive.dens(ign_acc)=NaN;
       explorer3_bydive.pitch_filter(ign_acc)=NaN;
       explorer3_bydive.oil(ign_acc)=NaN;
       explorer3_bydive.W_glider_filter(ign_acc)=NaN;
       explorer3_bydive.temp(ign_acc)=NaN;
       explorer3_bydive.time(ign_acc)=NaN;
       explorer3_bydive.s(ign_acc)=NaN;
      
    for k=1:explorer.last_dive-explorer.first_dive+1
    
       
        
        global pres dens pitch oil_vol Wglider temp mg
        pres=explorer3_bydive.pressure(:,k);
        dens=explorer3_bydive.dens(:,k);
        pitch=explorer3_bydive.pitch_filter(:,k);
        oil_vol=explorer3_bydive.oil(:,k);
        Wglider=explorer3_bydive.W_glider_filter(:,k)';
        temp=explorer3_bydive.temp(:,k);
        mg=explorer_all.M;
       aw=3.7;%initial
       Cd1w = 0.78;%initial
       
       param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
       options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
       [x_bydive] = fminsearch('cost',param0,options);
          
       global V0 eps Cd
        V0=x_bydive(1);
        eps=x_bydive(2);
        Cd=x_bydive(3);
        param0=[aw,Cd1w];
        [x_bydive2] = fminsearch('cost_wings',param0,options); 
       
        x_bydive1=[x_bydive x_bydive2];
          
        %[W_model3,att_deg] = flight_model(pres,dens,pitch,oil_vol,temp,x_bydive(1),x_bydive(2),x_bydive(3),mg);
         [W_model3,att_deg] = flight_model2(pres,dens,pitch,oil_vol,temp,x_bydive(1),x_bydive(2),x_bydive(3),mg,x_bydive2(1),x_bydive2(2));
         %W_model = fillmissing(W_model,'next');
         
         %[W_model] = flight_model(pres,dens,pitch,oil_vol,temp,param0(1),param0(2),param0(3),mg);
         Tri = [Tri x_bydive1'];
         W_model_all= [W_model_all W_model3'];
         W_glider_all = [W_glider_all explorer3_bydive.W_glider_filter(:,k)'];
         time_all = [time_all explorer_bydive.time(:,k)'];
         tab_att= [tab_att att_deg];
    
    end
    %[W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);
   
 
 end
if attack_angle_bydive == 1 
    pause(1.5)
    figure('Name','Attack angle','NumberTitle','off','Units','normalized','Position',[0,0.5,0.33,0.42]);
    histogram(tab_att,200)
    title('Attack angle bydive')
    xlabel('Attack angle °')
    xlim([-20 20])

end

%   W_glider_lisse = smoothdata(explorer_all.W_glider_filter,'SmoothingFactor',0.02);
%   Pitch_lisse = smoothdata(explorer_all.pitch_filter','SmoothingFactor',0.017);
%   Temp_lisse = smoothdata(explorer_all.temp,'SmoothingFactor',0.02);

 %%%%% Optimisation sur tout %%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%---------------------------------------------------------------------%
 
 if opt_all == 1
 
tab_all.tab=[];
 
tab_all.tab(:,1)=explorer_all.pressure; %Création du tableau
tab_all.tab(:,2)=explorer_all.dens;
tab_all.tab(:,3)=explorer_all.pitch_filter';
%tab_all.tab(:,3)=explorer_all.pitch; %non filtré
tab_all.tab(:,4)=explorer_all.oil;
tab_all.tab(:,5)=explorer_all.W_glider_filter';
%tab_all.tab(:,5)=explorer_all.W_glider;%non filtré
tab_all.tab(:,6)=explorer_all.temp;
tab_all.tab(:,7)=explorer_all.time;
tab_all.tab(:,8)=explorer_all.W_glider_acc*explorer_all.M;

ign1 = abs(tab_all.tab(:,8)) > 0.01; %filtre sur l'accélération 
tab_all.tab(ign1,:)=[];

ign2=tab_all.tab(:,1) > 570 | tab_all.tab(:,1) < 40;
tab_all.tab(ign2,:)=[];


        global pres dens pitch oil_vol Wglider temp mg
        pres=tab_all.tab(:,1);
        dens=tab_all.tab(:,2);
        pitch=tab_all.tab(:,3);
        oil_vol=tab_all.tab(:,4);
        Wglider=tab_all.tab(:,5)';
        temp=tab_all.tab(:,6);
        mg=explorer_all.M;
   
        param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
        options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
        [x_all] = fminsearch('cost',param0,options);
        Triplet=[Triplet x_all'];
         global V0 eps Cd
         V0=x_all(1);
         eps=x_all(2);
         Cd=x_all(3);
         aw=3.7;
         Cd1w=0.78;
         param0=[aw,Cd1w];
         [x_all2] = fminsearch('cost_wings',param0,options);
      
         [W_model_all,att_deg] = flight_model2(pres,dens,pitch,oil_vol,temp,x_all(1),x_all(2),x_all(3),mg,x_all2(1),x_all2(2));
        Ww_all = tab_all.tab(:,5)-W_model_all;
       % Ww_model_one_quintuplet= flight_model2(pres,dens,pitch,oil_vol,temp,0.59945,4.0265e-10,0.1040,mg,2.6,2.6);
        %Ww_all_one_quintuplet = tab_all.tab(:,5)-W_model_one_quintuplet;
 end
 %%%%% Optimisation sur descente  %%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%---------------------------------------------------------------------%
 if opt_diving == 1
     
tab_diving.tab=[];

tab_diving.tab(:,1)=explorer_all.pressure;
tab_diving.tab(:,2)=explorer_all.dens;
tab_diving.tab(:,3)=explorer_all.pitch_filter';
tab_diving.tab(:,4)=explorer_all.oil;
tab_diving.tab(:,5)=explorer_all.W_glider_filter';
tab_diving.tab(:,6)=explorer_all.temp;
tab_diving.tab(:,7)=explorer_all.time;
tab_diving.tab(:,8)=explorer_all.W_glider_acc*explorer_all.M;

ign1 = abs(tab_diving.tab(:,8)) > 0.01; %filtre sur l'accélération 
tab_diving.tab(ign1,:)=[];

ign = tab_diving.tab(:,3) > 0 | tab_diving.tab(:,4) > 0; %enleve les montées
tab_diving.tab(ign,:)=[];

ign2=tab_diving.tab(:,1) > 500 | tab_diving.tab(:,1) < 100;
tab_diving.tab(ign2,:)=[];


% Mean_glider= mean(tab(:,5));
% ign1=tab(:,5) > Mean_glider +0.05;
% tab(ign1,:)=[];
%  ign2=tab(:,1) > 500 | tab(:,1) < 100;
%  tab(ign2,:)=[];

     global pres dens pitch oil_vol Wglider temp mg
        pres=tab_diving.tab(:,1);
        dens=tab_diving.tab(:,2);
        pitch=tab_diving.tab(:,3);
        oil_vol=tab_diving.tab(:,4);
        Wglider=tab_diving.tab(:,5)';
        temp=tab_diving.tab(:,6);
        mg=explorer_all.M;

        param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
        options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
       
        [x_diving] = fminsearch('cost',param0,options);
        Triplet_diving=[Triplet_diving x_diving'];
     
% [W_model_diving] = flight_model2(pres,dens,pitch,oil_vol,temp,test(1),test(2),test(3),mg,test(4),test(5));
 [W_model_diving,att_deg_diving,U,att,Fg,Fb] = flight_model(pres,dens,pitch,oil_vol,temp,x_diving(1),x_diving(2),x_diving(3),mg);
 Ww_diving =tab_diving.tab(:,5)-W_model_diving;
 end
 
if attack_angle_all == 1 
    pause(1.5)
    figure('Name','Attack angle','NumberTitle','off','Units','normalized','Position',[0.33,0.5,0.33,0.42]);
    histogram(att_deg,200)
    title('Attack angle all')
    xlabel('Attack angle °')
    xlim([-20 20])

end

%%%%% Optimisation sur montée  %%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%---------------------------------------------------------------------%
 if opt_climbing == 1
     
tab_climbing.tab=[];

tab_climbing.tab(:,1)=explorer_all.pressure;
tab_climbing.tab(:,2)=explorer_all.dens;
tab_climbing.tab(:,3)=explorer_all.pitch_filter';
tab_climbing.tab(:,4)=explorer_all.oil;
tab_climbing.tab(:,5)=explorer_all.W_glider_filter';
tab_climbing.tab(:,6)=explorer_all.temp;
tab_climbing.tab(:,7)=explorer_all.time;
tab_climbing.tab(:,8)=explorer_all.W_glider_acc*explorer_all.M;

ign1 = abs(tab_climbing.tab(:,8)) > 0.01; %filtre sur l'accélération 
tab_climbing.tab(ign1,:)=[];

ign = tab_climbing.tab(:,3) < 0 | tab_climbing.tab(:,4) < 0; %enleve les descentes
tab_climbing.tab(ign,:)=[];

ign2=tab_climbing.tab(:,1) > 500 | tab_climbing.tab(:,1) < 100;
tab_climbing.tab(ign2,:)=[];

     global pres dens pitch oil_vol Wglider temp mg
        pres=tab_climbing.tab(:,1);
        dens=tab_climbing.tab(:,2);
        pitch=tab_climbing.tab(:,3);
        oil_vol=tab_climbing.tab(:,4);
        Wglider=tab_climbing.tab(:,5)';
        temp=tab_climbing.tab(:,6);
        mg=explorer_all.M;

        param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
        options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
       
        [x_climbing] = fminsearch('cost',param0,options);
        Triplet_climbing=[Triplet_climbing x_climbing'];
     
% [W_model_diving] = flight_model2(pres,dens,pitch,oil_vol,temp,test(1),test(2),test(3),mg,test(4),test(5));
 [W_model_climbing,att_deg_diving,U,att,Fg,Fb] = flight_model(pres,dens,pitch,oil_vol,temp,x_climbing(1),x_climbing(2),x_climbing(3),mg);
 Ww_climbing =tab_climbing.tab(:,5)-W_model_climbing;
 end
 
if attack_angle_all == 1 
    pause(1.5)
    figure('Name','Attack angle','NumberTitle','off','Units','normalized','Position',[0.33,0.5,0.33,0.42]);
    histogram(att_deg,200)
    title('Attack angle all')
    xlabel('Attack angle °')
    xlim([-20 20])

end
 %%%%% Optimisation sur 1ère descente  %%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%---------------------------------------------------------------------%

 if opt_descente_1 == 1 
   
 %%%Préparation des données %%%%%% 
 
   tab_descent_1.tab=[];
   tab_descent_1.pressure=[];  
   tab_descent_1.dens=[];
   tab_descent_1.pitch=[];
   tab_descent_1.oil=[];
   tab_descent_1.w_glider=[];
   tab_descent_1.temp=[];
   tab_descent_1.time=[];
   tab_descent_1.acc=[];
   
     for k=1:2:explorer.last_dive-explorer.first_dive+1
         
       
        tab_descent_1.pressure=[tab_descent_1.pressure explorer_bydive.pressure(:,k)']; % On les remplit avec montée et descente une plongée sur deux 
        tab_descent_1.dens=[tab_descent_1.dens explorer_bydive.dens(:,k)'];
        tab_descent_1.pitch=[tab_descent_1.pitch explorer_bydive.pitch_filter(:,k)'];
        tab_descent_1.oil=[tab_descent_1.oil explorer_bydive.oil(:,k)'];
        tab_descent_1.w_glider=[tab_descent_1.w_glider explorer_bydive.W_glider_filter(:,k)'];
        tab_descent_1.temp=[tab_descent_1.temp explorer_bydive.temp(:,k)'];
        tab_descent_1.time=[tab_descent_1.time explorer_bydive.time(:,k)'];
        tab_descent_1.acc=[tab_descent_1.acc (explorer_bydive.W_glider_acc(:,k).*explorer_all.M)'];


     end
     
tab_descent_1.tab(:,1)=tab_descent_1.pressure; %Création du tableau
tab_descent_1.tab(:,2)=tab_descent_1.dens;
tab_descent_1.tab(:,3)=tab_descent_1.pitch;
tab_descent_1.tab(:,4)=tab_descent_1.oil;
tab_descent_1.tab(:,5)=tab_descent_1.w_glider;
tab_descent_1.tab(:,6)=tab_descent_1.temp;
tab_descent_1.tab(:,7)=tab_descent_1.time;
tab_descent_1.tab(:,8)=tab_descent_1.acc;

ign = tab_descent_1.tab(:,4) > 0 | tab_descent_1.tab(:,3) > 0 ; %ignorer les montées donc pitch > 0
tab_descent_1.tab(ign,:)=[];
% 
ign1 = abs(tab_descent_1.tab(:,8)) > 0.005;
tab_descent_1.tab(ign1,:)=[];

ign2=tab_descent_1.tab(:,1) > 500 | tab_descent_1.tab(:,1) < 100;
tab_descent_1.tab(ign2,:)=[];
 
%%%Optimisation %%%%%% 

global pres dens pitch oil_vol Wglider temp mg 
pres=tab_descent_1.tab(:,1);
dens=tab_descent_1.tab(:,2);
pitch=tab_descent_1.tab(:,3);
oil_vol=tab_descent_1.tab(:,4);
Wglider=tab_descent_1.tab(:,5)';
temp=tab_descent_1.tab(:,6);
mg=explorer_all.M;


%%%Première optimisation
aw=3.7;
Cd1w = 0.78;
param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
[x_desc_1,precision_descent_1] = fminsearch('cost',param0,options); 
precision_descent_1=sqrt(precision_descent_1/length(pres));
%%%Deuxième optimisation
global V0 eps Cd
V0=x_desc_1(1);
eps=x_desc_1(2);
Cd=x_desc_1(3);
param0=[aw,Cd1w];
[x_desc_11,precision_descent_11] = fminsearch('cost_wings',param0,options);

Precision_d1=[Precision_d1 precision_descent_11];
Quintuplet_desc_1 =[Quintuplet_desc_1 [x_desc_1 x_desc_11]'];
[W_model_desc_1,att_deg,U,att,Fg,Fb_d1,Fl,Fd,vx_d1] = flight_model2(pres,dens,pitch,oil_vol,temp,x_desc_1(1),x_desc_1(2),x_desc_1(3),mg,x_desc_11(1),x_desc_11(2));
ind=Fb_d1 ==0;
Fb_d1(ind)=[];
Fbd1=nanmean(Fb_d1);
Fb_desc_1=[Fb_desc_1 Fbd1];

Ww_desc_1 = tab_descent_1.tab(:,5)-W_model_desc_1;

[W_model_desc_1_moy] = flight_model2(pres,dens,pitch,oil_vol,temp,Const.V0,Const.eps_d,Const.Cd_d,mg,Const.aw,Const.Cdw);
Ww_desc_1_moy = tab_descent_1.tab(:,5)-W_model_desc_1_moy;
 end
 
 %%%%% Optimisation sur 2ème descente  %%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%---------------------------------------------------------------------%
 
 if opt_descente_2 == 1 
    
     %%%Préparation des données %%%%%% 
 
   tab_descent_2.tab=[];
   tab_descent_2.pressure=[];  
   tab_descent_2.dens=[];
   tab_descent_2.pitch=[];
   tab_descent_2.oil=[];
   tab_descent_2.w_glider=[];
   tab_descent_2.temp=[];
   tab_descent_2.time=[];
   tab_descent_2.acc=[];
 
     for k=2:2:explorer.last_dive-explorer.first_dive+1
         
       
        tab_descent_2.pressure=[tab_descent_2.pressure explorer_bydive.pressure(:,k)']; % On les remplit avec montée et descente une plongée sur deux 
        tab_descent_2.dens=[tab_descent_2.dens explorer_bydive.dens(:,k)'];
        tab_descent_2.pitch=[tab_descent_2.pitch explorer_bydive.pitch_filter(:,k)'];
        tab_descent_2.oil=[tab_descent_2.oil explorer_bydive.oil(:,k)'];
        tab_descent_2.w_glider=[tab_descent_2.w_glider explorer_bydive.W_glider_filter(:,k)'];
        tab_descent_2.temp=[tab_descent_2.temp explorer_bydive.temp(:,k)'];
        tab_descent_2.time=[tab_descent_2.time explorer_bydive.time(:,k)'];
        tab_descent_2.acc=[tab_descent_2.acc (explorer_bydive.W_glider_acc(:,k).*explorer_all.M)'];


     end
     
tab_descent_2.tab(:,1)=tab_descent_2.pressure; %Création du tableau
tab_descent_2.tab(:,2)=tab_descent_2.dens;
tab_descent_2.tab(:,3)=tab_descent_2.pitch;
tab_descent_2.tab(:,4)=tab_descent_2.oil;
tab_descent_2.tab(:,5)=tab_descent_2.w_glider;
tab_descent_2.tab(:,6)=tab_descent_2.temp;
tab_descent_2.tab(:,7)=tab_descent_2.time;
tab_descent_2.tab(:,8)=tab_descent_2.acc;

ign = tab_descent_2.tab(:,4) > 0 | tab_descent_2.tab(:,3) > 0 ; %ignorer les montées donc pitch > 0
tab_descent_2.tab(ign,:)=[];

ign1 = abs(tab_descent_2.tab(:,8)) > 0.005; %ignorer les accélérations trop fortes 
tab_descent_2.tab(ign1,:)=[];

ign2=tab_descent_2.tab(:,1) > 500 | tab_descent_2.tab(:,1) < 100;
tab_descent_2.tab(ign2,:)=[];


%%%Optimisation %%%%%% 

global pres dens pitch oil_vol Wglider temp mg
pres=tab_descent_2.tab(:,1);
dens=tab_descent_2.tab(:,2);
pitch=tab_descent_2.tab(:,3);
oil_vol=tab_descent_2.tab(:,4);
Wglider=tab_descent_2.tab(:,5)';
temp=tab_descent_2.tab(:,6);
mg=explorer_all.M;

%%%Première optimisation
aw=3.7;
Cd1w = 0.78;
param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
[x_desc_2,precision_descent_2] = fminsearch('cost',param0,options); 
precision_descent_2=sqrt(precision_descent_2/length(pres));
%%%Deuxième optimisation
global V0 eps Cd
V0=x_desc_2(1);
eps=x_desc_2(2);
Cd=x_desc_2(3);
param0=[aw,Cd1w];
[x_desc_22,precision_descent_22] = fminsearch('cost_wings',param0,options);

Precision_d2=[Precision_d2 precision_descent_22];
Quintuplet_desc_2 =[Quintuplet_desc_2 [x_desc_2 x_desc_22]'];
[W_model_desc_2,att_deg,U,att,Fg,Fb_d1,Fl,Fd,vx_d2] = flight_model2(pres,dens,pitch,oil_vol,temp,x_desc_2(1),x_desc_2(2),x_desc_2(3),mg,x_desc_22(1),x_desc_22(2));
Ww_desc_2 = tab_descent_2.tab(:,5)-W_model_desc_2;

[W_model_desc_2_moy] = flight_model2(pres,dens,pitch,oil_vol,temp,Const.V0,Const.eps_d,Const.Cd_d,mg,Const.aw,Const.Cdw);
Ww_desc_2_moy = tab_descent_2.tab(:,5)-W_model_desc_2_moy;
 end
 
 


  %%%%% Optimisation sur 1ère montée  %%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%---------------------------------------------------------------------%
 
 if opt_montee_1 == 1 
    
     %%%Préparation des données %%%%%% 
 
   tab_montee_1.tab=[];
   tab_montee_1.pressure=[];  
   tab_montee_1.dens=[];
   tab_montee_1.pitch=[];
   tab_montee_1.oil=[];
   tab_montee_1.w_glider=[];
   tab_montee_1.temp=[];
   tab_montee_1.time=[];
   tab_montee_1.acc=[];
 
     for k=1:2:explorer.last_dive-explorer.first_dive+1
         
       
        tab_montee_1.pressure=[tab_montee_1.pressure explorer_bydive.pressure(:,k)']; % On les remplit avec montée et descente une plongée sur deux 
        tab_montee_1.dens=[tab_montee_1.dens explorer_bydive.dens(:,k)'];
        tab_montee_1.pitch=[tab_montee_1.pitch explorer_bydive.pitch_filter(:,k)'];
        tab_montee_1.oil=[tab_montee_1.oil explorer_bydive.oil(:,k)'];
        tab_montee_1.w_glider=[tab_montee_1.w_glider explorer_bydive.W_glider_filter(:,k)'];
        tab_montee_1.temp=[tab_montee_1.temp explorer_bydive.temp(:,k)'];
        tab_montee_1.time=[tab_montee_1.time explorer_bydive.time(:,k)'];
        tab_montee_1.acc=[tab_montee_1.acc (explorer_bydive.W_glider_acc(:,k).*explorer_all.M)'];


     end
     
tab_montee_1.tab(:,1)=tab_montee_1.pressure; %Création du tableau
tab_montee_1.tab(:,2)=tab_montee_1.dens;
tab_montee_1.tab(:,3)=tab_montee_1.pitch;
tab_montee_1.tab(:,4)=tab_montee_1.oil;
tab_montee_1.tab(:,5)=tab_montee_1.w_glider;
tab_montee_1.tab(:,6)=tab_montee_1.temp;
tab_montee_1.tab(:,7)=tab_montee_1.time;
tab_montee_1.tab(:,8)=tab_montee_1.acc;

 ign = tab_montee_1.tab(:,4) < 0 | tab_montee_1.tab(:,3) < 0 ; %ignorer les descentes donc pitch > 0
 tab_montee_1.tab(ign,:)=[];

ign1 = abs(tab_montee_1.tab(:,8)) > 0.005; %ignorer les accélérations trop fortes 
tab_montee_1.tab(ign1,:)=[];

ign2=tab_montee_1.tab(:,1) > 500 | tab_montee_1.tab(:,1) < 100;
tab_montee_1.tab(ign2,:)=[];

%%%Optimisation %%%%%% 

global pres dens pitch oil_vol Wglider temp mg
pres=tab_montee_1.tab(:,1);
dens=tab_montee_1.tab(:,2);
pitch=tab_montee_1.tab(:,3);
oil_vol=tab_montee_1.tab(:,4);
Wglider=tab_montee_1.tab(:,5)';
temp=tab_montee_1.tab(:,6);
mg=explorer_all.M;

%%%Première optimisation
aw=3.7;
Cd1w = 0.78;
param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
[x_mont_1] = fminsearch('cost',param0,options); 
%precision_montee_1=sqrt(precision_montee_1/length(pres));
%%%Deuxième optimisation
global V0 eps Cd
V0=x_mont_1(1);
eps=x_mont_1(2);
Cd=x_mont_1(3);
param0=[aw,Cd1w];
[x_mont_11,precision_montee_11] = fminsearch('cost_wings',param0,options);

Precision_m1=[Precision_m1 precision_montee_11];
Quintuplet_mont_1 =[Quintuplet_mont_1 [x_mont_1 x_mont_11]'];
[W_model_mont_1,att_deg,U,att,Fg,Fb_m1,Fl,Fd,vx_m1] = flight_model2(pres,dens,pitch,oil_vol,temp,x_mont_1(1),x_mont_1(2),x_mont_1(3),mg,x_mont_11(1),x_mont_11(2));
Ww_mont_1 = tab_montee_1.tab(:,5)-W_model_mont_1;

[W_model_mont_1_moy] = flight_model2(pres,dens,pitch,oil_vol,temp,Const.V0,Const.eps_m,Const.Cd_m,mg,Const.aw,Const.Cdw);
Ww_mont_1_moy = tab_montee_1.tab(:,5)-W_model_mont_1_moy;
 end

   %%%%% Optimisation sur 2ème montée  %%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%---------------------------------------------------------------------%
 
 if opt_montee_2 == 1 
    
     %%%Préparation des données %%%%%% 
 
   tab_montee_2.tab=[];
   tab_montee_2.pressure=[];  
   tab_montee_2.dens=[];
   tab_montee_2.pitch=[];
   tab_montee_2.oil=[];
   tab_montee_2.w_glider=[];
   tab_montee_2.temp=[];
   tab_montee_2.time=[];
   tab_montee_2.acc=[];
 
     for k=2:2:explorer.last_dive-explorer.first_dive+1
         
       
        tab_montee_2.pressure=[tab_montee_2.pressure explorer_bydive.pressure(:,k)']; % On les remplit avec montée et descente une plongée sur deux 
        tab_montee_2.dens=[tab_montee_2.dens explorer_bydive.dens(:,k)'];
        tab_montee_2.pitch=[tab_montee_2.pitch explorer_bydive.pitch_filter(:,k)'];
        tab_montee_2.oil=[tab_montee_2.oil explorer_bydive.oil(:,k)'];
        tab_montee_2.w_glider=[tab_montee_2.w_glider explorer_bydive.W_glider_filter(:,k)'];
        tab_montee_2.temp=[tab_montee_2.temp explorer_bydive.temp(:,k)'];
        tab_montee_2.time=[tab_montee_2.time explorer_bydive.time(:,k)'];
        tab_montee_2.acc=[tab_montee_2.acc (explorer_bydive.W_glider_acc(:,k).*explorer_all.M)'];


     end
     
tab_montee_2.tab(:,1)=tab_montee_2.pressure; %Création du tableau
tab_montee_2.tab(:,2)=tab_montee_2.dens;
tab_montee_2.tab(:,3)=tab_montee_2.pitch;
tab_montee_2.tab(:,4)=tab_montee_2.oil;
tab_montee_2.tab(:,5)=tab_montee_2.w_glider;
tab_montee_2.tab(:,6)=tab_montee_2.temp;
tab_montee_2.tab(:,7)=tab_montee_2.time;
tab_montee_2.tab(:,8)=tab_montee_2.acc;

ign = tab_montee_2.tab(:,4) < 0 | tab_montee_2.tab(:,3) < 0 ; %ignorer les descentes donc pitch > 0
tab_montee_2.tab(ign,:)=[];

ign1 = abs(tab_montee_2.tab(:,8)) > 0.005; %ignorer les accélérations trop fortes 
tab_montee_2.tab(ign1,:)=[];

ign2=tab_montee_2.tab(:,1) > 500 | tab_montee_2.tab(:,1) < 100;
tab_montee_2.tab(ign2,:)=[];


%%%Optimisation %%%%%% 

global pres dens pitch oil_vol Wglider temp mg
pres=tab_montee_2.tab(:,1);
dens=tab_montee_2.tab(:,2);
pitch=tab_montee_2.tab(:,3);
oil_vol=tab_montee_2.tab(:,4);
Wglider=tab_montee_2.tab(:,5)';
temp=tab_montee_2.tab(:,6);
mg=explorer_all.M;

%%%Première optimisation
aw=3.7;
Cd1w = 0.78;
param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
[x_mont_2,precision_montee_2] = fminsearch('cost',param0,options); 
precision_montee_2=sqrt(precision_montee_2/length(pres));
%%%Deuxième optimisation
global V0 eps Cd
V0=x_mont_2(1);
eps=x_mont_2(2);
Cd=x_mont_2(3);
param0=[aw,Cd1w];
[x_mont_22,precision_montee_22] = fminsearch('cost_wings',param0,options);

Precision_m2=[Precision_m2 precision_montee_22];
Quintuplet_mont_2 =[Quintuplet_mont_2 [x_mont_2 x_mont_22]'];
[W_model_mont_2,att_deg,U,att,Fg,Fb_d1,Fl,Fd,vx_m2] = flight_model2(pres,dens,pitch,oil_vol,temp,x_mont_2(1),x_mont_2(2),x_mont_2(3),mg,x_mont_22(1),x_mont_22(2));

%histogramme
Ww_mont_2 = tab_montee_2.tab(:,5)-W_model_mont_2;

[W_model_mont_2_moy] = flight_model2(pres,dens,pitch,oil_vol,temp,Const.V0,Const.eps_m,Const.Cd_m,mg,Const.aw,Const.Cdw);
Ww_mont_2_moy = tab_montee_2.tab(:,5)-W_model_mont_2_moy;
 end

 
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------%


if W_glider_W_Model == 1 
    pause(1.5)

  if opt_bydive == 1       
  figure('Name','Wglider and Wmodel bydive','NumberTitle','off','Units','normalized','Position',[0.66,0.5,0.33,0.42]);
  plot(time_all, W_glider_all,'LineWidth',1.5)
  %plot( W_glider_all,'LineWidth',1.5)
  title('Vitesses verticales du glider by dive')
  hold on 
  %plot(W_model_all,'LineWidth',1.5)
  plot(time_all,W_model_all,'LineWidth',1.5)
 % datetick('x',0,'keepticks')
  legend('W\_glider','W\_model')
  pause(1)
  end
  
  if opt_all ==1
  figure('Name','Wglider and Wmodel all','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38])
 % plot(explorer_all.time, explorer_all.W_glider_filter,'LineWidth',1.5)
  plot(tab_all.tab(:,7), tab_all.tab(:,5),'LineWidth',1.5)
  title('Vitesses verticales du glider all')
  hold on 
  plot(tab_all.tab(:,7),W_model_all,'LineWidth',1.5)
  xticks(explorer_all.time(1):0.2:explorer_all.time(end))
  xticklabels(datestr(explorer_all.time(1):0.2:explorer_all.time(end)))
  %datetick('x',0,'keepticks')
  legend('W\_glider','W\_model')
  
  figure()
  histogram(Ww_all,30,'BinLimits',[-0.025 0.025])
 % xlim([-0.02 0.02])
  %ylim([0 200])
  title('Ww\_all')
  
   figure()
  histogram(att_deg,50,'BinLimits',[-5 5])
 % xlim([-0.02 0.02])
  %ylim([0 200])
  title('Attack angle')
  end
  
  if opt_diving ==1
  figure('Name','Wglider and Wmodel all test','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38])
  plot(tab_diving.tab(:,7), tab_diving.tab(:,5)','+','LineWidth',1)
  title('Vitesses verticales du glider diving')
  hold on 
  plot(tab_diving.tab(:,7),W_model_diving,'+','LineWidth',1)
  ylim([-0.25 0.1])
  datetick('x',0,'keepticks')
  legend('W\_glider','W\_model')
  
  figure()
  histogram(Ww_diving,30,'BinLimits',[-0.025 0.025])
  %xlim([-0.02 0.02])
  %ylim([0 200])
  title('Ww\_diving')
  end
  
  if opt_climbing ==1
  figure('Name','Wglider and Wmodel all test','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38])
  plot(tab_climbing.tab(:,7), tab_climbing.tab(:,5)','+','LineWidth',1)
  title('Vitesses verticales du glider climbing')
  hold on 
  plot(tab_climbing.tab(:,7),W_model_climbing,'+','LineWidth',1)
  ylim([0 0.3])
  datetick('x',0,'keepticks')
  legend('W\_glider','W\_model')
  
  figure()
  histogram(Ww_climbing,30,'BinLimits',[-0.025 0.025])
  title('Ww\_climbing')
  end
  
  if opt_descente_1 ==1
  figure('Name','Wglider and Wmodel descent_1','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38])
  plot(tab_descent_1.tab(:,7),tab_descent_1.tab(:,5)','+','LineWidth',1)
  title('Vitesses verticales du glider descent\_1')
  hold on 
  plot(tab_descent_1.tab(:,7),W_model_desc_1,'+','LineWidth',1)
  datetick('x',15)
  legend('W\_glider','W\_model')
  
  figure()
  histogram(Ww_desc_1,30,'BinLimits',[-0.02 0.02])
%   xlim([-0.02 0.02])
%   ylim([0 200])
  title('Ww\_desc\_1')
  end
  
  if opt_descente_2 ==1
  figure('Name','Wglider and Wmodel descent_2','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38])
  plot(tab_descent_2.tab(:,7),tab_descent_2.tab(:,5)','+','LineWidth',1)
  title('Vitesses verticales du glider descent\_2')
  hold on 
  plot(tab_descent_2.tab(:,7),W_model_desc_2,'+','LineWidth',1)
  ylim([-0.5 0.1])
  datetick('x',0)
  legend('W\_glider','W\_model')
  
  figure()
  histogram(Ww_desc_2,30,'BinLimits',[-0.02 0.02])
%    xlim([-0.02 0.02])
%   ylim([0 200])
  title('Ww\_desc\_2')
  end
  
  if opt_montee_1 ==1
  figure('Name','Wglider and Wmodel montee_1','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38])
  plot(tab_montee_1.tab(:,7),tab_montee_1.tab(:,5)','+','LineWidth',1)
  title('Vitesses verticales du glider montee\_1')
  hold on 
  plot(tab_montee_1.tab(:,7),W_model_mont_1,'+','LineWidth',1)
  datetick('x',15)
  legend('W\_glider','W\_model')
  
  figure()
  histogram(Ww_mont_1,30,'BinLimits',[-0.02 0.02])
%    xlim([-0.02 0.02])
%   ylim([0 200])
  title('Ww\_mont\_1')
  end
  
  if opt_montee_2 ==1
  figure('Name','Wglider and Wmodel montee_2','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38])
  plot(tab_montee_2.tab(:,7),tab_montee_2.tab(:,5)','+','LineWidth',1)
  title('Vitesses verticales du glider montee\_2')
  hold on 
  plot(tab_montee_2.tab(:,7),W_model_mont_2,'+','LineWidth',1)
  ylim([0 0.25])
  datetick('x',15)
  legend('W\_glider','W\_model')
  
  figure()
  histogram(Ww_mont_2,30,'BinLimits',[-0.02 0.02])
%    xlim([-0.02 0.02])
%   ylim([0 200])
  title('Ww\_mont\_2')
  end
end 


%% %% ----- Display Parameters evolution ----- %%
% ------------------------------------- %


if Parameters_evolution == 1
    pause(1.5)
 
    Triplet = [];   %%%%%évolution des paramètres
    for p=explorer.first:2:explorer.last    %%fenetre glissante 
        
        
        first = p;
        last = p+11;
        i_dive= tableau(:,1) >= first & tableau(:,1) <= last ;
        tab = tableau(i_dive,:);   

        %Read Data
        explorer_evolution = read_EXPLORER(tab,explorer,Const);
            
tab_all2.tab=[];
 
tab_all2.tab(:,1)=explorer.pressure; %Création du tableau
tab_all2.tab(:,2)=explorer.dens;
tab_all2.tab(:,3)=explorer.pitch_filter';
%tab_all.tab(:,3)=explorer_all.pitch; %non filtré
tab_all2.tab(:,4)=explorer.oil;
tab_all2.tab(:,5)=explorer.W_glider_filter';
%tab_all.tab(:,5)=explorer_all.W_glider;%non filtré
tab_all2.tab(:,6)=explorer.temp;
tab_all2.tab(:,7)=explorer.time;
tab_all2.tab(:,8)=explorer.W_glider_acc*explorer_all.M;

ign1 = abs(tab_all2.tab(:,8)) > 0.01; %filtre sur l'accélération 
tab_all2.tab(ign1,:)=[];


        global pres dens pitch oil_vol Wglider temp mg
        pres=tab_all2.tab(:,1);
        dens=tab_all2.tab(:,2);
        pitch=tab_all2.tab(:,3);
        oil_vol=tab_all2.tab(:,4);
        Wglider=tab_all2.tab(:,5)';
        temp=tab_all2.tab(:,6);
        mg=explorer_all.M;
   
        param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
        options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
        [x] = fminsearch('cost',param0,options);
    
    Triplet = [Triplet x'];

    end
    
    figure('Name','Parameters evolution','NumberTitle','off','Units','normalized','Position',[0.33,0.05,0.33,0.38]);
    subplot(3,1,1)
    plot([explorer.first:2:explorer.last],Triplet(3,:),'Color','#0072BD','LineWidth',2)
    title('Paramètres optimisés avec fenetre glissante de 3 plongées')
    xlabel('Numéro de plongée')
    ylabel('Cd')
    subplot(3,1,2)
    plot([explorer.first:2:explorer.last],Triplet(1,:),'Color','#D95319','LineWidth',2)
    xlabel('Numéro de plongée')
    ylabel('V0')
    subplot(3,1,3)
    plot([explorer.first:2:explorer.last],Triplet(2,:),'Color','#77AC30','LineWidth',2)
    xlabel('Numéro de plongée')
    ylabel('eps')

end

if Parameters_evolution_bydive == 1
    pause(1.5)
    
    figure('Name','Parameters evolution','NumberTitle','off','Units','normalized','Position',[0.66,0.05,0.33,0.38]);
    subplot(3,1,1)
    plot([explorer.first_dive:explorer.last_dive],Tri(3,:),'Color','#0072BD','LineWidth',2)
    title('Paramètres optimisés by dive')
    xlabel('Numéro de plongée')
    ylabel('Cd')
    subplot(3,1,2)
    plot([explorer.first_dive:explorer.last_dive],Tri(1,:),'Color','#D95319','LineWidth',2)
    xlabel('Numéro de plongée')
    ylabel('V0')
    subplot(3,1,3)
    plot([explorer.first_dive:explorer.last_dive],Tri(2,:),'Color','#77AC30','LineWidth',2)
    xlabel('Numéro de plongée')
    ylabel('eps')
    
end

%% %% ----- Display water vertical velocities ----- %%
% ------------------------------------- %



if water_velocity_descent_bydive == 1 
    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww1=[];
for k=1:explorer.last_dive-explorer.first_dive+1 
      
   global pres dens pitch oil_vol Wglider temp mg
  pres=explorer3_bydive.pressure(:,k);
  dens=explorer3_bydive.dens(:,k);
  pitch=explorer3_bydive.pitch_filter(:,k);
  oil_vol=explorer3_bydive.oil(:,k);
  Wglider=explorer3_bydive.W_glider_filter(:,k)';%prime avec flightmodel2
  temp=explorer3_bydive.temp(:,k);
  mg=explorer_all.M;

  % [W_model] = flight_model(pres,dens,pitch,oil_vol,temp,Tri(1,k),Tri(2,k),Tri(3,k),mg);
   [W_model] = flight_model2(pres,dens,pitch,oil_vol,temp,Tri(1,k),Tri(2,k),Tri(3,k),mg,Tri(4,k),Tri(5,k));
    %W_model = fillmissing(W_model,'previous');
 
      explorer2_bydive.temp=explorer3_bydive.temp(:,k);
      explorer2_bydive.dens=explorer3_bydive.dens(:,k);
      explorer2_bydive.W_glider_filter=explorer3_bydive.W_glider_filter(:,k);
       explorer2_bydive.time=explorer3_bydive.time(:,k);
       explorer2_bydive.pressure=explorer3_bydive.pressure(:,k);
       explorer2_bydive.s=explorer3_bydive.s(:,k);

%Select descent
 max_ind = find(explorer_bydive.pressure(:,k) == max(explorer_bydive.pressure(:,k)));
 to_ign = [max_ind:length(explorer2_bydive.pressure)];
 explorer2_bydive.p_sorted = explorer2_bydive.pressure;
 explorer2_bydive.p_sorted(to_ign)=[];
 explorer2_bydive.temp(to_ign)=[];
 explorer2_bydive.dens(to_ign)=[];
 explorer2_bydive.s(to_ign)=[];
 explorer2_bydive.W_glider_filter(to_ign)=[];
 W_model(to_ign)=[];
 explorer2_bydive.time(to_ign)=[];
 
 E=NaN(length(explorer2_bydive.p_sorted),7);
 E(:,1)= explorer2_bydive.p_sorted;
 E(:,2)= explorer2_bydive.temp;
 E(:,3)= explorer2_bydive.s;
 E(:,4)= explorer2_bydive.dens;
 E(:,5)= explorer2_bydive.W_glider_filter;
 E(:,6)= W_model;
 E(:,7)= explorer2_bydive.time;
 
 F=sortrows(E,1);
 
   to_igno = [];
        for i=1:length(F)-1
         if  F(i,1) == F(i+1,1) % delete doublon
           to_igno = [to_igno i+1];
         end
        end
        F(to_igno,:)=[];

 ign_nan = isnan(F(:,1));
 F(ign_nan,:)=[];
 
 
    
    Ww=F(:,5)-F(:,6);
  
    
    pi=[0:1:max(explorer_all.pressure)];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
    ti=interp1(F(:,1),F(:,2),pi); 
    si=interp1(F(:,1),F(:,3),pi);
    di=interp1(F(:,1),F(:,4),pi);
    wg=interp1(F(:,1),F(:,5),pi);
    wm=interp1(F(:,1),F(:,6),pi);
    ww=interp1(F(:,1),Ww,pi);
    
         %Filter on ww 
    ind = find(abs(ww) > 0.02);
    ww(ind)=NaN;
   
    tableau_ti=[tableau_ti ti'];
    tableau_si=[tableau_si si'];
    tableau_di=[tableau_di di'];
    tableau_w_glider=[tableau_w_glider wg'];
    tableau_w_model=[tableau_w_model wm'];
    tableau_ww1=[tableau_ww1 ww'];
    
    
end

figure('Name','Water vertical velocity descent','NumberTitle','off','Units','normalized','Position',[0,0.5,0.33,0.42]);

if sbplot == 1

    subplot(2,1,1)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ti);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'       (°C)','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Température')

  
    subplot(2,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww1);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww')
else 
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww1);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    colormap(blue_red_cmap)
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_bydive\_descent')
    
end
if hist_desc == 1 
    pause(1.5)
    figure('Name','Velocity','NumberTitle','off','Units','normalized','Position',[0.33,0.5,0.33,0.42]);
    histogram(tableau_ww1,100)
    title('velocities')
    xlabel('velocities m/s')
    %xlim([-4 5])

end
end

if water_velocity_descent_all == 1 


  
    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww2=[];
e=0;
for j= explorer.first_dive:explorer.last_dive 
     if j==explorer.first + e 
       explorer = by_dive(tableau,j,explorer,Const,'d1');%data by dive
     else 
       explorer = by_dive(tableau,j,explorer,Const,'d2');%data by dive
          e=e+2;
     end

       
   global pres dens pitch oil_vol Wglider temp mg
   pres=explorer.pressure;
   dens=explorer.dens;
   pitch=explorer.pitch_filter;
   oil_vol=explorer.oil;
   Wglider=explorer.W_glider;
   temp=explorer.temp;
   mg=explorer_all.M;

   %[W_model] = flight_model(pres,dens,pitch,oil_vol,temp,x_all(1),x_all(2),x_all(3),mg);
   [W_model] = flight_model(pres,dens,pitch,oil_vol,temp,0.059945,4.0265e-10,0.1040,mg);

 max_ind = find(explorer.pressure == max(explorer.pressure));
 to_ign = [max_ind:length(explorer.pressure)];
 explorer.p_sorted = explorer.pressure;
 explorer.p_sorted(to_ign)=[];
 explorer.temp(to_ign)=[];
 explorer.s(to_ign)=[];
 explorer.dens(to_ign)=[];
 explorer.W_glider(to_ign)=[];
 W_model(to_ign)=[];
 explorer.time(to_ign)=[];
 
 E=NaN(length(explorer.p_sorted),7);
 E(:,1)= explorer.p_sorted;
 E(:,2)= explorer.temp;
 E(:,3)= explorer.s;
 E(:,4)= explorer.dens;
 E(:,5)= explorer.W_glider';
 E(:,6)= W_model;
 E(:,7)= explorer.time;
 
 F=sortrows(E,1);
 q=size(F);
   to_igno = [];
        for i=1:q(1)-1
         if  F(i,1) == F(i+1,1) % delete doublon
           to_igno = [to_igno i+1];
         end
        end
        F(to_igno,:)=[];
   
    Ww=F(:,5)-F(:,6);
  
  
    pi=[0:1:640];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
    ti=interp1(F(:,1),F(:,2),pi); 
    si=interp1(F(:,1),F(:,3),pi);
    di=interp1(F(:,1),F(:,4),pi);
    wg=interp1(F(:,1),F(:,5),pi);
    wm=interp1(F(:,1),F(:,6),pi);
    ww=interp1(F(:,1),Ww,pi);
    
         %Filter on ww 
    ind = find(abs(ww) > 0.03);
    ww(ind)=NaN;
   
    tableau_ti=[tableau_ti ti'];
    tableau_si=[tableau_si si'];
    tableau_di=[tableau_di di'];
    tableau_w_glider=[tableau_w_glider wg'];
    tableau_w_model=[tableau_w_model wm'];
    tableau_ww2=[tableau_ww2 ww'];
    
end
figure('Name','Water vertical velocity descent','NumberTitle','off','Units','normalized','Position',[0.66,0.5,0.33,0.42]);

if sbplot == 1

    subplot(2,1,1)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ti);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'       (°C)','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Température')

   
    subplot(2,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww2);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww')
else 
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww2);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
     colormap(blue_red_cmap)
    H=colorbar;
    caxis([-0.025 0.025])
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_all\_descent')
    
end
if hist_desc == 1 
    pause(1.5)
    figure('Name','Velocity','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38]);
    histogram(tableau_ww2,100)
    title('velocities')
    xlabel('velocities m/s')
    %xlim([-4 5])

end
end
%-----------------------------------------------------------------------%
if water_velocity_diving == 1 


  
    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww6=[];
for j= explorer.first_dive:explorer.last_dive 
       explorer = by_dive(tableau,j,explorer,Const,'d1');%data by dive
 %%% Opt on diving 
sz=length(explorer.pressure);
tab1=zeros(sz,7);
tab1(:,1)=explorer.pressure;
tab1(:,2)=explorer.dens;
tab1(:,3)=explorer.pitch;
tab1(:,4)=explorer.oil;
tab1(:,5)=explorer.W_glider;
tab1(:,6)=explorer.temp;
tab1(:,7)=explorer.time;


ign = tab1(:,3) > 0;
tab1(ign,:)=[];
Mean_glider= mean(tab1(:,5));
ign1=tab1(:,5) > Mean_glider +0.05;
tab1(ign1,:)=[];
ign2=tab1(:,1) > 500 | tab1(:,1) < 100;
tab1(ign2,:)=[];

     global pres dens pitch oil_vol Wglider temp mg
        pres=tab1(:,1);
        dens=tab1(:,2);
        pitch=tab1(:,3);
        oil_vol=tab1(:,4);
        Wglider=tab1(:,5)';
        temp=tab1(:,6);
        mg=explorer_all.M;

   [W_model] = flight_model(pres,dens,pitch,oil_vol,temp,test(1),test(2),test(3),mg);


 G=sortrows(tab1,1);
 
   to_igno = [];
        for i=1:length(G)-1
         if  G(i,1) == G(i+1,1) % delete doublon
           to_igno = [to_igno i+1];
         end
        end
        G(to_igno,:)=[];
  
    
    Ww_test=G(:,5)-W_model;
  
  
    
    pi=[0:1:600];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
    ti=interp1(G(:,1),G(:,2),pi); 
    si=interp1(G(:,1),G(:,3),pi);
    di=interp1(G(:,1),G(:,4),pi);
    wg=interp1(G(:,1),G(:,5),pi);
    wm=interp1(G(:,1),G(:,6),pi);
    ww_test=interp1(G(:,1),Ww_test,pi);
    
         %Filter on ww 
    ind = find(abs(ww_test) > 0.03);
    ww_test(ind)=NaN;
   
    tableau_ti=[tableau_ti ti'];
    tableau_si=[tableau_si si'];
    tableau_di=[tableau_di di'];
    tableau_w_glider=[tableau_w_glider wg'];
    tableau_w_model=[tableau_w_model wm'];
    tableau_ww6=[tableau_ww6 ww_test'];
    

end
figure('Name','Water vertical velocity diving','NumberTitle','off','Units','normalized','Position',[0.66,0.5,0.33,0.42]);

if sbplot == 1

    subplot(2,1,1)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ti);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'       (°C)','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Température')

   
    subplot(2,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww2);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww')
else 
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww6);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    colormap(blue_red_cmap)
    H=colorbar;
    caxis([-0.025 0.025])
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_diving')
    
end
if hist_desc == 1 
    pause(1.5)
    figure('Name','Velocity','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38]);
    histogram(tableau_ww6,100)
    title('velocities diving')
    xlabel('velocities m/s')
    %xlim([-4 5])

end
end
%-----------------------------------------------------------------------%
if water_velocity_ascent_all == 1 

   
    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww3=[];
for j= explorer.first_dive:explorer.last_dive 
       explorer = by_dive(tableau,j,explorer,Const,'d1');%data by dive
       
   
   global pres dens pitch oil_vol Wglider temp mg
   pres=explorer.pressure;
   dens=explorer.dens;
   pitch=explorer.pitch;
   oil_vol=explorer.oil;
   Wglider=explorer.W_glider;
   temp=explorer.temp;
   mg=explorer_all.M;

%   [W_model] = flight_model(pres,dens,pitch,oil_vol,temp,x_all(1),x_all(2),x_all(3),mg);
   [W_model] = flight_model(pres,dens,pitch,oil_vol,temp,0.059945,4.0265e-10,0.995,mg);


 max_ind = find(explorer.pressure == max(explorer.pressure));
 to_ign = [1:max_ind];
 explorer.p_sorted = explorer.pressure;
 explorer.p_sorted(to_ign)=[];
 explorer.temp(to_ign)=[];
 explorer.s(to_ign)=[];
 explorer.dens(to_ign)=[];
 explorer.W_glider(to_ign)=[];
 W_model(to_ign)=[];
 explorer.time(to_ign)=[];
 
 E=NaN(length(explorer.p_sorted),7);
 E(:,1)= explorer.p_sorted;
 E(:,2)= explorer.temp;
 E(:,3)= explorer.s;
 E(:,4)= explorer.dens;
 E(:,5)= explorer.W_glider';
 E(:,6)= W_model;
 E(:,7)= explorer.time;
 
 F=sortrows(E,1);
 q=size(F);
   to_igno = [];
        for i=1:q(1)-1
         if  F(i,1) == F(i+1,1) % delete doublon
           to_igno = [to_igno i+1];
         end
        end
        F(to_igno,:)=[];
  
    Ww=F(:,5)-F(:,6);
  
  
    pi=[0:1:640];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
    ti=interp1(F(:,1),F(:,2),pi); 
    si=interp1(F(:,1),F(:,3),pi);
    di=interp1(F(:,1),F(:,4),pi);
    wg=interp1(F(:,1),F(:,5),pi);
    wm=interp1(F(:,1),F(:,6),pi);
    ww=interp1(F(:,1),Ww,pi);
    
         %Filter on ww 
    ind = find(abs(ww) > 0.03);
    ww(ind)=NaN;
   
    tableau_ti=[tableau_ti ti'];
    tableau_si=[tableau_si si'];
    tableau_di=[tableau_di di'];
    tableau_w_glider=[tableau_w_glider wg'];
    tableau_w_model=[tableau_w_model wm'];
    tableau_ww3=[tableau_ww3 ww'];
    
end

 figure('Name','Water vertical velocity ascent','NumberTitle','off','Units','normalized','Position',[0.33,0.05,0.33,0.38]);

 if sbplot == 1
 
     subplot(2,1,1)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ti);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'       (°C)','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Température')

  
    subplot(2,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww3);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
 
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww')
 else
     pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww3);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    colormap(blue_red_cmap)
    H=colorbar;
    caxis([-0.025 0.025])
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_all\_ascent')
 end
 
 if hist_asc == 1 
    pause(1.5)
    figure('Name','Velocity','NumberTitle','off','Units','normalized','Position',[0.66,0.05,0.33,0.38]);
    histogram(tableau_ww3,100)
    title('velocities')
    xlabel('velocities m/s')
    %xlim([-4 5])

end
 
end

if water_velocity_ascent_by_dive == 1 

    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww4=[];
for k=1:explorer.last_dive-explorer.first_dive+1 
       
   
   global pres dens pitch oil_vol Wglider temp mg
  pres=explorer_bydive.pressure(:,k);
  dens=explorer_bydive.dens(:,k);
  pitch=explorer_bydive.pitch_filter(:,k);
  oil_vol=explorer_bydive.oil(:,k);
  Wglider=explorer_bydive.W_glider_filter(:,k);
  temp=explorer_bydive.temp(:,k);
  mg=explorer_all.M;

   [W_model] = flight_model(pres,dens,pitch,oil_vol,temp,Tri(1,k),Tri(2,k),Tri(3,k),mg);
   W_model = fillmissing(W_model,'previous');
      %explorer2_bydive=explorer_bydive;
      explorer2_bydive.temp=explorer_bydive.temp(:,k);
      explorer2_bydive.s=explorer_bydive.s(:,k);
      explorer2_bydive.dens=explorer_bydive.dens(:,k);
      explorer2_bydive.W_glider_filter=explorer_bydive.W_glider_filter(:,k);
       explorer2_bydive.time=explorer_bydive.time(:,k);
       explorer2_bydive.pressure=explorer_bydive.pressure(:,k);
       
 max_ind = find(explorer2_bydive.pressure == max(explorer2_bydive.pressure));
 to_ign = [1:max_ind];
 explorer2_bydive.p_sorted = explorer2_bydive.pressure;
 explorer2_bydive.p_sorted(to_ign)=[];
 explorer2_bydive.temp(to_ign)=[];
 explorer2_bydive.s(to_ign)=[];
 explorer2_bydive.dens(to_ign)=[];
 explorer2_bydive.W_glider_filter(to_ign)=[];
 W_model(to_ign)=[];
 explorer2_bydive.time(to_ign)=[];
 
 E=NaN(length(explorer2_bydive.p_sorted),7);
 E(:,1)= explorer2_bydive.p_sorted;
 E(:,2)= explorer2_bydive.temp;
 E(:,3)= explorer2_bydive.s;
 E(:,4)= explorer2_bydive.dens;
 E(:,5)= explorer2_bydive.W_glider_filter;
 E(:,6)= W_model;
 E(:,7)= explorer2_bydive.time;
 
 F=sortrows(E,1);
 
   to_igno = [];
        for i=1:length(F)-1
         if  F(i,1) == F(i+1,1) % delete doublon
           to_igno = [to_igno i+1];
         end
        end
        F(to_igno,:)=[];
        
    %remove NaN 
    igno=isnan(F(:,1));
    F(igno,:)=[];
  
  
%        
%        % remontada
%         to_ign=[];
%         max_ind = find(explorer2_bydive.pressure == max(explorer2_bydive.pressure));
%         i=max_ind;
%         while i ~= length(explorer2_bydive.pressure)
%           if  not(explorer2_bydive.pressure(i)>explorer2_bydive.pressure(i+1)) % delete variations
%             to_ign = [to_ign i+1];
%           end
%           i=i+1;
%         end
% 
%        to_ign = [to_ign 1:max_ind]; % taking only remontada 
%        explorer2_bydive.p_sorted = explorer2_bydive.pressure;
%        explorer2_bydive.p_sorted(to_ign)=[];
%        explorer2_bydive.temp(to_ign)=[];
%        explorer2_bydive.s(to_ign)=[];
%        explorer2_bydive.dens(to_ign)=[];
%        explorer2_bydive.W_glider_filter(to_ign)=[];
%        W_model(to_ign)=[];
%        explorer2_bydive.time(to_ign)=[];
% 
% 
%         %remontada
%         ind1=[];
%         o=length(explorer2_bydive.p_sorted);
%         while explorer2_bydive.p_sorted(o) < explorer2_bydive.p_sorted(o-1) && o > 2
%           % delete variations2
%             ind1 = [ind1 o];
%           o=o-1;
%         end
%         %remontada
%         explorer2_bydive.p_sorted=explorer2_bydive.p_sorted(ind1);
%         explorer2_bydive.temp=explorer2_bydive.temp(ind1);
%         explorer2_bydive.s=explorer2_bydive.s(ind1);
%         explorer2_bydive.dens=explorer2_bydive.dens(ind1);
%         explorer2_bydive.W_glider_filter=explorer2_bydive.W_glider_filter(ind1);
%         W_model=W_model(ind1);
%         explorer2_bydive.time=explorer2_bydive.time(ind1);
% 
%     to_igno = [];
%     for i=1:length(explorer2_bydive.p_sorted)-1
%      if  explorer2_bydive.p_sorted(i) == explorer2_bydive.p_sorted(i+1) % delete doublon
%        to_igno = [to_igno i+1];
%      end
%     end
%     
%     explorer2_bydive.p_sorted(to_igno)=[];
%     explorer2_bydive.temp(to_igno)=[];
%     explorer2_bydive.s(to_igno)=[];
%     explorer2_bydive.dens(to_igno)=[];
%     explorer2_bydive.W_glider_filter(to_igno)=[];
%     W_model(to_igno)=[];
%     explorer2_bydive.time(to_igno)=[];
    
    Ww= F(:,5)- F(:,6);
  
  
    
    pi=[0:1:600];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
    ti=interp1(F(:,1),F(:,2),pi); 
    si=interp1(F(:,1),F(:,3),pi);
    di=interp1(F(:,1),F(:,4),pi);
    wg=interp1(F(:,1),F(:,5),pi);
    wm=interp1(F(:,1),F(:,6),pi);
    ww=interp1(F(:,1),Ww,pi);
    
         %Filter on ww 
    ind = find(abs(ww) > 0.02);
    ww(ind)=NaN;
   
    tableau_ti=[tableau_ti ti'];
    tableau_si=[tableau_si si'];
    tableau_di=[tableau_di di'];
    tableau_w_glider=[tableau_w_glider wg'];
    tableau_w_model=[tableau_w_model wm'];
    tableau_ww4=[tableau_ww4 ww'];
    
%     figure()
%     plot(explorer.time,explorer.W_glider)
%     hold on 
%     plot(explorer.time,W_model)
    
end

 figure('Name','Water vertical velocity ascent','NumberTitle','off','Units','normalized','Position',[0,0.5,0.33,0.42]);

 if sbplot == 1
 
     subplot(2,1,1)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ti);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'       (°C)','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Température')

    
    subplot(2,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww4);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww')
 else
     pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww4);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    colormap(blue_red_cmap)
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_bydive\_ascent')
 end
 
 if hist_asc == 1 
    pause(1.5)
    figure('Name','Velocity','NumberTitle','off','Units','normalized','Position',[0.33,0.5,0.33,0.42]);
    histogram(tableau_ww4,100)
    title('velocities')
    xlabel('velocities m/s')
    %xlim([-4 5])

end
 
end

if sub_asc_desc == 1 
    
figure()
subplot(2,1,1)
 pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww1);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    colormap(blue_red_cmap)
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_descent')
    
subplot(2,1,2)
  pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww4);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    colormap(blue_red_cmap)
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_ascent')
    
end

%------------------------------------------------------------------------%
if water_velocity_all_bydive == 1 


tableau_ww5=zeros(length(tableau_ww1),min(size(tableau_ww1)));
a=1;
b=1;
for i =1:min(size(tableau_ww1))*2
    
    if mod(i,2)~=0 % i impair
        tableau_ww5(:,i)=tableau_ww1(:,a);
        a=a+1;
    else
        tableau_ww5(:,i)=tableau_ww4(:,b);
        b=b+1;
    end
end

figure('Name','Water vertical velocity ','NumberTitle','off','Units','normalized','Position',[0.66,0.5,0.33,0.42])
    pcolor([explorer.first_dive:(explorer.last_dive-explorer.first_dive+1)*2+explorer.first_dive-1],-pi,tableau_ww5);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    colormap(blue_red_cmap)
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_all')



end


%% %% ----- Descente_1 ----- %%
% ------------------------------------- %
if Descente_1==1
if valeurs_opt == 1 
    
q=size(tab_descent_1.tab);
Ww_d1=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
while i < q(1)+1
    if i == 1 
        Ww_1(j)=Ww_desc_1(i);
        pression(j)=tab_descent_1.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if isnan(tab_descent_1.tab(i,1)) || abs(tab_descent_1.tab(i-1,1)-tab_descent_1.tab(i,1)) > max(tab_descent_1.tab(:,1))/2 || i == q(1)
            B=zeros(length(pression),2);
            B(:,1)=pression;
            B(:,2)=Ww_1;
            C=sortrows(B,1);
            to_igno = [];
            for k=1:length(C)-1
                if  C(k,1) == C(k+1,1) % delete doublon
                to_igno = [to_igno k+1];
                end
            end
            C(to_igno,:)=[];
            
           ww_1=interp1(C(:,1),C(:,2),pi);
            disp(i)
           
           Ww_d1=[Ww_d1 ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           while (isnan(tab_descent_1.tab(i,1)) || tab_descent_1.tab(i,1)<tab_descent_1.tab(i-1,1)) && i < q(1)
              i=i+1; 
           %   disp(i)
           end
           while i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_desc_1(i);
           pression(j)=tab_descent_1.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
    
end
Ww_d1_final =[Ww_d1_final Ww_d1]; 

figure()
pcolor([1:1:6],-pi,Ww_d1);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
 caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_d1')
end
    if valeurs_moyennes == 1
        q=size(tab_descent_1.tab);
Ww_d1=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
while i < q(1)+1
    if i == 1 
        Ww_1(j)=Ww_desc_1_moy(i);
        pression(j)=tab_descent_1.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if isnan(tab_descent_1.tab(i,1)) || tab_descent_1.tab(i,1)<tab_descent_1.tab(i-1,1) || i == q(1)
            disp(i)
           ww_1=interp1(pression,Ww_1,pi);
           Ww_d1=[Ww_d1 ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           while (isnan(tab_descent_1.tab(i,1)) || tab_descent_1.tab(i,1)<tab_descent_1.tab(i-1,1)) && i < q(1)
              i=i+1; 
           %   disp(i)
           end
           while i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_desc_1_moy(i);
           pression(j)=tab_descent_1.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
    
end

figure()
pcolor([1:1:24],-pi,Ww_d1);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
 caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_d1')
        
    end
end
%% %% ----- Montee_1 ----- %%
% ------------------------------------- %
if Montee_1 ==1

    if valeurs_opt == 1
q=size(tab_montee_1.tab);
Ww_m1=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
Ww_1=[];
pression=[];
while i < q(1)+1
    if i == 1 
        Ww_1(j)=Ww_mont_1(i);
        pression(j)=tab_montee_1.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if isnan(tab_montee_1.tab(i,1))|| abs(tab_montee_1.tab(i-1,1)-tab_montee_1.tab(i,1)) > max(tab_montee_1.tab(:,1))/2 || i == q(1)
            
           B=zeros(length(pression),2);
            B(:,1)=pression;
            B(:,2)=Ww_1;
            C=sortrows(B,1);
            to_igno = [];
            for k=1:length(C)-1
                if  C(k,1) == C(k+1,1) % delete doublon
                to_igno = [to_igno k+1];
                end
            end
            C(to_igno,:)=[];
            
           ww_1=interp1(C(:,1),C(:,2),pi);
           
           Ww_m1=[Ww_m1 ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           while (isnan(tab_montee_1.tab(i,1)) || tab_montee_1.tab(i,1)>tab_montee_1.tab(i-1,1)) && i < q(1)
              i=i+1; 
           end
            while i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_mont_1(i);
           pression(j)=tab_montee_1.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
end

Ww_m1_final =[Ww_m1_final Ww_m1];

figure()
pcolor([1:1:6],-pi,Ww_m1);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_m1')
    end
    
   if valeurs_moyennes == 1
       
       q=size(tab_montee_1.tab);
Ww_m1=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
Ww_1=[];
pression=[];
while i < q(1)+1
    if i == 1 
        Ww_1(j)=Ww_mont_1_moy(i);
        pression(j)=tab_montee_1.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if isnan(tab_montee_1.tab(i,1)) || tab_montee_1.tab(i,1)>tab_montee_1.tab(i-1,1) || i == q(1)

           ww_1=interp1(pression,Ww_1,pi);
           Ww_m1=[Ww_m1 ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           while (isnan(tab_montee_1.tab(i,1)) || tab_montee_1.tab(i,1)>tab_montee_1.tab(i-1,1)) && i < q(1)
              i=i+1; 
           end
            while i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_mont_1_moy(i);
           pression(j)=tab_montee_1.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
end

figure()
pcolor([1:1:24],-pi,Ww_m1);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_m1')
   end
end

%% %% ----- Descente_2 ----- %%
% ------------------------------------- %
if Descente_2==1
    if valeurs_opt == 1
q=size(tab_descent_2.tab);
Ww_d2=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
Ww_1=[];
pression=[];
while i < q(1)+1
    if i == 1 
        Ww_1(j)=Ww_desc_2(i);
        pression(j)=tab_descent_2.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if isnan(tab_descent_2.tab(i,1))|| abs(tab_descent_2.tab(i-1,1)-tab_descent_2.tab(i,1)) > max(tab_descent_2.tab(:,1))/2  || i == q(1)
            B=zeros(length(pression),2);
            B(:,1)=pression;
            B(:,2)=Ww_1;
            C=sortrows(B,1);
            to_igno = [];
            for k=1:length(C)-1
                if  C(k,1) == C(k+1,1) % delete doublon
                to_igno = [to_igno k+1];
                end
            end
            C(to_igno,:)=[];
            
           ww_1=interp1(C(:,1),C(:,2),pi);
         
           Ww_d2=[Ww_d2 ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           while (isnan(tab_descent_2.tab(i,1)) || tab_descent_2.tab(i,1)<tab_descent_2.tab(i-1,1)) && i < q(1)
              i=i+1; 
           end
           while i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_desc_2(i);
           pression(j)=tab_descent_2.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
end

Ww_d2_final =[Ww_d2_final Ww_d2];

figure()
pcolor([1:1:6],-pi,Ww_d2);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_d2')
    end
    
    if valeurs_moyennes == 1 
        q=size(tab_descent_2.tab);
Ww_d2=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
Ww_1=[];
pression=[];
while i < q(1)+1
    if i == 1 
        Ww_1(j)=Ww_desc_2_moy(i);
        pression(j)=tab_descent_2.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if isnan(tab_descent_2.tab(i,1)) || tab_descent_2.tab(i,1)<tab_descent_2.tab(i-1,1) || i == q(1)

           ww_1=interp1(pression,Ww_1,pi);
           Ww_d2=[Ww_d2 ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           while (isnan(tab_descent_2.tab(i,1)) || tab_descent_2.tab(i,1)<tab_descent_2.tab(i-1,1)) && i < q(1)
              i=i+1; 
           end
           while i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_desc_2_moy(i);
           pression(j)=tab_descent_2.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
end

figure()
pcolor([1:1:24],-pi,Ww_d2);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_d2')
    end
end
%% %% ----- Montee_2----- %%
% ------------------------------------- %
if Montee_2==1
    if valeurs_opt == 1 
q=size(tab_montee_2.tab);
Ww_m2=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
Ww_1=[];
pression=[];
while i < q(1)+1
    if i == 1 
        Ww_1(j)=Ww_mont_2(i);
        pression(j)=tab_montee_2.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if isnan(tab_montee_2.tab(i,1)) || abs(tab_montee_2.tab(i-1,1)-tab_montee_2.tab(i,1)) > max(tab_montee_2.tab(:,1))/2 || i == q(1)
            
           B=zeros(length(pression),2);
            B(:,1)=pression;
            B(:,2)=Ww_1;
            C=sortrows(B,1);
            to_igno = [];
            for k=1:length(C)-1
                if  C(k,1) == C(k+1,1) % delete doublon
                to_igno = [to_igno k+1];
                end
            end
            C(to_igno,:)=[];
            
           ww_1=interp1(C(:,1),C(:,2),pi);
           
           Ww_m2=[Ww_m2 ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           while (isnan(tab_montee_2.tab(i,1)) || tab_montee_2.tab(i,1)>tab_montee_2.tab(i-1,1)) && i < q(1)
              i=i+1; 
           end
           while i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_mont_2(i);
           pression(j)=tab_montee_2.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
end

Ww_m2_final =[Ww_m2_final Ww_m2];

figure()
pcolor([1:1:6],-pi,Ww_m2);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_m2')
    end
    
    if valeurs_moyennes ==1
        q=size(tab_montee_2.tab);
Ww_m2=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
Ww_1=[];
pression=[];
while i < q(1)+1
    if i == 1 
        Ww_1(j)=Ww_mont_2_moy(i);
        pression(j)=tab_montee_2.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if isnan(tab_montee_2.tab(i,1)) || tab_montee_2.tab(i,1)>tab_montee_2.tab(i-1,1) || i == q(1)

           ww_1=interp1(pression,Ww_1,pi);
           Ww_m2=[Ww_m2 ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           while (isnan(tab_montee_2.tab(i,1)) || tab_montee_2.tab(i,1)>tab_montee_2.tab(i-1,1)) && i < q(1)
              i=i+1; 
           end
           while i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_mont_2_moy(i);
           pression(j)=tab_montee_2.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
end

figure()
pcolor([1:1:24],-pi,Ww_m2);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_m2')
    end
end

end %fenetre glissante 
%% %% ----- Together 

if Together ==1
q=size(Ww_d1_final);
Ww_4parts=zeros(q(1),q(2)*4);
a=1;
b=1;
for i =1:q(2)*4
    
    if i == b 
        Ww_4parts(:,i)=Ww_d1_final(:,a);
        
        
    end
    if i == b+1
         Ww_4parts(:,i)=Ww_m1_final(:,a);
       
    end
    if i == b+2
         Ww_4parts(:,i)=Ww_d2_final(:,a);
        
    end
    if i == b+3
         Ww_4parts(:,i)=Ww_m2_final(:,a);
        b=b+4;
        a=a+1;
    end
end
 ind = find(abs(Ww_4parts) > 0.018);
 Ww_4parts(ind)=0;
 size_Ww_4parts=size(Ww_4parts);
 
 %xlabel 
 debut = 7.375465128009259e+05; 
fin = 7.375496856018519e+05;

A=linspace(debut,fin,7);
B=datestr(A);
 
figure('Name','Water vertical velocity ','NumberTitle','off','Units','normalized','Position',[0.66,0.5,0.33,0.42])
    pcolor([1:1:size_Ww_4parts(2)],-pi,Ww_4parts); % utilisation de la fonction pcolor + "shading interp"
    shading interp
    colormap(blue_red_cmap)
    H=colorbar;
    caxis([-0.025 0.025])
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xticklabels({B})
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_together')



end
%% %% ----- Display speed up and speed down ----- %%
% ------------------------------------- %
% speed_up_tot = [];
% speed_down_tot = [];
% explorer.time_tot=[];
% for j= explorer.first_dive:explorer.last_dive %data by dive
%     i_dive1 = tableau(:,1) == j;
%     tableau_dive = tableau(i_dive1,:);
%     explorer.depth_i = tableau_dive(:,9);
%     explorer.time_i = tableau_dive(:,8);
%     explorer.time_tot = [explorer.time_tot explorer.time_i(1:end-1)'];
%   
%     ind = find(explorer.depth_i == max(explorer.depth_i));
%     speed_down = zeros(1,length(explorer.depth_i)-1);
%     speed_up = zeros(1,length(explorer.depth_i)-1);
%     for k =1:length(explorer.depth_i)-1
%         if k < ind 
%         speed_down(k) = (explorer.depth_i(k)-explorer.depth_i(k+1))/(explorer.time_i(k+1)-explorer.time_i(k));
%         speed_up(k)=NaN;
%         else
%         speed_up(k) = (explorer.depth_i(k)-explorer.depth_i(k+1))/(explorer.time_i(k+1)-explorer.time_i(k));
%         speed_down(k)=NaN;
%         end
%     end
%     speed_down = speed_down./Const.d2s;
%     speed_up = speed_up./Const.d2s;
%     speed_down_tot = [speed_down_tot speed_down];
%     speed_up_tot = [speed_up_tot speed_up];
% end
% 
% figure()
% 
% plot(explorer.time_tot,speed_up_tot)
% hold on
% plot(explorer.time_tot,speed_down_tot)
% title('speed')
% datetick('x',0,'keepticks')
% legend('speed\_up','speed\_down')


%% %% ----- Display pitch and speed ----- %%
% ------------------------------------- %
% % 
% 
% tableau_pitch=[];
% tableau_W_glider=[];
% b1=[];
% for j= explorer.first_dive:explorer.last_dive 
%        explorer = by_dive(tableau,j);%data by dive
%     
%      pi=[0:1:600];
%      ind = find(explorer.pressure == max(explorer.pressure));
%      explorer.pressure=explorer.pressure(50:ind-10);
%      explorer.pitch=explorer.pitch(50:ind-10);
%      explorer.time=explorer.time(50:ind-10);
% 
%     W_glider = zeros(1,length(explorer.pressure)-10);
%     for k=1:length(explorer.pressure)-10
%         W_glider(k) = (explorer.pressure(k)-explorer.pressure(k+10))/(explorer.time(k+10)-explorer.time(k));
%     end
%     W_glider = W_glider./Const.d2s;
%     
%     
%      pitch=interp1(explorer.pressure,explorer.pitch,pi);
%      W_glider = interp1(explorer.pressure(1:end-10),W_glider,pi);
%      tableau_pitch=[tableau_pitch pitch'];
%      tableau_W_glider=[tableau_W_glider W_glider'];
%      
% end
% 
% for i=1:min(size(tableau_W_glider))
%  tab_pitch1=tableau_pitch(:,i);
%  tab_glider1=tableau_W_glider(:,i);
% to_ign=isnan(tab_pitch1);
% to_ign1=isnan(tab_glider1);
% to_ign2=to_ign+to_ign1;
% tab_pitch1(to_ign2 ~= 0)=[];
% 
% pression_p=[0:600/(length(tab_pitch1)-1):600];
% tab_pitch1_1=detrend(tab_pitch1)+mean(tab_pitch1); 
% 
% tab_glider1(to_ign2 ~= 0)=[];
% tab_glider1_1=detrend(tab_glider1)+mean(tab_glider1);
% 
% % tab_pitch1_1=smoothdata(tab_pitch1_1,'movmedian','SmoothingFactor',0.05);
%  %tab_glider1_1=smoothdata(tab_glider1_1,'movmedian','SmoothingFactor',0.05);
% %Regression
% b1 = tab_pitch1_1\tab_glider1_1;
% X = [ones(length(tab_pitch1_1),1) tab_pitch1_1];
% b = X\tab_glider1_1;
% p = polyfit(tab_pitch1_1,tab_glider1_1,2);
% 
% 
% figure(i)
% subplot(1,3,1)
% plot(tab_pitch1_1,-pression_p)
% title('pitch')
% subplot(1,3,2)
% plot(tab_glider1_1,-pression_p,'b')
% title('glider')
% subplot(1,3,3)
% plot(tab_pitch1_1,tab_glider1_1,'+','Color','#4DBEEE')
% hold on 
% ylabel('W_glider')
% xlabel('Pitch')
% y1=b1*tab_pitch1_1;
% y2 = X*b;
% y3 = polyval(p,tab_pitch1_1);
% plot(tab_pitch1_1,y1,'r','LineWidth',2);
% plot(tab_pitch1_1,y2,'--','LineWidth',2)
% plot(tab_pitch1_1,y3,'--','LineWidth',2,'Color','#7E2F8E')
% Rsq1 = 1 - sum((tab_glider1_1 - y1).^2)/sum((tab_glider1_1 - mean(tab_glider1_1)).^2);
% Rsq2 = 1 - sum((tab_glider1_1 - y2).^2)/sum((tab_glider1_1 - mean(tab_glider1_1)).^2);
% Rsq3 = 1 - sum((tab_glider1_1 - y3).^2)/sum((tab_glider1_1 - mean(tab_glider1_1)).^2);
% legend('W\_glider',['R\_lineaire= ' num2str(Rsq1)],['R\_affine= ' num2str(Rsq2)],['R\_polynome= ' num2str(Rsq3)])
% end
% 
% 
%% %% ----- All ----- %%
% ------------------------------------- %
if All==1
     
q=size(tab_all.tab);
Ww_all_d=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
while i < q(1)+1
    if i == 1 
        Ww_1(j)= Ww_all(i);
        pression(j)=tab_all.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if tab_all.tab(i,3)*tab_all.tab(i-1,3) <0 || i == q(1)
            disp(i)
            B=zeros(length(pression),2);
            B(:,1)=pression;
            B(:,2)=Ww_1;
            C=sortrows(B,1);
            to_igno = [];
            for k=1:length(C)-1
                if  C(k,1) == C(k+1,1) % delete doublon
                to_igno = [to_igno k+1];
                end
            end
            C(to_igno,:)=[];
            
           ww_1=interp1(C(:,1),C(:,2),pi);
           Ww_all_d=[Ww_all_d ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           i=i+1;
           if i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_all(i);
           pression(j)=tab_all.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
    
end
ind = find(abs(Ww_all_d) > 0.025);
Ww_all_d(ind)=0;
size_Ww_all_d=size(Ww_all_d);

figure()
pcolor([1:size_Ww_all_d(2)],-pi,Ww_all_d);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_All')
end

%% %% ----- Diving ----- %%
% ------------------------------------- %
if Diving==1
     
q=size(tab_diving.tab);
Ww_dive=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
while i < q(1)+1
    if i == 1 
        Ww_1(j)= Ww_diving(i);
        pression(j)=tab_diving.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if tab_diving.tab(i-1,1)-tab_diving.tab(i,1) > max(tab_diving.tab(:,1))/2 || i == q(1)
            disp(i)
            B=zeros(length(pression),2);
            B(:,1)=pression;
            B(:,2)=Ww_1;
            C=sortrows(B,1);
            to_igno = [];
            for k=1:length(C)-1
                if  C(k,1) == C(k+1,1) % delete doublon
                to_igno = [to_igno k+1];
                end
            end
            C(to_igno,:)=[];
            
           ww_1=interp1(C(:,1),C(:,2),pi);
           Ww_dive=[Ww_dive ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           i=i+1;
           if i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_diving(i);
           pression(j)=tab_diving.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
    
end
ind = find(abs(Ww_dive) > 0.025);
Ww_dive(ind)=0;
size_Ww_dive =size(Ww_dive);

figure()
pcolor([1:size_Ww_dive(2)],-pi,Ww_dive);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_Diving')
end

%% %% ----- Climbing ----- %%
% ------------------------------------- %
if Climbing==1
     
q=size(tab_climbing.tab);
Ww_climb=[];
a=1;
j=1;
pi=[0:1:640];
i=1;
while i < q(1)+1
    if i == 1 
        Ww_1(j)= Ww_climbing(i);
        pression(j)=tab_climbing.tab(i,1);
        j=j+1;
        i=i+1;
    else
    
       if abs(tab_climbing.tab(i-1,1)-tab_climbing.tab(i,1)) > max(tab_climbing.tab(:,1))/2 || i == q(1)
            disp(i)
            B=zeros(length(pression),2);
            B(:,1)=pression;
            B(:,2)=Ww_1;
            C=sortrows(B,1);
            to_igno = [];
            for k=1:length(C)-1
                if  C(k,1) == C(k+1,1) % delete doublon
                to_igno = [to_igno k+1];
                end
            end
            C(to_igno,:)=[];
            
           ww_1=interp1(C(:,1),C(:,2),pi);
           Ww_climb=[Ww_climb ww_1'];
           Ww_1=[];
           pression=[];
           j=1;
           i=i+1;
           if i==q(1)
              i=i+1; 
           end
       else
           Ww_1(j)=Ww_climbing(i);
           pression(j)=tab_climbing.tab(i,1);
           j=j+1;
           i=i+1;
       end
    end
    
end
ind = find(abs(Ww_climb) > 0.025);
Ww_climb(ind)=0;
size_Ww_climb =size(Ww_climb);

figure()
pcolor([1:size_Ww_climb(2)],-pi,Ww_climb);   % utilisation de la fonction pcolor + "shading interp"
shading interp
colormap(blue_red_cmap)
caxis([-0.025 0.025])
H=colorbar;
ylabel(H,'          ','FontSize',12,'Rotation',0);
grid on
ax = gca;
ax.Layer='Top';
xlabel('Numéro de plongée')
ylabel('- Pression (dbar)')
title('Ww\_Climbing')
end
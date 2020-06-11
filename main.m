%% -- MAIN -- %% 

%% ----- Initialisation ----- %%
% ------------------------------------- %
close all;clc;clear;clear global;
%Data path
addpath('Data_Stage');
addpath('Fonction');
addpath('SeaWater');

%%% Choose what dives you want
explorer.first_dive = 10;
explorer.last_dive = 11;
explorer.first=explorer.first_dive;
explorer.last=explorer.last_dive;

%%% Choose what you want to see
yes=1;
no=0;

Display_Web_Map = no;
Display_Map_Figure = no;
Lat_Lon = no;
Depth_oil_pitch = yes;
Colored_Temperature_Salinity = no;
Display_3D_trip = no; Animation_3D = no;
Temperature_Salinity_Profiles = no;

%For vertical velocity 
vertical_velocities = yes; % always yes in order to have vertical velocities
W_glider_W_Model =yes;  three_loops_method = no;  fminsearch_method = yes; fminsearch_method_bydive=yes;
Parameters_evolution = no;
attack_angle = no; fminsearch_method = yes; %both yes to show angle_attack 
water_velocity_descent_bydive = no; hist_desc = no;
water_velocity_ascent_by_dive=no; 
water_velocity_descent_all = no;
water_velocity_ascent_all = no; hist_asc=no; sbplot = no;
sub_asc_desc = no;
water_velocity_all_bydive = no;


%Constant
Const.d2s = 86400;
Const.g=9.81; % gravité
Const.R = 6371;



load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
%load('BioswotSeaExplorer_Nav&CTD.mat');
load('blue_red_cmap.mat')
%% ----- Read Data ----- %%
% ------------------------------------- %

%dives wanted
i_dive= tableau(:,1) >= explorer.first_dive & tableau(:,1) <= explorer.last_dive ;
tableau = tableau(i_dive,:);   

%Read Data
explorer_all = read_EXPLORER(tableau,explorer,Const);
[explorer_bydive,mx]=by_dive2(tableau,explorer,Const);

%% ----- Display Map Web Browser ----- %%
% ------------------------------------- %

if Display_Web_Map == 1
    webmap('Ocean Basemap')
    indices = linspace(1,explorer.size,50);% nombre de points à afficher
    indices=fix(indices);

    %wmmarker(explorer.lat(indices),explorer.lon(indices));
    wmmarker(explorer.lat(1),explorer.lon(1),'Color','green');
    wmmarker(explorer.lat(end),explorer.lon(end),'Color','red');

    wmline(explorer.lat(indices),explorer.lon(indices));
end

%% ----- Display Map Figure + gif -------------------------- %%
% ----------------------------------------------------------- %
if Display_Map_Figure == 1
    
    pause(1.5)
    %h = figure;
    %axis off
    %filename = 'trajectoire.gif';
    figure('Name','Map','NumberTitle','off','Units','normalized','Position',[0,0.55,0.33,0.42]);
    load coastlines
    axesm('ortho','origin',[45 0]);
    axesm('mercator','MapLatLimit',[42 44],'MapLonLimit',[7 10])
    axis off;
    gridm off;
    framem on;
    %mlabel('equator')
    %plabel('fontweight','bold')

    hold on 
    plotm(coastlat,coastlon)
    plotm(explorer.lat(1),explorer.lon(1),'+g')
     for i=2:700:explorer.size-100
         plotm(explorer.lat(i),explorer.lon(i),'.r')
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

    plotm(explorer.lat(end),explorer.lon(end),'+b')
    hold off

end
%% ----- Display latitude and longitude ----- %%
% ----------------------------------------------- %


if Lat_Lon == 1
    
    pause(1.5)
    figure('Name','Latitude and Longitude','NumberTitle','off','Units','normalized','Position',[0.33,0.55,0.33,0.42]);
    subplot(2,1,1)
    plot(explorer.time,explorer.lat,'LineWidth',2)
    title('Evolution de la latitude')
    ylabel('latitude en °')
    datetick('x',2)

    subplot(2,1,2)
    plot(explorer.time,explorer.lon,'LineWidth',2)
    title('Evolution de la longitude')
    ylabel('longitude en °')
    datetick('x',2)

end

%% %% ----- Display depth, oil volume, pitch ----- %%
% ------------------------------------- %

if Depth_oil_pitch == 1
    

    pause(1.5)
    figure('Name','depth oil volume and pitch','NumberTitle','off','Units','normalized','Position',[0.66,0.55,0.33,0.42]);
   % yyaxis left
    %plot(explorer.time,-explorer.depth,'LineWidth',1)
    hold on 
    %plot(explorer.time,explorer.oil,'-+','LineWidth',1)
    ylabel('profondeur(m), huile (ml)')
    yyaxis right
    plot(explorer_all.time,explorer_all.pitch_filter,'-+','LineWidth',1)
    
    hold off
    title('Profondeur, huile, pitch')
    ylabel('pitch (°)')
    datetick('x',0,'keepticks')
    legend('profondeur','volume huile','pitch')

end

%% %% ----- Display temperature and salinity ----- %%
% ------------------------------------- %

% %partition in dive

if Colored_Temperature_Salinity == 1
    pause(1.5)
    
    tab_ti=[];
    tab_si=[];

    for j= explorer.first_dive:explorer.last_dive %data by dive
       explorer = by_dive(tableau,j,explorer,Const);%data by dive

        to_ign=[];
        max_ind = find(explorer.pressure == max(explorer.pressure));
        i=1;
        while explorer.pressure(i) ~= max(explorer.pressure)
          if  not(explorer.pressure(i)<explorer.pressure(i+1)) % delete variations
            to_ign = [to_ign i+1];
          end
          i=i+1;
        end
        to_ign = [to_ign max_ind:length(explorer.pressure)]; % taking only descent 
        explorer.p_sorted = explorer.pressure;
        explorer.p_sorted(to_ign)=[];
        explorer.temp(to_ign)=[];
        explorer.s(to_ign)=[];
        to_igno = [];
        for i=1:length(explorer.p_sorted)-1
         if  explorer.p_sorted(i) == explorer.p_sorted(i+1) % delete doublon
           to_igno = [to_igno i+1];
         end
        end

        explorer.p_sorted(to_igno)=[];
        explorer.temp(to_igno)=[];
        explorer.s(to_igno)=[];
        pi=[0:1:600];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
        ti=interp1(explorer.p_sorted,explorer.temp,pi); 
        si=interp1(explorer.p_sorted,explorer.s,pi);
        tab_ti=[tab_ti ti'];
        tab_si=[tab_si si'];

    end

    figure('Name','Temperature and Salinity','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.42]);
    subplot(2,1,1)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tab_ti);   % utilisation de la fonction pcolor + "shading interp"
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
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tab_si);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'   psu','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Salinité')

end


  
%% %% ----- Display 3D ----- %%
% ------------------------------------- %

if Display_3D_trip == 1
    pause(0.1)
    explorer = read_EXPLORER(tableau,explorer);
    ind = [];
    for j = 1:explorer.size-1
        if explorer.dive(j) ~= explorer.dive(j+1) && mod(explorer.dive(j),2)==0 %collect first points of dives 
           ind  = [ind j+1]; 
        end
    end

    explorer.lat_i=explorer.lat(ind);
    explorer.time_i=explorer.time(ind);
    time_interp = linspace(explorer.time_i(1),explorer.time_i(end),(explorer.last_dive-explorer.first_dive)*150);
    %figure()
    explorer.lat_interp = interp1(explorer.time_i,explorer.lat_i,time_interp,'spline');
    %plot(explorer.time_i,explorer.lat_i,'o',time_interp,explorer.lat_interp,':.');

    explorer.lon_i=explorer.lon(ind);
    %figure()
    explorer.lon_interp = interp1(explorer.time_i,explorer.lon_i,time_interp,'spline');
    %plot(explorer.time_i,explorer.lon_i,'o',time_interp,explorer.lon_interp,':.');
    
    %Filter on explorer.time
    to_ign=[];
    for i=1:length(explorer.time)-1
        if explorer.time(i)> explorer.time(i+1) || explorer.time(i)==explorer.time(i+1)
            to_ign=[to_ign i];
        end
    end
    explorer.time(to_ign)=[];
    explorer.depth(to_ign)=[];
    explorer.temp(to_ign)=[];

    explorer.depth_interp = interp1(explorer.time,explorer.depth,time_interp,'spline');
    explorer.temp_interp = interp1(explorer.time,explorer.temp,time_interp,'pchip');
    %plot(explorer.time,explorer.temp,'o',time_interp,explorer.temp_interp,':.')
    
    X = explorer.lat_interp;
    Y = explorer.lon_interp;
    Z = -explorer.depth_interp;
    C = explorer.temp_interp;
    Z(end)=NaN;
    C(end)=NaN;
    figure('Name','3D_Trip_Temperature','NumberTitle','off','Units','normalized','Position',[0.33,0.05,0.33,0.42]);

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
    explorer = read_EXPLORER(tableau,explorer);
    
    figure('Name','Temperature_Salinity_Density_Profiles','NumberTitle','off','Units','normalized','Position',[0.66,0.05,0.33,0.42]);
    subplot(1,3,1)
    plot(explorer.temp,-explorer.pressure,'b')
    title('Température')
    ylabel('Pression (dbar)')
    xlabel('°C')
    subplot(1,3,2)
    plot(explorer.s,-explorer.pressure,'r')
    title('Salinité')
    xlabel('psu')
    subplot(1,3,3)
    plot(explorer.dens,-explorer.pressure,'k')
    title('Masse volumique')
    xlabel('kg/m^3')

end

%% %% ----- Display vertical velocities ----- %%
% ------------------------------------- %

% W : vertical water velocity
% W_model : velocity from the model
% W_glider : Pressure/time
% W = W_glider - W_model
if vertical_velocities == 1
    pause(1.5)
  


%%%%%%%%%%%%%%%% Méthode 3 boucles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------%

if three_loops_method == 1 
    
    F=[];
    Triplet = [];
    dt=25;
    V0_min=0.0572;
    V0_max=0.0575;
    alpha_min=2e-10;
    alpha_max=5e-10;
    Cd_min=0.08;
    Cd_max=0.12;
    W_model_tot=[];
    A_tot=[];
    for V0=V0_min:(V0_max-V0_min)/dt:V0_max
        for alpha=alpha_min:(alpha_max-alpha_min)/dt:alpha_max 
            for Cd=Cd_min:(Cd_max-Cd_min)/dt:Cd_max
               [W_model] = flight_model(explorer.pressure,explorer.dens,...
                             explorer.pitch,explorer.oil,explorer.temp,V0,alpha,Cd,explorer.M);
                %B=abs(mean((explorer.W_glider'.^2-W_model(1:end-5).^2)));
               A=explorer.W_glider'.^2-W_model(1:end-5).^2;
                B=sum(abs(A));       
                disp(B)
                F=[F B];
                Triplet = [Triplet [V0 alpha Cd]'];
            end 
        end
    end
    % 
    opt =find(F == min(abs(F)) | F == -min(abs(F)));
    disp(min(abs(F)))
    explorer.V0 = Triplet(1,opt);
    explorer.alpha = Triplet(2,opt);
    explorer.Cd=Triplet(3,opt);

    [W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(explorer.pressure,explorer.dens,...
                             explorer.pitch,explorer.oil,explorer.temp,explorer.V0,explorer.alpha,explorer.Cd,explorer.M);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%% Méthode fminsearch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---------------------------------------------------------------------%

 if fminsearch_method_bydive == 1
     W_model_all=[];
     W_glider_all=[];
     time_all=[];
     Tri=[];
    for k=1:explorer.last_dive-explorer.first_dive+1
    
        global pres dens pitch oil_vol Wglider temp mg
        pres=explorer_bydive.pressure(:,k);
        dens=explorer_bydive.dens(:,k);
        pitch=explorer_bydive.pitch_filter(:,k);
        oil_vol=explorer_bydive.oil(:,k);
        Wglider=explorer_bydive.W_glider_filter(:,k);
        temp=explorer_bydive.temp(:,k);
        mg=explorer_all.M;
    %     aw=3.7;
    %     ah=2.4;
    %     Cd1w = 0.78;
    %     Cd1h = 2.1;
    %     alphat = 7.05e-5;
        param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
        options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
        [x,fval,exitflag,output] = fminsearch('cost_bydive',param0,options);
        [W_model] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);
         %W_model = fillmissing(W_model,'next');
         
         %[W_model] = flight_model(pres,dens,pitch,oil_vol,temp,param0(1),param0(2),param0(3),mg);
         Tri = [Tri x'];
         W_model_all= [W_model_all W_model'];
         W_glider_all = [W_glider_all explorer_bydive.W_glider_filter(:,k)'];
         time_all = [time_all explorer_bydive.time(:,k)'];
    
    end
    %[W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);
   
 
%     %Remove spikes
%     ind = find(W_model<-0.3);
%     W_model(ind)=NaN;
%     
if attack_angle == 1 
    pause(1.5)
    figure('Name','Attack angle','NumberTitle','off','Units','normalized','Position',[0,0.55,0.33,0.42]);
    histogram(att_deg,200)
    title('Attack angle')
    xlabel('Attack angle °')
    xlim([-4 5])

end


 end
 if fminsearch_method == 1
  explorer_all=read_EXPLORER(tableau,explorer,Const);
    
        global pres dens pitch oil_vol Wglider temp mg
        pres=explorer_all.pressure;
        dens=explorer_all.dens;
        pitch=explorer_all.pitch_filter';
        oil_vol=explorer_all.oil;
        Wglider=explorer_all.W_glider_filter;
        temp=explorer_all.temp;
        mg=explorer_all.M;
    %     aw=3.7;
    %     ah=2.4;
    %     Cd1w = 0.78;
    %     Cd1h = 2.1;
    %     alphat = 7.05e-5;
        param0 = [explorer_all.V0,explorer_all.alpha,explorer_all.Cd];
        options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
        [x,fval,exitflag,output] = fminsearch('cost',param0,options);
        %[W_model] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);
         %W_model = fillmissing(W_model,'next');
         
         %[W_model] = flight_model(pres,dens,pitch,oil_vol,temp,param0(1),param0(2),param0(3),mg);
        
    
end
    [W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);
   
 
%     %Remove spikes
%     ind = find(W_model<-0.3);
%     W_model(ind)=NaN;
%     
if attack_angle == 1 
    pause(1.5)
    figure('Name','Attack angle','NumberTitle','off','Units','normalized','Position',[0,0.55,0.33,0.42]);
    histogram(att_deg,200)
    title('Attack angle')
    xlabel('Attack angle °')
    xlim([-4 5])

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------%


if W_glider_W_Model == 1 
    pause(1.5)

%W = explorer.W_glider'-W_model(1:end-5);                    
  figure('Name','Wglider and Wmodel','NumberTitle','off','Units','normalized','Position',[0.33,0.55,0.33,0.42]);
  plot(time_all, W_glider_all,'LineWidth',1.5)
  title('Vitesses verticales du glider')
  hold on 
  plot(time_all,W_model_all,'LineWidth',1.5)
  %plot(explorer.time(1:end-5),W,'LineWidth',1.5)
  datetick('x',0,'keepticks')
   legend('W\_glider','W\_model')
  figure()
  plot(explorer_all.time, explorer_all.W_glider_filter,'LineWidth',1.5)
  title('Vitesses verticales du glider')
  hold on 
  plot(explorer_all.time,W_model,'LineWidth',1.5)
  datetick('x',0,'keepticks')
  legend('W\_glider','W\_model')
end 

if Parameters_evolution == 1
    pause(1.5)
 
    Triplet = [];   %%%%%évolution des paramètres
    for p=explorer.first:explorer.last-3    %%fenetre glissante 
        
        load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
        explorer.first = p;
        explorer.last = p+3;
        i_dive= tableau(:,1) >= explorer.first & tableau(:,1) <= explorer.last ;
        tableau = tableau(i_dive,:);   

        %Read Data
        explorer = read_EXPLORER(tableau,explorer,Const);
            
        
    global pres dens pitch oil_vol Wglider temp mg
    pres=explorer.pressure;
    dens=explorer.dens;
    pitch=explorer.pitch_filter';
    oil_vol=explorer.oil;
    Wglider=explorer.W_glider_filter;
    temp=explorer.temp;
    mg=explorer.M;

    param0 = [explorer.V0,explorer.alpha,explorer.Cd];
    options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
    [x,fval,exitflag,output] = fminsearch('cost',param0,options);
    
    Triplet = [Triplet x'];

    end
    figure('Name','Parameters evolution','NumberTitle','off','Units','normalized','Position',[0.66,0.55,0.33,0.42]);
    title('Evolution des paramètres optimisés')
    subplot(3,1,1)
    plot([explorer.first_dive:explorer.last_dive-3],Triplet(3,:),'Color','#0072BD','LineWidth',2)
    xlabel('Numéro de plongée')
    ylabel('Cd')
    subplot(3,1,2)
    plot([explorer.first_dive:explorer.last_dive-3],Triplet(1,:),'Color','#D95319','LineWidth',2)
    xlabel('Numéro de plongée')
    ylabel('V0')
    subplot(3,1,3)
    plot([explorer.first_dive:explorer.last_dive-3],Triplet(2,:),'Color','#77AC30','LineWidth',2)
    xlabel('Numéro de plongée')
    ylabel('eps')

end

end

%% %% ----- Display water vertical velocities ----- %%
% ------------------------------------- %

if water_velocity_descent_bydive == 1 


    load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww1=[];
for k=1:explorer.last_dive-explorer.first_dive+1 
       %explorer = by_dive(tableau,j,explorer,Const);%data by dive
   
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
       %descent
       to_ign=[];
        max_ind = find(explorer2_bydive.pressure == max(explorer2_bydive.pressure));
        i=1;
        while explorer2_bydive.pressure(i) ~= max(explorer2_bydive.pressure)
          if  not(explorer2_bydive.pressure(i)<explorer2_bydive.pressure(i+1)) % delete variations
            to_ign = [to_ign i+1];
          end
          i=i+1;
        end
        to_ign = [to_ign max_ind:length(explorer2_bydive.pressure)]; % taking only descent 
        explorer2_bydive.p_sorted = explorer2_bydive.pressure;
        explorer2_bydive.p_sorted(to_ign)=[];
        explorer2_bydive.temp(to_ign)=[];
        explorer2_bydive.s(to_ign)=[];
        explorer2_bydive.dens(to_ign)=[];
        explorer2_bydive.W_glider_filter(to_ign)=[];
        W_model(to_ign)=[];
        explorer2_bydive.time(to_ign)=[];
    
 

    to_igno = [];
    for i=1:length(explorer2_bydive.p_sorted)-1
     if  explorer2_bydive.p_sorted(i) == explorer2_bydive.p_sorted(i+1) % delete doublon
       to_igno = [to_igno i+1];
     end
    end
    
    explorer2_bydive.p_sorted(to_igno)=[];
    explorer2_bydive.temp(to_igno)=[];
    explorer2_bydive.s(to_igno)=[];
    explorer2_bydive.dens(to_igno)=[];
    explorer2_bydive.W_glider_filter(to_igno)=[];
    W_model(to_igno)=[];
    explorer2_bydive.time(to_igno)=[];
    
    Ww=explorer2_bydive.W_glider_filter-W_model;
  
  
    
    pi=[0:1:600];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
    ti=interp1(explorer2_bydive.p_sorted,explorer2_bydive.temp,pi); 
    si=interp1(explorer2_bydive.p_sorted,explorer2_bydive.s,pi);
    di=interp1(explorer2_bydive.p_sorted,explorer2_bydive.dens,pi);
    wg=interp1(explorer2_bydive.p_sorted,explorer2_bydive.W_glider_filter,pi);
    wm=interp1(explorer2_bydive.p_sorted,W_model,pi);
    ww=interp1(explorer2_bydive.p_sorted,Ww,pi);
    
         %Filter on ww 
    ind = find(abs(ww) > 0.02);
    ww(ind)=NaN;
   
    tableau_ti=[tableau_ti ti'];
    tableau_si=[tableau_si si'];
    tableau_di=[tableau_di di'];
    tableau_w_glider=[tableau_w_glider wg'];
    tableau_w_model=[tableau_w_model wm'];
    tableau_ww1=[tableau_ww1 ww'];
    
%     figure()
%     plot(explorer.time,explorer.W_glider)
%     hold on 
%     plot(explorer.time,W_model)
    
end

figure('Name','Water vertical velocity descent','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.42]);

if sbplot == 1

    subplot(4,1,1)
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

    subplot(4,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_si);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'     psu','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Salinité')

    subplot(4,1,3)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_di);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          kg/m^3','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Densité')

    % subplot(4,1,4)
    % pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_w_glider);   % utilisation de la fonction pcolor + "shading interp"
    % shading interp
    % H=colorbar;
    % ylabel(H,'          ','FontSize',12,'Rotation',0);
    % grid on
    % ax = gca;
    % ax.Layer='Top';
    % xlabel('Numéro de plongée')
    % ylabel('- Pression (dbar)')
    % title('W\_glider')

    % subplot(4,1,3)
    % pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_w_model);   % utilisation de la fonction pcolor + "shading interp"
    % shading interp
    % H=colorbar;
    % ylabel(H,'          ','FontSize',12,'Rotation',0);
    % grid on
    % ax = gca;
    % ax.Layer='Top';
    % xlabel('Numéro de plongée')
    % ylabel('- Pression (dbar)')
    % title('W\_model')
    % 
    subplot(4,1,4)
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
    figure('Name','Attack angle','NumberTitle','off','Units','normalized','Position',[0,0.55,0.33,0.42]);
    histogram(tableau_ww1,100)
    title('velocities')
    xlabel('velocities m/s')
    %xlim([-4 5])

end
end

if water_velocity_descent_all == 1 


    load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww2=[];
for j= explorer.first_dive:explorer.last_dive 
       explorer = by_dive(tableau,j,explorer,Const);%data by dive
   
   global pres dens pitch oil_vol Wglider temp mg
   pres=explorer.pressure;
   dens=explorer.dens;
   pitch=explorer.pitch;
   oil_vol=explorer.oil;
   Wglider=explorer.W_glider;
   temp=explorer.temp;
   mg=explorer_all.M;

   [W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);

      
       %descent
       to_ign=[];
        max_ind = find(explorer.pressure == max(explorer.pressure));
        i=1;
        while explorer.pressure(i) ~= max(explorer.pressure)
          if  not(explorer.pressure(i)<explorer.pressure(i+1)) % delete variations
            to_ign = [to_ign i+1];
          end
          i=i+1;
        end
        to_ign = [to_ign max_ind:length(explorer.pressure)]; % taking only descent 
        explorer.p_sorted = explorer.pressure;
        explorer.p_sorted(to_ign)=[];
        explorer.temp(to_ign)=[];
        explorer.s(to_ign)=[];
        explorer.dens(to_ign)=[];
        explorer.W_glider(to_ign)=[];
        W_model(to_ign)=[];
        explorer.time(to_ign)=[];
    
 

    to_igno = [];
    for i=1:length(explorer.p_sorted)-1
     if  explorer.p_sorted(i) == explorer.p_sorted(i+1) % delete doublon
       to_igno = [to_igno i+1];
     end
    end
    
    explorer.p_sorted(to_igno)=[];
    explorer.temp(to_igno)=[];
    explorer.s(to_igno)=[];
    explorer.dens(to_igno)=[];
    explorer.W_glider(to_igno)=[];
    W_model(to_igno)=[];
    explorer.time(to_igno)=[];
    
    Ww=explorer.W_glider'-W_model;
  
  
    
    pi=[0:1:600];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
    ti=interp1(explorer.p_sorted,explorer.temp,pi); 
    si=interp1(explorer.p_sorted,explorer.s,pi);
    di=interp1(explorer.p_sorted,explorer.dens,pi);
    wg=interp1(explorer.p_sorted,explorer.W_glider',pi);
    wm=interp1(explorer.p_sorted,W_model,pi);
    ww=interp1(explorer.p_sorted,Ww,pi);
    
         %Filter on ww 
    ind = find(abs(ww) > 0.03);
    ww(ind)=NaN;
   
    tableau_ti=[tableau_ti ti'];
    tableau_si=[tableau_si si'];
    tableau_di=[tableau_di di'];
    tableau_w_glider=[tableau_w_glider wg'];
    tableau_w_model=[tableau_w_model wm'];
    tableau_ww2=[tableau_ww2 ww'];
    
%     figure()
%     plot(explorer.time,explorer.W_glider)
%     hold on 
%     plot(explorer.time,W_model)
    
end
figure('Name','Water vertical velocity descent','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.42]);

if sbplot == 1

    subplot(4,1,1)
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

    subplot(4,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_si);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'     psu','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Salinité')

    subplot(4,1,3)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_di);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          kg/m^3','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Densité')

    % subplot(4,1,4)
    % pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_w_glider);   % utilisation de la fonction pcolor + "shading interp"
    % shading interp
    % H=colorbar;
    % ylabel(H,'          ','FontSize',12,'Rotation',0);
    % grid on
    % ax = gca;
    % ax.Layer='Top';
    % xlabel('Numéro de plongée')
    % ylabel('- Pression (dbar)')
    % title('W\_glider')

    % subplot(4,1,3)
    % pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_w_model);   % utilisation de la fonction pcolor + "shading interp"
    % shading interp
    % H=colorbar;
    % ylabel(H,'          ','FontSize',12,'Rotation',0);
    % grid on
    % ax = gca;
    % ax.Layer='Top';
    % xlabel('Numéro de plongée')
    % ylabel('- Pression (dbar)')
    % title('W\_model')
    % 
    subplot(4,1,4)
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
    H=colorbar;
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
    figure('Name','Attack angle','NumberTitle','off','Units','normalized','Position',[0,0.55,0.33,0.42]);
    histogram(tableau_ww1,100)
    title('velocities')
    xlabel('velocities m/s')
    %xlim([-4 5])

end
end
%-----------------------------------------------------------------------%
if water_velocity_ascent_all == 1 

     load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww3=[];
for j= explorer.first_dive:explorer.last_dive 
       explorer = by_dive(tableau,j,explorer,Const);%data by dive
       
   
   global pres dens pitch oil_vol Wglider temp mg
   pres=explorer.pressure;
   dens=explorer.dens;
   pitch=explorer.pitch;
   oil_vol=explorer.oil;
   Wglider=explorer.W_glider;
   temp=explorer.temp;
   mg=explorer_all.M;

   [W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);

       % remontada
        to_ign=[];
        max_ind = find(explorer.pressure == max(explorer.pressure));
        i=max_ind;
        while explorer.pressure(i) ~= explorer.pressure(end)
          if  not(explorer.pressure(i)>explorer.pressure(i+1)) % delete variations
            to_ign = [to_ign i+1];
          end
          i=i+1;
        end

       to_ign = [to_ign 1:max_ind]; % taking only remontada 
       explorer.p_sorted = explorer.pressure;
       explorer.p_sorted(to_ign)=[];
       explorer.temp(to_ign)=[];
       explorer.s(to_ign)=[];
       explorer.dens(to_ign)=[];
       explorer.W_glider(to_ign)=[];
       W_model(to_ign)=[];
       explorer.time(to_ign)=[];


        %remontada
        ind1=[];
        o=length(explorer.p_sorted);
        while explorer.p_sorted(o) < explorer.p_sorted(o-1) && o > 2
          % delete variations2
            ind1 = [ind1 o];
          o=o-1;
        end
        %remontada
        explorer.p_sorted=explorer.p_sorted(ind1);
        explorer.temp=explorer.temp(ind1);
        explorer.s=explorer.s(ind1);
        explorer.dens=explorer.dens(ind1);
        explorer.W_glider=explorer.W_glider(ind1);
        W_model=W_model(ind1);
        explorer.time=explorer.time(ind1);

    to_igno = [];
    for i=1:length(explorer.p_sorted)-1
     if  explorer.p_sorted(i) == explorer.p_sorted(i+1) % delete doublon
       to_igno = [to_igno i+1];
     end
    end
    
    explorer.p_sorted(to_igno)=[];
    explorer.temp(to_igno)=[];
    explorer.s(to_igno)=[];
    explorer.dens(to_igno)=[];
    explorer.W_glider(to_igno)=[];
    W_model(to_igno)=[];
    explorer.time(to_igno)=[];
    
    Ww=explorer.W_glider'-W_model;
  
  
    
    pi=[0:1:600];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
    ti=interp1(explorer.p_sorted,explorer.temp,pi); 
    si=interp1(explorer.p_sorted,explorer.s,pi);
    di=interp1(explorer.p_sorted,explorer.dens,pi);
    wg=interp1(explorer.p_sorted,explorer.W_glider',pi);
    wm=interp1(explorer.p_sorted,W_model,pi);
    ww=interp1(explorer.p_sorted,Ww,pi);
    
         %Filter on ww 
    ind = find(abs(ww) > 0.03);
    ww(ind)=NaN;
   
    tableau_ti=[tableau_ti ti'];
    tableau_si=[tableau_si si'];
    tableau_di=[tableau_di di'];
    tableau_w_glider=[tableau_w_glider wg'];
    tableau_w_model=[tableau_w_model wm'];
    tableau_ww3=[tableau_ww3 ww'];
    
%     figure()
%     plot(explorer.time,explorer.W_glider)
%     hold on 
%     plot(explorer.time,W_model)
    
end

 figure('Name','Water vertical velocity ascent','NumberTitle','off','Units','normalized','Position',[0.33,0.05,0.33,0.42]);

 if sbplot == 1
 
     subplot(4,1,1)
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

    subplot(4,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_si);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'     psu','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Salinité')

    subplot(4,1,3)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_di);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          kg/m^3','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Densité')

    % subplot(4,1,4)
    % pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_w_glider);   % utilisation de la fonction pcolor + "shading interp"
    % shading interp
    % H=colorbar;
    % ylabel(H,'          ','FontSize',12,'Rotation',0);
    % grid on
    % ax = gca;
    % ax.Layer='Top';
    % xlabel('Numéro de plongée')
    % ylabel('- Pression (dbar)')
    % title('W\_glider')

    % subplot(4,1,3)
    % pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_w_model);   % utilisation de la fonction pcolor + "shading interp"
    % shading interp
    % H=colorbar;
    % ylabel(H,'          ','FontSize',12,'Rotation',0);
    % grid on
    % ax = gca;
    % ax.Layer='Top';
    % xlabel('Numéro de plongée')
    % ylabel('- Pression (dbar)')
    % title('W\_model')
    % 
    subplot(4,1,4)
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
    %colormap(blue_red_cmap)
    H=colorbar;
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
    figure('Name','Attack angle','NumberTitle','off','Units','normalized','Position',[0,0.55,0.33,0.42]);
    histogram(tableau_ww2,100)
    title('velocities')
    xlabel('velocities m/s')
    %xlim([-4 5])

end
 
end

if water_velocity_ascent_by_dive == 1 

     load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
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
       
       % remontada
        to_ign=[];
        max_ind = find(explorer2_bydive.pressure == max(explorer2_bydive.pressure));
        i=max_ind;
        while i ~= length(explorer2_bydive.pressure)
          if  not(explorer2_bydive.pressure(i)>explorer2_bydive.pressure(i+1)) % delete variations
            to_ign = [to_ign i+1];
          end
          i=i+1;
        end

       to_ign = [to_ign 1:max_ind]; % taking only remontada 
       explorer2_bydive.p_sorted = explorer2_bydive.pressure;
       explorer2_bydive.p_sorted(to_ign)=[];
       explorer2_bydive.temp(to_ign)=[];
       explorer2_bydive.s(to_ign)=[];
       explorer2_bydive.dens(to_ign)=[];
       explorer2_bydive.W_glider_filter(to_ign)=[];
       W_model(to_ign)=[];
       explorer2_bydive.time(to_ign)=[];


        %remontada
        ind1=[];
        o=length(explorer2_bydive.p_sorted);
        while explorer2_bydive.p_sorted(o) < explorer2_bydive.p_sorted(o-1) && o > 2
          % delete variations2
            ind1 = [ind1 o];
          o=o-1;
        end
        %remontada
        explorer2_bydive.p_sorted=explorer2_bydive.p_sorted(ind1);
        explorer2_bydive.temp=explorer2_bydive.temp(ind1);
        explorer2_bydive.s=explorer2_bydive.s(ind1);
        explorer2_bydive.dens=explorer2_bydive.dens(ind1);
        explorer2_bydive.W_glider_filter=explorer2_bydive.W_glider_filter(ind1);
        W_model=W_model(ind1);
        explorer2_bydive.time=explorer2_bydive.time(ind1);

    to_igno = [];
    for i=1:length(explorer2_bydive.p_sorted)-1
     if  explorer2_bydive.p_sorted(i) == explorer2_bydive.p_sorted(i+1) % delete doublon
       to_igno = [to_igno i+1];
     end
    end
    
    explorer2_bydive.p_sorted(to_igno)=[];
    explorer2_bydive.temp(to_igno)=[];
    explorer2_bydive.s(to_igno)=[];
    explorer2_bydive.dens(to_igno)=[];
    explorer2_bydive.W_glider_filter(to_igno)=[];
    W_model(to_igno)=[];
    explorer2_bydive.time(to_igno)=[];
    
    Ww=explorer2_bydive.W_glider_filter-W_model;
  
  
    
    pi=[0:1:600];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
    ti=interp1(explorer2_bydive.p_sorted,explorer2_bydive.temp,pi); 
    si=interp1(explorer2_bydive.p_sorted,explorer2_bydive.s,pi);
    di=interp1(explorer2_bydive.p_sorted,explorer2_bydive.dens,pi);
    wg=interp1(explorer2_bydive.p_sorted,explorer2_bydive.W_glider_filter,pi);
    wm=interp1(explorer2_bydive.p_sorted,W_model,pi);
    ww=interp1(explorer2_bydive.p_sorted,Ww,pi);
    
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

 figure('Name','Water vertical velocity ascent','NumberTitle','off','Units','normalized','Position',[0.33,0.05,0.33,0.42]);

 if sbplot == 1
 
     subplot(4,1,1)
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

    subplot(4,1,2)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_si);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'     psu','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Salinité')

    subplot(4,1,3)
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_di);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          kg/m^3','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Densité')

    % subplot(4,1,4)
    % pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_w_glider);   % utilisation de la fonction pcolor + "shading interp"
    % shading interp
    % H=colorbar;
    % ylabel(H,'          ','FontSize',12,'Rotation',0);
    % grid on
    % ax = gca;
    % ax.Layer='Top';
    % xlabel('Numéro de plongée')
    % ylabel('- Pression (dbar)')
    % title('W\_glider')

    % subplot(4,1,3)
    % pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_w_model);   % utilisation de la fonction pcolor + "shading interp"
    % shading interp
    % H=colorbar;
    % ylabel(H,'          ','FontSize',12,'Rotation',0);
    % grid on
    % ax = gca;
    % ax.Layer='Top';
    % xlabel('Numéro de plongée')
    % ylabel('- Pression (dbar)')
    % title('W\_model')
    % 
    subplot(4,1,4)
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
    figure('Name','Attack angle','NumberTitle','off','Units','normalized','Position',[0,0.55,0.33,0.42]);
    histogram(tableau_ww2,100)
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

% 
% load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
% pause(1.5)
% tableau_ti=[];
% tableau_si=[];
% tableau_di=[];
% tableau_pi=[];
% tableau_w_glider=[];
% tableau_w_model=[];
% tableau_ww5=[];
% for k=1:explorer.last_dive-explorer.first_dive+1 
%        
%    
%    global pres dens pitch oil_vol Wglider temp mg
%   pres=explorer_bydive.pressure(:,k);
%   dens=explorer_bydive.dens(:,k);
%   pitch=explorer_bydive.pitch_filter(:,k);
%   oil_vol=explorer_bydive.oil(:,k);
%   Wglider=explorer_bydive.W_glider_filter(:,k);
%   temp=explorer_bydive.temp(:,k);
%   mg=explorer_all.M;
% 
%    [W_model] = flight_model(pres,dens,pitch,oil_vol,temp,Tri(1,k),Tri(2,k),Tri(3,k),mg);
%    W_model = fillmissing(W_model,'previous');
%       %explorer2_bydive=explorer_bydive;
%       explorer2_bydive.temp=explorer_bydive.temp(:,k);
%       explorer2_bydive.s=explorer_bydive.s(:,k);
%       explorer2_bydive.dens=explorer_bydive.dens(:,k);
%       explorer2_bydive.W_glider_filter=explorer_bydive.W_glider_filter(:,k);
%        explorer2_bydive.time=explorer_bydive.time(:,k);
%        explorer2_bydive.pressure=explorer_bydive.pressure(:,k);
%        
% %     
% 
% 
% to_igno = [];
%     for i=1:length(explorer2_bydive.time)-1
%      if  explorer2_bydive.time(i) == explorer2_bydive.time(i+1) % delete doublon
%        to_igno = [to_igno i+1];
%      end
%     end
%     
%     explorer2_bydive.time(to_igno)=[];
%     explorer2_bydive.temp(to_igno)=[];
%     explorer2_bydive.pressure(to_igno)=[];
%     explorer2_bydive.s(to_igno)=[];
%     explorer2_bydive.dens(to_igno)=[];
%     explorer2_bydive.W_glider_filter(to_igno)=[];
%     W_model(to_igno)=[];
%  
%     
%     Ww=explorer2_bydive.W_glider_filter-W_model;
%   
%   
%     
%    % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
%    
%     %Remove NaN
%     ind=isnan(explorer2_bydive.time);
%     explorer2_bydive.time(ind)=[];
%     explorer2_bydive.temp(ind)=[];
%     explorer2_bydive.s(ind)=[];
%     explorer2_bydive.pressure(ind)=[];
%     explorer2_bydive.dens(ind)=[];
%     explorer2_bydive.W_glider_filter(ind)=[];
%     W_model(ind)=[];
%     Ww(ind)=[];
%     %explorer2_bydive.temp=fillmissing(explorer2_bydive.temp,'previous');
%     
%      T=[explorer2_bydive.time(1):(explorer2_bydive.time(end)-explorer2_bydive.time(1))/1000:explorer2_bydive.time(end)];  
%     
%     
%     ti=interp1(explorer2_bydive.time,explorer2_bydive.temp,T); 
%     si=interp1(explorer2_bydive.time,explorer2_bydive.s,T);
%     pi=interp1(explorer2_bydive.time,explorer2_bydive.pressure,T);
%     di=interp1(explorer2_bydive.time,explorer2_bydive.dens,T);
%     wg=interp1(explorer2_bydive.time,explorer2_bydive.W_glider_filter,T);
%     wm=interp1(explorer2_bydive.time,W_model,T);
%     ww=interp1(explorer2_bydive.time,Ww,T);
%     
%          %Filter on ww 
%     ind = find(abs(ww) > 0.02);
%     ww(ind)=NaN;
%    
%     tableau_ti=[tableau_ti ti'];
%     tableau_si=[tableau_si si'];
%     tableau_di=[tableau_di di'];
%     tableau_pi=[tableau_pi pi'];
%     tableau_w_glider=[tableau_w_glider wg'];
%     tableau_w_model=[tableau_w_model wm'];
%     tableau_ww5=[tableau_ww5 ww'];
%     
% %     figure()
% %     plot(explorer.time,explorer.W_glider)
% %     hold on 
% %     plot(explorer.time,W_model)
%     
% end
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

figure()
    pcolor([explorer.first_dive:explorer.last_dive*2],-pi,tableau_ww5);   % utilisation de la fonction pcolor + "shading interp"
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



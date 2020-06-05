%% -- MAIN -- %% 

%% ----- Initialisation ----- %%
% ------------------------------------- %
close all;clc;clear;clear global;
%Data path
addpath('Data_Stage');
addpath('Fonction');
addpath('SeaWater');

%%% Choose what dives you want
explorer.first_dive = 1;
explorer.last_dive = 63;

%%% Choose what you want to see
yes=1;
no=0;

Display_Web_Map = no;
Display_Map_Figure = no;
Lat_Lon = no;
Depth_oil_pitch = no;
Colored_Temperature_Salinity = no;
Display_3D_trip = no; Animation_3D = no;
Temperature_Salinity_Profiles = no;

%For vertical velocity 
vertical_velocities = yes; % always yes in order to have vertical velocities
W_glider_W_Model = yes;  three_loops_method = no;  fminsearch_method = yes;
Parameters_evolution = yes;
attack_angle = yes; fminsearch_method = yes; %both yes to show angle_attack 
water_velocity_descent = yes;
water_velocity_ascent = yes; sbplot = no;


%Constant
Const.d2s = 86400;
Const.g=9.81; % gravité
Const.R = 6371;



load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');

%% ----- Read Data ----- %%
% ------------------------------------- %

%dives wanted
i_dive= tableau(:,1) >= explorer.first_dive & tableau(:,1) <= explorer.last_dive ;
tableau = tableau(i_dive,:);   

%Read Data
explorer = read_EXPLORER(tableau,explorer);


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
    figure('Name','Map','NumberTitle','off','Units','centimeters','Position',[1,10,11,8]);
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
    figure('Name','Latitude and Longitude','NumberTitle','off','Units','centimeters','Position',[13,10,11,8]);
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
    
    
%     TF = isoutlier(explorer.pitch);
%     ind = find(TF)
%     Aoutlier = explorer.pitch(ind);
%     Afill = filloutliers(explorer.pitch,'next');
    pause(1.5)
    figure('Name','depth oil volume and pitch','NumberTitle','off','Units','centimeters','Position',[25,10,11,8]);
    yyaxis left
    plot(explorer.time,-explorer.depth,'LineWidth',1)
    hold on 
    plot(explorer.time,explorer.oil,'LineWidth',1)
    ylabel('profondeur(m), huile (ml)')
    yyaxis right
    plot(explorer.time,explorer.pitch,'LineWidth',1)
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
    
    tableau_ti=[];
    tableau_si=[];

    for j= explorer.first_dive:explorer.last_dive %data by dive
       explorer = by_dive(tableau,j,explorer);%data by dive

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
        tableau_ti=[tableau_ti ti'];
        tableau_si=[tableau_si si'];

    end

    figure('Name','Temperature and Salinity','NumberTitle','off','Units','centimeters','Position',[1,1,11,7]);
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
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_si);   % utilisation de la fonction pcolor + "shading interp"
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
    figure('Name','3D_Trip_Temperature','NumberTitle','off','Units','centimeters','Position',[13,1,11,7]);

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
    
    figure('Name','Temperature_Salinity_Density_Profiles','NumberTitle','off','Units','centimeters','Position',[25,1,11,7]);
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
    explorer = read_EXPLORER(tableau,explorer);

W_glider = zeros(1,explorer.size-5);
    for k=1:explorer.size-5
    W_glider(k) = (explorer.pressure(k)-explorer.pressure(k+5))/(explorer.time(k+5)-explorer.time(k));
    end

   explorer.W_glider = W_glider./Const.d2s;
 
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

if fminsearch_method == 1

    global pres dens pitch oil_vol Wglider temp mg
    pres=explorer.pressure;
    dens=explorer.dens;
    pitch=explorer.pitch;
    oil_vol=explorer.oil;
    Wglider=explorer.W_glider;
    temp=explorer.temp;
    mg=explorer.M;

    param0 = [explorer.V0,explorer.alpha,explorer.Cd];
    options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
    [x,fval,exitflag,output] = fminsearch('cost',param0,options);
    [W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);


    %Remove spikes
    ind = find(W_model<-0.3);
    W_model(ind)=NaN;
    
if attack_angle == 1 
    pause(1.5)
    figure('Name','Attack angle','NumberTitle','off','Units','centimeters','Position',[1,10,11,7]);
    histogram(att_deg,200)
    title('Attack angle')
    xlabel('Attack angle °')
    xlim([-4 5])

end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------%


if W_glider_W_Model == 1 
    pause(1.5)

%W = explorer.W_glider'-W_model(1:end-5);                    
  figure('Name','Wglider and Wmodel','NumberTitle','off','Units','centimeters','Position',[13,10,11,7]);
  plot(explorer.time(1:end-5), explorer.W_glider,'LineWidth',1.5)
  title('Vitesses verticales du glider')
  hold on 
  plot(explorer.time,W_model,'LineWidth',1.5)
  %plot(explorer.time(1:end-5),W,'LineWidth',1.5)
  datetick('x',0,'keepticks')
  legend('W\_glider','W\_model')
  
end 

if Parameters_evolution == 1
    pause(1.5)
  explorer.first=explorer.first_dive;
  explorer.last=explorer.last_dive;
    Triplet = [];   %%%%%évolution des paramètres
    for p=explorer.first:explorer.last-3    %%fenetre glissante 
        
        load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
        explorer.first = p;
        explorer.last = p+3;
        i_dive= tableau(:,1) >= explorer.first & tableau(:,1) <= explorer.last ;
        tableau = tableau(i_dive,:);   

        %Read Data
        explorer = read_EXPLORER(tableau,explorer);
            
        W_glider = zeros(1,explorer.size-5);
        for k=1:explorer.size-5
            W_glider(k) = (explorer.pressure(k)-explorer.pressure(k+5))/(explorer.time(k+5)-explorer.time(k));
        end
        explorer.W_glider = W_glider./Const.d2s;
        
    global pres dens pitch oil_vol Wglider temp mg
    pres=explorer.pressure;
    dens=explorer.dens;
    pitch=explorer.pitch;
    oil_vol=explorer.oil;
    Wglider=explorer.W_glider;
    temp=explorer.temp;
    mg=explorer.M;

    param0 = [explorer.V0,explorer.alpha,explorer.Cd];
    options = optimset('Display','iter','MaxFunEvals',8000,'MaxIter',8000);
    [x,fval,exitflag,output] = fminsearch('cost',param0,options);
    
    Triplet = [Triplet x'];

    end
    figure('Name','Parameters evolution','NumberTitle','off','Units','centimeters','Position',[25,10,11,7]);
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

if water_velocity_descent == 1 
    
    load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww=[];
for j= explorer.first_dive:explorer.last_dive 
       explorer = by_dive(tableau,j,explorer);%data by dive
       W_glider = zeros(1,explorer.size-5);
    for k=1:explorer.size-5
    W_glider(k) = (explorer.pressure(k)-explorer.pressure(k+5))/(explorer.time(k+5)-explorer.time(k));
    end
    
   explorer.W_glider = W_glider./Const.d2s;
   
   global pres dens pitch oil_vol Wglider temp mg
   pres=explorer.pressure;
   dens=explorer.dens;
   pitch=explorer.pitch;
   oil_vol=explorer.oil;
   Wglider=explorer.W_glider;
   temp=explorer.temp;
   mg=explorer.M;

   [W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);

      %descent
       explorer.pressure=explorer.pressure(1:end-5);
       explorer.temp=explorer.temp(1:end-5);
       explorer.s=explorer.s(1:end-5);
       explorer.dens=explorer.dens(1:end-5);
       W_model=W_model(1:end-5);
       explorer.time=explorer.time(1:end-5);
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
    tableau_ww=[tableau_ww ww'];
    
%     figure()
%     plot(explorer.time,explorer.W_glider)
%     hold on 
%     plot(explorer.time,W_model)
    
end

figure('Name','Water vertical velocity descent','NumberTitle','off','Units','centimeters','Position',[1,1,11,7]);

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
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww);   % utilisation de la fonction pcolor + "shading interp"
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
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww')
    
end

end


if water_velocity_ascent == 1 
     load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');
    pause(1.5)
tableau_ti=[];
tableau_si=[];
tableau_di=[];
tableau_w_glider=[];
tableau_w_model=[];
tableau_ww=[];
for j= explorer.first_dive:explorer.last_dive 
       explorer = by_dive(tableau,j,explorer);%data by dive
       
       W_glider = zeros(1,explorer.size-5);
    for k=1:explorer.size-5
    W_glider(k) = (explorer.pressure(k)-explorer.pressure(k+5))/(explorer.time(k+5)-explorer.time(k));
    end
    
   explorer.W_glider = W_glider./Const.d2s;
   
   global pres dens pitch oil_vol Wglider temp mg
   pres=explorer.pressure;
   dens=explorer.dens;
   pitch=explorer.pitch;
   oil_vol=explorer.oil;
   Wglider=explorer.W_glider;
   temp=explorer.temp;
   mg=explorer.M;

   [W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(pres,dens,pitch,oil_vol,temp,x(1),x(2),x(3),mg);

   
      % remontada
       explorer.pressure=explorer.pressure(6:end);
       explorer.temp=explorer.temp(6:end);
       explorer.s=explorer.s(6:end);
       explorer.dens=explorer.dens(6:end);
       W_model=W_model(6:end);
       explorer.time=explorer.time(6:end);



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
    tableau_ww=[tableau_ww ww'];
    
%     figure()
%     plot(explorer.time,explorer.W_glider)
%     hold on 
%     plot(explorer.time,W_model)
    
end

 figure('Name','Water vertical velocity ascent','NumberTitle','off','Units','centimeters','Position',[13,1,11,7]);

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
    pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww);   % utilisation de la fonction pcolor + "shading interp"
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
     pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww')
 end
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



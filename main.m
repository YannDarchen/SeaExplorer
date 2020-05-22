%% -- MAIN -- %% 

%% ----- Initialisation ----- %%
% ------------------------------------- %

%Clean Data
close all;clc;clear;clear global;

%Data path
addpath('..\Data_Stage');
load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');

%Constant
Const.d2s = 86400;
Const.g=9.81; % gravité

%% ----- Read Data ----- %%
% ------------------------------------- %

%Choose dives wanted
first = 11;
last = 22;
i_dive= tableau(:,1) >=first & tableau(:,1) <= last;
tableau = tableau(i_dive,:);   

%Read Data
explorer.lat = tableau(:,10)/100;
explorer.lon = tableau(:,11)/100;
explorer.time = tableau(:,8);
explorer.depth = tableau(:,9);
explorer.oil = tableau(:,13);
explorer.pitch = tableau(:,15);
explorer.pressure = tableau(:,17);
explorer.temp = tableau(:,18);



%Remove NaN
to_ignore_lat=isnan(explorer.lat);
to_ignore_lon=isnan(explorer.lon);
to_ignore_time=isnan(explorer.time);
to_ignore_depth=isnan(explorer.depth);
to_ignore_oil=isnan(explorer.oil);
to_ignore_pitch=isnan(explorer.pitch);
to_ignore_pressure=isnan(explorer.pressure);
to_ignore_temp=isnan(explorer.temp);

to_ignore=to_ignore_lat+to_ignore_lon+to_ignore_time+to_ignore_depth+to_ignore_oil+to_ignore_pitch+to_ignore_pressure+to_ignore_temp; 


explorer.lat(to_ignore ~= 0 )=[];
explorer.lon(to_ignore ~= 0 )=[];
explorer.time(to_ignore ~= 0 )=[];
explorer.depth(to_ignore ~= 0 )=[];
explorer.oil(to_ignore ~= 0 )=[];
explorer.pitch(to_ignore ~= 0 )=[];
explorer.pressure(to_ignore ~= 0 )=[];
explorer.temp(to_ignore ~= 0 )=[];


dens = sw_smow(explorer.temp);

s=length(explorer.lat);

explorer.M = 59; % Mass in kg
explorer.V0 = 0.051358; %Volume
explorer.alpha = 3*10^-6; %compressibility
explorer.beta =  1.09*10^-4; %coefficient d'expansion thermique 
explorer.T0=20; % température
explorer.a=3; % angle d'attaque degrés 
explorer.S= 0.032365; % surface du glider 
explorer.Cd=0.4; % coefficient de trainée 

   

%% ----- Display Map Web Browser ----- %%
% ------------------------------------- %
% 
% webmap('Ocean Basemap')
% latitude = 43.44; 
% longitude = 8.2; 
% zoom = 9; 
% indices = linspace(1,s,1000);% nombre de points à afficher
% indices=fix(indices);
% 
% wmmarker(explorer.lat(indices),explorer.lon(indices));
% wmmarker(explorer.lat(1),explorer.lon(1),'Color','green');
% wmmarker(explorer.lat(end),explorer.lon(end),'Color','black');
% wmcenter(latitude,longitude,zoom) 

%% ----- Display Map Figure + gif ----- %%
% ------------------------------------- %
% h = figure;
% axis off
% filename = 'trajectoire.gif';
% load coastlines
% axesm('ortho','origin',[45 0]);
% axesm('mercator','MapLatLimit',[42 44],'MapLonLimit',[7 10])
% axis off;
% gridm off;
% framem on;
% %mlabel('equator')
% %plabel('fontweight','bold')
% 
% hold on 
% plotm(coastlat,coastlon)
% plotm(explorer.lat(1),explorer.lon(1),'+g')
%  for i=2:700:s-100
%      plotm(explorer.lat(i),explorer.lon(i),'.r')
%      %title(['Trajectoire du glider jusquau ' num2str(position(i,3)) '/' num2str(position(i,4)) '/' num2str(position(i,5)) '  ' num2str(position(i,6)) ':' num2str(position(i,7))])
%      drawnow
%      % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if i == 2 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05); 
%       end 
% end
% 
% plotm(explorer.lat(end),explorer.lon(end),'+b')
% hold off

%% %% ----- Display latitude and longitude ----- %%
% ------------------------------------- %
% 
% figure()
% subplot(2,1,1)
% plot(explorer.time,explorer.lat)
% title('Evolution de la latitude')
% ylabel('latitude en °')
% datetick('x',2)
% 
% subplot(2,1,2)
% plot(explorer.time,explorer.lon)
% title('Evolution de la longitude')
% ylabel('longitude en °')
% datetick('x',2)

%% %% ----- Display depth, oil volume, pitch ----- %%
% ------------------------------------- %
% 
% figure()
% yyaxis left
% plot(explorer.time,-explorer.depth,'LineWidth',1)
% hold on 
% plot(explorer.time,explorer.oil,'LineWidth',1)
% ylabel('profondeur(m), huile (ml)')
% yyaxis right
% plot(explorer.time,explorer.pitch,'LineWidth',1)
% hold off
% title('Profondeur, huile, pitch')
% ylabel('pitch (°)')
% datetick('x',0,'keepticks')
% legend('profondeur','volume huile','pitch')
% 


%% %% ----- Display temperature ----- %%
% ------------------------------------- %

% %partition in dive
% tableau_ti=[];
% 
% for j= first:last 
%     i_dive1 = tableau(:,1) == j;
%     tableau_dive = tableau(i_dive1,:);
%     explorer.p = tableau_dive(:,17);
%     explorer.t = tableau_dive(:,18);
%     to_ignore_t=isnan(explorer.t);
%     explorer.t(to_ignore_t ~= 0 )=[];
%     explorer.p(to_ignore_t ~= 0 )=[];
%     
%     to_ign=[];
%     max_ind = find(explorer.p == max(explorer.p));
%     i=1;
%     while explorer.p(i) ~= max(explorer.p)
%       if  not(explorer.p(i)<explorer.p(i+1)) % taking only descent 
%         to_ign = [to_ign i+1];
%       end
%       i=i+1;
%     end
%     to_ign = [to_ign max_ind:length(explorer.p)];
%     explorer.p_sorted = explorer.p;
%     explorer.p_sorted(to_ign)=[];
%     explorer.t(to_ign)=[];
%     to_igno = [];
%     for i=1:length(explorer.p_sorted)-1
%      if  explorer.p_sorted(i) == explorer.p_sorted(i+1) % taking only descent 
%        to_igno = [to_igno i+1];
%      end
%     end
%     
%     explorer.p_sorted(to_igno)=[];
%     explorer.t(to_igno)=[];
%     pi=[0:1:600];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
%     ti=interp1(explorer.p_sorted,explorer.t,pi); 
%     tableau_ti=[tableau_ti ti'];
%     
% end
% 
% figure()
% pcolor([first:last],-pi,tableau_ti);   % utilisation de la fonction pcolor + "shading interp"
% shading interp
% H=colorbar;
% ylabel(H,'   (°C)','FontSize',12,'Rotation',0);
% grid on
% ax = gca;
% ax.Layer='Top'
% xlabel('Numéro de plongée')
% ylabel('- Pression (dbar)')
% title('Température')
% 
% 

%% %% ----- Display vertical velocities ----- %%
% ------------------------------------- %

% W : vertical water velocity
% W_model : velocity from the model
% W_glider : Pressure/time
% W = W_glider - W_model

W_glider = zeros(1,s-1);
    for k=1:s-1
    W_glider(k) = (explorer.depth(k)-explorer.depth(k+1))/(explorer.time(k+1)-explorer.time(k));
    end
    W_glider = W_glider./Const.d2s;
 
    
 explorer.V = 5*10^-4 +  explorer.V0*(1-explorer.alpha*explorer.pressure+explorer.beta*(explorer.temp-explorer.T0));% volume
 W_model=zeros(s,1);
 for b=1:s
    W_model(b) = -sqrt((2*(explorer.M-dens(b)*explorer.V(b))*Const.g*sind(explorer.pitch(b)+explorer.a)^3)/(dens(b)*explorer.S*explorer.Cd));
 end
 
figure()
plot(explorer.time(1:end-1), W_glider)
hold on 
plot(explorer.time,W_model)
datetick('x',0,'keepticks')
legend('W\_glider','W\_model')
hold off


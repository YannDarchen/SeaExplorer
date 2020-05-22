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
Const.R = 6371000;

%% ----- Read Data ----- %%
% ------------------------------------- %

%Choose dives wanted
explorer.first_dive = 10;
explorer.last_dive = 20;
i_dive= tableau(:,1) >= explorer.first_dive & tableau(:,1) <= explorer.last_dive;
tableau = tableau(i_dive,:);   

%Read Data
explorer = read_EXPLORER(tableau);


%% ----- Display Map Web Browser ----- %%
% ------------------------------------- %

% webmap('Ocean Basemap')
% indices = linspace(1,explorer.size,50);% nombre de points à afficher
% indices=fix(indices);
% 
% %wmmarker(explorer.lat(indices),explorer.lon(indices));
% wmmarker(explorer.lat(1),explorer.lon(1),'Color','green');
% wmmarker(explorer.lat(end),explorer.lon(end),'Color','red');
% 
% wmline(explorer.lat(indices),explorer.lon(indices));


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
figure()
subplot(2,1,1)
plot(explorer.time,explorer.lat,'+')
title('Evolution de la latitude')
ylabel('latitude en °')
datetick('x',0)

subplot(2,1,2)
plot(explorer.time,explorer.lon)
title('Evolution de la longitude')
ylabel('longitude en °')
datetick('x',2)

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

% W_glider = zeros(1,s-1);
%     for k=1:s-1
%     W_glider(k) = (explorer.depth(k)-explorer.depth(k+1))/(explorer.time(k+1)-explorer.time(k));
%     end
%     W_glider = W_glider./Const.d2s;
%  
%     
%  explorer.V = 5*10^-4 +  explorer.V0*(1-explorer.alpha*explorer.pressure+explorer.beta*(explorer.temp-explorer.T0));% volume
%  W_model=zeros(s,1);
%  for b=1:s
%     W_model(b) = -sqrt((2*(explorer.M-dens(b)*explorer.V(b))*Const.g*sind(explorer.pitch(b)+explorer.a)^3)/(dens(b)*explorer.S*explorer.Cd));
%  end
%  
% figure()
% plot(explorer.time(1:end-1), W_glider)
% hold on 
% plot(explorer.time,W_model)
% datetick('x',0,'keepticks')
% legend('W\_glider','W\_model')
% hold off

%% %% ----- Display 3D ----- %%
% ------------------------------------- %


ind = [];
for j = 1:explorer.size-1
    if explorer.dive(j) ~= explorer.dive(j+1) && mod(explorer.dive(j),2)==0 
       ind  = [ind j+1]; 
    end
end

explorer.lat_i=explorer.lat(ind);
explorer.time_i=explorer.time(ind);
time_interp = linspace(explorer.time_i(1),explorer.time_i(end),2000);
figure()
explorer.lat_interp = interp1(explorer.time_i,explorer.lat_i,time_interp,'spline');
plot(explorer.time_i,explorer.lat_i,'o',time_interp,explorer.lat_interp,':.');

explorer.lon_i=explorer.lon(ind);
figure()
explorer.lon_interp = interp1(explorer.time_i,explorer.lon_i,time_interp,'spline');
plot(explorer.time_i,explorer.lon_i,'o',time_interp,explorer.lon_interp,':.');


explorer.depth_interp = interp1(explorer.time,explorer.depth,time_interp,'spline');


figure()
plot3(explorer.lat_interp,explorer.lon_interp,explorer.depth_interp,'w')
hold on 
for j= 1:length(explorer.lat_interp)
    disp(j)
    plot3(explorer.lat_interp(j),explorer.lon_interp(j),explorer.depth_interp(j),'+k')
    drawnow
end
hold off






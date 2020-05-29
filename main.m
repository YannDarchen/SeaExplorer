%% -- MAIN -- %% 

%% ----- Initialisation ----- %%
% ------------------------------------- %

%Clean Data
close all;clc;clear;clear global;

%Data path
addpath('Data_Stage');
%addpath('Fonction');
addpath('SeaWater');
load('Fumseck_SeaExplorer_Nav&CTD_PROVISOIRE.mat');

%Constant
Const.d2s = 86400;
Const.g=9.81; % gravité
Const.R = 6371;

%% ----- Read Data ----- %%
% ------------------------------------- %

%Choose dives wanted
explorer.first_dive = 10;
explorer.last_dive = 15;
i_dive= tableau(:,1) >= explorer.first_dive & tableau(:,1) <= explorer.last_dive;
tableau = tableau(i_dive,:);   

%Read Data
explorer = read_EXPLORER(tableau,explorer);


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
% figure()
% subplot(2,1,1)
% plot(explorer.time,explorer.lat,'+')
% title('Evolution de la latitude')
% ylabel('latitude en °')
% datetick('x',0)
% 
% subplot(2,1,2)
% plot(explorer.time,explorer.lon)
% title('Evolution de la longitude')
% ylabel('longitude en °')
% datetick('x',2)

%% %% ----- Display depth, oil volume, pitch ----- %%
% ------------------------------------- %

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


%% %% ----- Display temperature and salinity ----- %%
% ------------------------------------- %

% %partition in dive
% tableau_ti=[];
% tableau_si=[];
% 
% for j= explorer.first_dive:explorer.last_dive %data by dive
%     i_dive1 = tableau(:,1) == j;
%     tableau_dive = tableau(i_dive1,:);
%     explorer.p = tableau_dive(:,17); %pression
%     explorer.t = tableau_dive(:,18); %température
%     explorer.c = tableau_dive(:,19); % conductivité
%     
%     %Remove NaN
%     to_ignore_t=isnan(explorer.t);
%     to_ignore_c=isnan(explorer.c);
%     to_ignore = to_ignore_t + to_ignore_c;
%     explorer.t(to_ignore ~= 0 )=[];
%     explorer.c(to_ignore ~= 0 )=[];
%     explorer.p(to_ignore ~= 0 )=[];
%     
%     %Conversion conductivité -> salinité 
%     %Coefficients 
%     a0 = 0.008;
%     a1 = -0.1692;
%     a2 = 25.3851;
%     a3 = 14.0941;
%     a4 = -7.0261;
%     a5 = 2.7081;
%     b0 = 0.0005;
%     b1 = -0.0056;
%     b2 = -0.0066;
%     b3 = -0.0375;
%     b4 = 0.0636;
%     b5 = -0.0144;
%     c0 = 0.6766097;
%     c1 = 0.0200564;
%     c2 = 0.0001104259;
%     c3 = -6.9698e-07;
%     c4 = 0.000000001;
%     k=0.0162;
%     
%  
%     D = ((explorer.c*10)/42.914)./(c0 + c1*explorer.t + c2*(explorer.t.^2)+c3*(explorer.t.^3)+c4*(explorer.t.^4));
%     explorer.s = a0 + a1*D.^0.5 + a2*D + a3*D.^1.5 + a4*D.^2 + a5*D.^2.5 +((explorer.t-15)./(1+k*(explorer.t-15))).*(b0+b1*D.^0.5+b2*D+b3*D.^1.5+b4*D.^2+b5*D.^2.5);
%     
%     to_ign=[];
%     max_ind = find(explorer.p == max(explorer.p));
%     i=1;
%     while explorer.p(i) ~= max(explorer.p)
%       if  not(explorer.p(i)<explorer.p(i+1)) % delete variations
%         to_ign = [to_ign i+1];
%       end
%       i=i+1;
%     end
%     to_ign = [to_ign max_ind:length(explorer.p)]; % taking only descent 
%     explorer.p_sorted = explorer.p;
%     explorer.p_sorted(to_ign)=[];
%     explorer.t(to_ign)=[];
%     explorer.s(to_ign)=[];
%     to_igno = [];
%     for i=1:length(explorer.p_sorted)-1
%      if  explorer.p_sorted(i) == explorer.p_sorted(i+1) % delete doublon
%        to_igno = [to_igno i+1];
%      end
%     end
%     
%     explorer.p_sorted(to_igno)=[];
%     explorer.t(to_igno)=[];
%     explorer.s(to_igno)=[];
%     pi=[0:1:600];  % construction du vecteur pression régulier pi avec un pas de 0.5 dbar qui va servir de base à l'interpolation
%     ti=interp1(explorer.p_sorted,explorer.t,pi); 
%     si=interp1(explorer.p_sorted,explorer.s,pi);
%     tableau_ti=[tableau_ti ti'];
%     tableau_si=[tableau_si si'];
%     
% end
% 
% figure()
% subplot(2,1,1)
% pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ti);   % utilisation de la fonction pcolor + "shading interp"
% shading interp
% H=colorbar;
% ylabel(H,'   (°C)','FontSize',12,'Rotation',0);
% grid on
% ax = gca;
% ax.Layer='Top';
% xlabel('Numéro de plongée')
% ylabel('- Pression (dbar)')
% title('Température')
% 
% subplot(2,1,2)
% pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_si);   % utilisation de la fonction pcolor + "shading interp"
% shading interp
% H=colorbar;
% ylabel(H,'   (°C)','FontSize',12,'Rotation',0);
% grid on
% ax = gca;
% ax.Layer='Top';
% xlabel('Numéro de plongée')
% ylabel('- Pression (dbar)')
% title('Salinité')


%% %% ----- Display vertical velocities ----- %%
% ------------------------------------- %

% W : vertical water velocity
% W_model : velocity from the model
% W_glider : Pressure/time
% W = W_glider - W_model
% 
% W_glider = zeros(1,explorer.size-5);
%     for k=1:explorer.size-5
%     W_glider(k) = (explorer.pressure(k)-explorer.pressure(k+5))/(explorer.time(k+5)-explorer.time(k));
%     end
%    % W_glider1 =smoothdata(W_glider./Const.d2s,'movmedian');
%    W_glider = W_glider./Const.d2s;
% %     SPRAY MODEL
% %  explorer.V = 5*10^-4 +  explorer.V0*(1-explorer.alpha*explorer.pressure+explorer.beta*(explorer.temp-explorer.T0));% volume
% %  W_model=zeros(explorer.size,1);
% %  for b=1:explorer.size
% %     W_model(b) = -sqrt((2*(explorer.M-explorer.dens(b)*explorer.V(b))*Const.g*sind(explorer.pitch(b)+explorer.a)^3)/(explorer.dens(b)*explorer.S*explorer.Cd));
% %  end
% %  
% F=[];
% Triplet = [];
% dt=30;
% V0_min=0.0571;
% V0_max=0.0575;
% alpha_min=2e-10;
% alpha_max=5e-10;
% Cd_min=0.08;
% Cd_max=0.12;
% W_model_tot=[];
% A_tot=[];
% for V0=V0_min:(V0_max-V0_min)/dt:V0_max
%     for alpha=alpha_min:(alpha_max-alpha_min)/dt:alpha_max 
%         for Cd=Cd_min:(Cd_max-Cd_min)/dt:Cd_max
%            [W_model] = flight_model(explorer.pressure,explorer.dens,...
%                          explorer.pitch,explorer.oil,explorer.temp,V0,alpha,Cd,explorer.M);
%             A=(W_glider'.^2-W_model(1:end-5).^2);
%             B=sum(abs(A));       
%             disp(B)
%             F=[F B];
%             Triplet = [Triplet [V0 alpha Cd]'];
%         end 
%     end
% end
% % 
% opt =find(F == min(abs(F)) | F == -min(abs(F)));
% disp(min(abs(F)))
% explorer.V0 = Triplet(1,opt);
% explorer.alpha = Triplet(2,opt);
% explorer.Cd=Triplet(3,opt);
% 
% [W_model,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(explorer.pressure,explorer.dens,...
%                          explorer.pitch,explorer.oil,explorer.temp,explorer.V0,explorer.alpha,explorer.Cd,explorer.M);
% 
% 
% 
% %  %W = W_glider'-W_model(1:end-5);                    
%  figure()
%  plot(explorer.time(1:end-5), W_glider,'LineWidth',1.5)
%  hold on 
% % 
%  plot(explorer.time,W_model,'LineWidth',1.5)
% % %plot(explorer.time(1:end-10),W,'LineWidth',1.5)
%  datetick('x',0,'keepticks')
%  legend('W\_glider','W\_model','W')
% hold off

%% %% ----- Display 3D ----- %%
% ------------------------------------- %

% 
% ind = [];
% for j = 1:explorer.size-1
%     if explorer.dive(j) ~= explorer.dive(j+1) && mod(explorer.dive(j),2)==0 
%        ind  = [ind j+1]; 
%     end
% end
% 
% explorer.lat_i=explorer.lat(ind);
% explorer.time_i=explorer.time(ind);
% time_interp = linspace(explorer.time_i(1),explorer.time_i(end),500);
% figure()
% explorer.lat_interp = interp1(explorer.time_i,explorer.lat_i,time_interp,'spline');
% plot(explorer.time_i,explorer.lat_i,'o',time_interp,explorer.lat_interp,':.');
% 
% explorer.lon_i=explorer.lon(ind);
% figure()
% explorer.lon_interp = interp1(explorer.time_i,explorer.lon_i,time_interp,'spline');
% plot(explorer.time_i,explorer.lon_i,'o',time_interp,explorer.lon_interp,':.');
% 
% 
% explorer.depth_interp = interp1(explorer.time,explorer.depth,time_interp,'spline');
% explorer.temp_interp = interp1(explorer.time,explorer.temp,time_interp,'spline');
% 
% X = explorer.lat_interp;
% Y = explorer.lon_interp;
% Z = -explorer.depth_interp;
% C = explorer.temp_interp;
% Z(end)=NaN;
% C(end)=NaN;
% figure()
% 
% fill3(X,Y,Z,C,'EdgeColor','w','LineWidth',1)
% c = colorbar;
% c.Label.String = 'Température °C';
% title('Trajectoire du glider')
% xlabel('latitude')
% ylabel('longitude')
% zlabel('profondeur (m)')
% hold on 
% for j= 1:length(explorer.lat_interp)
% 
%     fill3(X(j),Y(j),Z(j),C(j),'EdgeColor','interp','Marker','o','MarkerFaceColor','flat')
%     drawnow
% end
% hold off

%% %% ----- Display temperature and salinity profile ----- %%
% ------------------------------------- %

% figure()
% subplot(1,3,1)
% plot(explorer.temp,-explorer.pressure,'b')
% title('Température')
% ylabel('Pression (dbar)')
% xlabel('°C')
% subplot(1,3,2)
% plot(explorer.s,-explorer.pressure,'r')
% title('Salinité')
% subplot(1,3,3)
% plot(explorer.dens,-explorer.pressure,'k')
% title('Masse volumique')
% xlabel('kg/m^3')
% 


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

tableau_pitch=[];
tableau_W_glider=[];
b1=[];
for j= explorer.first_dive:explorer.last_dive 
       explorer = by_dive(tableau,j);%data by dive
    
     pi=[0:1:600];
     ind = find(explorer.pressure == max(explorer.pressure));
     explorer.pressure=explorer.pressure(50:ind-10);
     explorer.pitch=explorer.pitch(50:ind-10);
     explorer.time=explorer.time(50:ind-10);

    W_glider = zeros(1,length(explorer.pressure)-10);
    for k=1:length(explorer.pressure)-10
        W_glider(k) = (explorer.pressure(k)-explorer.pressure(k+10))/(explorer.time(k+10)-explorer.time(k));
    end
    W_glider = W_glider./Const.d2s;
    
    
     pitch=interp1(explorer.pressure,explorer.pitch,pi);
     W_glider = interp1(explorer.pressure(1:end-10),W_glider,pi);
     tableau_pitch=[tableau_pitch pitch'];
     tableau_W_glider=[tableau_W_glider W_glider'];
     
end

%%Linear Regression%%%%%
%Calcul coeff directeur

b1=[];
b=[];
p=[];
for i=1:min(size(tableau_pitch))
b1 = [b1 tableau_pitch(70:end-40,i)\tableau_W_glider(70:end-40,i)];
X = [ones(length(tableau_pitch(70:end-40,i)),1) tableau_pitch(70:end-40,i)];
b = [b X\tableau_W_glider(70:end-40,i)];
p = [p polyfit(tableau_pitch(70:end-40,i),tableau_W_glider(70:end-40,i),2)'];
end

% 
for i=1:min(size(tableau_W_glider))
 figure(i) 
 subplot(1,3,1)
 plot(tableau_pitch(:,i),-pi)
 title('pitch')
 ylabel('pression')
 subplot(1,3,2)
 plot(tableau_W_glider(:,i),-pi,'b')
 title('W\_glider')
 subplot(1,3,3)
 plot(tableau_pitch(:,i),tableau_W_glider(:,i),'+','Color','#4DBEEE')
 hold on 
 ylabel('W_glider')
 xlabel('Pitch')
 y1=b1(i)*tableau_pitch(70:end-40,i);
 X = [ones(length(tableau_pitch(70:end-40,i)),1) tableau_pitch(70:end-40,i)];
 y2 = X*b(:,i);
 y3 = polyval(p(:,i),tableau_pitch(70:end-40,i));
% plot(tableau_pitch(70:end-40,i),y1,'r','LineWidth',2);
 plot(tableau_pitch(70:end-40,i),y2,'--','LineWidth',2)
  plot(tableau_pitch(70:end-40,i),y3,'--','LineWidth',2,'Color','#7E2F8E')
% Rsq1 = 1 - sum((tableau_W_glider(70:end-40,i) - y1).^2)/sum((tableau_W_glider(70:end-40,i) - mean(tableau_W_glider(70:end-40,i))).^2);
 Rsq2 = 1 - sum((tableau_W_glider(70:end-40,i) - y2).^2)/sum((tableau_W_glider(70:end-40,i) - mean(tableau_W_glider(70:end-40,i))).^2);
 Rsq3 = 1 - sum((tableau_W_glider(70:end-40,i) - y3).^2)/sum((tableau_W_glider(70:end-40,i) - mean(tableau_W_glider(70:end-40,i))).^2);
 legend('W\_glider',['R\_affine= ' num2str(Rsq2)],['R\_polynome= ' num2str(Rsq3)])
 
end






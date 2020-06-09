function [explorer] = read_EXPLORER(tableau,explorer,Const)

%INPUT : tableau of Data 
         %explorer structure

 

%OUTPUT : explorer.dive           number of dive
%           explorer.lat          latitude
%           explorer.lon          longitude
%           explorer.time         time
%           explorer.depth        depth
%           explorer.oil          oil volume
%           explorer.pitch        pitch
%           explorer.pressure     pressure
%           explorer.temp         temperature
%           explorer.c            conductivity
%           explorer.s            salinity
%           explorer.dens         density


%%  one column for all dives 
explorer.dive = tableau(:,1);
explorer.lat = tableau(:,10)/100;
explorer.lon = tableau(:,11)/100;
explorer.time = tableau(:,8);
explorer.depth = tableau(:,9);
explorer.oil = tableau(:,13);
explorer.pitch = tableau(:,15);
explorer.pressure = tableau(:,17);
explorer.temp = tableau(:,18);
explorer.c = tableau(:,19);


%Remove NaN
to_ignore.lat=isnan(explorer.lat);
to_ignore.lon=isnan(explorer.lon);
to_ignore.time=isnan(explorer.time);
to_ignore.depth=isnan(explorer.depth);
to_ignore.oil=isnan(explorer.oil);
to_ignore.pitch=isnan(explorer.pitch);
to_ignore.pressure=isnan(explorer.pressure);
to_ignore.temp=isnan(explorer.temp);
to_ignore.c=isnan(explorer.c);

to_ignore.tot = to_ignore.lat+to_ignore.lon+to_ignore.time+to_ignore.depth+to_ignore.oil+to_ignore.pitch+to_ignore.pressure+to_ignore.temp+to_ignore.c; 


explorer.lat(to_ignore.tot ~= 0 )=[];
explorer.lon(to_ignore.tot ~= 0 )=[];
explorer.time(to_ignore.tot ~= 0 )=[];
explorer.depth(to_ignore.tot ~= 0 )=[];
explorer.oil(to_ignore.tot ~= 0 )=[];
explorer.pitch(to_ignore.tot ~= 0 )=[];
explorer.pressure(to_ignore.tot ~= 0 )=[];
explorer.temp(to_ignore.tot ~= 0 )=[];
explorer.dive(to_ignore.tot ~= 0 )=[];
explorer.c(to_ignore.tot ~= 0 )=[];

explorer.s = sw_salt(explorer.c*10/42.914,explorer.temp,explorer.pressure); 
 
explorer.dens = sw_dens(explorer.s,explorer.temp,explorer.pressure); 
 
explorer.size=length(explorer.lat); 
 

W_glider = zeros(1,explorer.size);
    for k=1:explorer.size-5
    W_glider(k) = (explorer.pressure(k)-explorer.pressure(k+5))/(explorer.time(k+5)-explorer.time(k));
    end
    for k=explorer.size-5:explorer.size
       W_glider(k)=NaN; 
    end
   W_glider=fillmissing(W_glider,'previous');
   
   explorer.W_glider = W_glider./Const.d2s;




explorer.pitch_filter=[];
explorer.W_glider_filter=[];
for j=explorer.first:explorer.last
    
   explorer1=by_dive(tableau,j,explorer,Const);
   %Filter pitch 
   explorer.pitch_filter = [explorer.pitch_filter explorer1.pitch'];
   
   %Filter W_glider 
   explorer.W_glider_filter = [explorer.W_glider_filter explorer1.W_glider];
end


%% SeaExplorer caracteristics 
 
explorer.M = 59; % Mass in kg 
explorer.V0 = 0.0573; %Volume 
explorer.alpha = 3*10^-10; %compressibility %Pa^-1 
explorer.beta =  1.09*10^-4; %coefficient d'expansion thermique  
explorer.T0=20; % température 
explorer.a=3; % angle d'attaque degrés  
explorer.S= 0.032365; % surface du glider  
explorer.Cd=0.1; % coefficient de trainée  

 
end


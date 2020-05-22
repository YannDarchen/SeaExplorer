function [explorer] = read_EXPLORER(tableau)


explorer.lat = tableau(:,10)/100;
explorer.lon = tableau(:,11)/100;
explorer.time = tableau(:,8);
explorer.depth = tableau(:,9);
explorer.oil = tableau(:,13);
explorer.pitch = tableau(:,15);
explorer.pressure = tableau(:,17);
explorer.temp = tableau(:,18);



%Remove NaN
to_ignore.lat=isnan(explorer.lat);
to_ignore.lon=isnan(explorer.lon);
to_ignore.time=isnan(explorer.time);
to_ignore.depth=isnan(explorer.depth);
to_ignore.oil=isnan(explorer.oil);
to_ignore.pitch=isnan(explorer.pitch);
to_ignore.pressure=isnan(explorer.pressure);
to_ignore.temp=isnan(explorer.temp);

to_ignore.tot = to_ignore.lat+to_ignore.lon+to_ignore.time+to_ignore.depth+to_ignore.oil+to_ignore.pitch+to_ignore.pressure+to_ignore.temp; 


explorer.lat(to_ignore.tot ~= 0 )=[];
explorer.lon(to_ignore.tot ~= 0 )=[];
explorer.time(to_ignore.tot ~= 0 )=[];
explorer.depth(to_ignore.tot ~= 0 )=[];
explorer.oil(to_ignore.tot ~= 0 )=[];
explorer.pitch(to_ignore.tot ~= 0 )=[];
explorer.pressure(to_ignore.tot ~= 0 )=[];
explorer.temp(to_ignore.tot ~= 0 )=[];


explorer.dens = sw_smow(explorer.temp);

explorer.size=length(explorer.lat);

explorer.M = 59; % Mass in kg
explorer.V0 = 0.051358; %Volume
explorer.alpha = 3*10^-6; %compressibility
explorer.beta =  1.09*10^-4; %coefficient d'expansion thermique 
explorer.T0=20; % température
explorer.a=3; % angle d'attaque degrés 
explorer.S= 0.032365; % surface du glider 
explorer.Cd=0.4; % coefficient de trainée 
end


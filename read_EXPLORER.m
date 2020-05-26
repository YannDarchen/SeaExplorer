function [explorer] = read_EXPLORER(tableau,explorer)

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


%Conversion conductivité -> salinité 
    %Coefficients 
    a0 = 0.008;
    a1 = -0.1692;
    a2 = 25.3851;
    a3 = 14.0941;
    a4 = -7.0261;
    a5 = 2.7081;
    b0 = 0.0005;
    b1 = -0.0056;
    b2 = -0.0066;
    b3 = -0.0375;
    b4 = 0.0636;
    b5 = -0.0144;
    c0 = 0.6766097;
    c1 = 0.0200564;
    c2 = 0.0001104259;
    c3 = -6.9698e-07;
    c4 = 0.000000001;
    k=0.0162;
    
D = ((explorer.c*10)/42.914)./(c0 + c1*explorer.temp + c2*(explorer.temp.^2)+c3*(explorer.temp.^3)+c4*(explorer.temp.^4));
explorer.s = a0 + a1*D.^0.5 + a2*D + a3*D.^1.5 + a4*D.^2 + a5*D.^2.5 + ...
            ((explorer.temp-15)./(1+k*(explorer.temp-15))).*(b0+b1*D.^0.5+b2*D+b3*D.^1.5+b4*D.^2+b5*D.^2.5);

explorer.dens = sw_dens(explorer.s,explorer.temp,explorer.pressure);

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


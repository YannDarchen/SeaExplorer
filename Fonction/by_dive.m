function [explorer] = by_dive(tableau,i,explorer,Const) 
 
%INPUT : tableau of Data  
%        number of Dive #i 
  
 
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
 
i_dive1 = tableau(:,1) == i; 
  tableau_dive = tableau(i_dive1,:); 
 
       
explorer.dive = tableau_dive(:,1); 
explorer.lat = tableau_dive(:,10)/100; 
explorer.lon = tableau_dive(:,11)/100; 
explorer.time = tableau_dive(:,8); 
explorer.depth = tableau_dive(:,9); 
explorer.oil = tableau_dive(:,13); 
explorer.pitch = tableau_dive(:,15); 
explorer.pressure = tableau_dive(:,17); 
explorer.temp = tableau_dive(:,18); 
explorer.c = tableau_dive(:,19); 
explorer.pitch=explorer.pitch;
 
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


   %Filter Pitch
  ind_moins = explorer.pitch <0;
   pitch_negatif = explorer.pitch(ind_moins);
   Mean_moins =mean(pitch_negatif);
   to_ign = explorer.pitch < Mean_moins -0.5;
  explorer.pitch(to_ign)=NaN;
   explorer.pitch = fillmissing(explorer.pitch,'next');
   
   

   %Filter W_glider 
   ind_neg = explorer.W_glider < 0;
   negatif = explorer.W_glider(ind_neg);
   [TF] = isoutlier(negatif,'percentiles',[3 100]);
   explorer.W_glider(TF)=NaN;
   explorer.W_glider = fillmissing(explorer.W_glider,'next');

 explorer.W_glider = smoothdata(explorer.W_glider,'SmoothingFactor',0.02);
 explorer.pitch = smoothdata(explorer.pitch,'SmoothingFactor',0.017);
  explorer.temp = smoothdata(explorer.temp,'SmoothingFactor',0.02);
end 
 


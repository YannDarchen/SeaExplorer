function [explorer,mx] = by_dive2(tableau,explorer,Const) 
 
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


%Find the longest dive
S =[];
s=1;
for j = 1:length(tableau)-1
    
   if tableau(j,1) == tableau(j+1,1)
       s=s+1;
   else 
       S=[S s];
       s=1;
   end
end

mx = max(S);
explorer.dive = NaN(mx,explorer.last_dive-explorer.first_dive+1); 
explorer.lat = NaN(mx,explorer.last_dive-explorer.first_dive+1);
explorer.lon =NaN(mx,explorer.last_dive-explorer.first_dive+1);
explorer.time = NaN(mx,explorer.last_dive-explorer.first_dive+1);
explorer.depth = NaN(mx,explorer.last_dive-explorer.first_dive+1);
explorer.oil = NaN(mx,explorer.last_dive-explorer.first_dive+1);
explorer.pitch = NaN(mx,explorer.last_dive-explorer.first_dive+1);
explorer.pressure = NaN(mx,explorer.last_dive-explorer.first_dive+1);
explorer.temp = NaN(mx,explorer.last_dive-explorer.first_dive+1);
explorer.c =NaN(mx,explorer.last_dive-explorer.first_dive+1);

a=0;
i=1;
for j = 1:length(tableau)-1
    
    if tableau(j,1) == tableau(j+1,1)
        explorer.lat(i,1+a) = tableau(j,10)/100;
        explorer.dive(i,1+a) = tableau(j,1); 
        explorer.lon(i,1+a) =tableau(j,11)/100;
        explorer.time(i,1+a) = tableau(j,8);
        explorer.depth(i,1+a) = tableau(j,9);
        explorer.oil(i,1+a) = tableau(j,13);
        explorer.pitch(i,1+a) = tableau(j,15);
        explorer.pressure(i,1+a) =tableau(j,17);
        explorer.temp(i,1+a) = tableau(j,18);
        explorer.c(i,1+a) = tableau(j,19);
        explorer.s(i,1+a) = sw_salt(explorer.c(i,1+a)*10/42.914,explorer.temp(i,1+a),explorer.pressure(i,1+a)); 
        explorer.dens(i,1+a) = sw_dens(explorer.s(i,1+a),explorer.temp(i,1+a),explorer.pressure(i,1+a)); 
    else
        explorer.lat(i,1+a)= tableau(j,10)/100;
        explorer.dive(i,1+a) = tableau(j,1); 
        explorer.lon(i,1+a) =tableau(j,11)/100;
        explorer.time(i,1+a) = tableau(j,8);
        explorer.depth(i,1+a) = tableau(j,9);
        explorer.oil(i,1+a) = tableau(j,13);
        explorer.pitch(i,1+a) = tableau(j,15);
        explorer.pressure(i,1+a) =tableau(j,17);
        explorer.temp(i,1+a) = tableau(j,18);
        explorer.c(i,1+a) = tableau(j,19);
        explorer.s(i,1+a) = sw_salt(explorer.c(i,1+a)*10/42.914,explorer.temp(i,1+a),explorer.pressure(i,1+a)); 
        explorer.dens(i,1+a) = sw_dens(explorer.s(i,1+a),explorer.temp(i,1+a),explorer.pressure(i,1+a)); 
        
        a=a+1;
        i=0;
    end
    i=i+1;
end

explorer.size=length(explorer.lat); 


W_glider = NaN(explorer.size,min(size(explorer.lat)));
    for k=1:min(size(explorer.lat)) % column
        for j = 1:explorer.size-5  %row
            W_glider(j,k) = (explorer.pressure(j,k)-explorer.pressure(j+5,k))/(explorer.time(j+5,k)-explorer.time(j,k));
        end
    end
    for k = 1:min(size(explorer.lat))
        for j=explorer.size-5:explorer.size
            W_glider(j,k)=NaN; 
        end
    end
   W_glider=fillmissing(W_glider,'previous');
   
   explorer.W_glider = W_glider./Const.d2s;

  explorer.pitch_filter=explorer.pitch;
   %Filter Pitch
   for k=1:min(size(explorer.lat))  %each column
       
       ind_moins = explorer.pitch_filter(:,k) <0;
       pitch_negatif = explorer.pitch_filter(ind_moins,k);
       [TF] = isoutlier(pitch_negatif,'percentiles',[4 100]);
       explorer.pitch_filter(TF,k)=NaN;
       explorer.pitch_filter = fillmissing(explorer.pitch_filter,'next');
    
   end
   explorer.W_glider_filter=explorer.W_glider;
   %Filter W_glider 
   for k=1:min(size(explorer.lat))
       
       ind_neg = explorer.W_glider_filter(:,k) < 0;
       negatif = explorer.W_glider_filter(ind_neg,k);
       [TF] = isoutlier(negatif,'percentiles',[3 100]);
       explorer.W_glider_filter(TF,k)=NaN;
       explorer.W_glider_filter = fillmissing(explorer.W_glider_filter,'next');
    end
   

end 
 


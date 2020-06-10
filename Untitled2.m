
A=0;
for i = 1:length(explorer.lat)-1

    dist(i)= distancelonlat(explorer.lat(i),explorer.lon(i),explorer.lat(i+1),explorer.lon(i+1),Const);
    A=A+dist(i);
    
end
ind=dist==0;
dist(ind)=[];
X = explorer.lat_interp;
Y = explorer.lon_interp;
Z = -explorer.depth_interp;
C = explorer.temp_interp;
Z(end)=NaN;
C(end)=NaN;
figure

fill3(X,Y,Z,C,'EdgeColor','w','LineWidth',1)
colorbar
hold on 
for j= 1:length(explorer.lat_interp)

    
    fill3(X(j),Y(j),Z(j),C(j),'EdgeColor','interp','Marker','o','MarkerFaceColor','flat')
    drawnow
end
hold off


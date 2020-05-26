x = linspace(1,10,100);
y = sin(x);
y(end) = NaN;
c = y;

figure
patch(x,y,c,'EdgeColor','w','Marker','o','MarkerFaceColor','flat');
colorbar;
for j=1:length(x)
    patch(x(j),y(j),c(j),'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
    drawnow
end

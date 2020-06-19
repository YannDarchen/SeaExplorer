figure()
subplot(3,1,1)
 pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww6);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
     colormap(blue_red_cmap)
    H=colorbar;
    caxis([-0.02 0.02])
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_diving')
subplot(3,1,2)
  pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww2);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
     colormap(blue_red_cmap)
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_all\_descent')
subplot(3,1,3)
   pcolor([explorer.first_dive:explorer.last_dive],-pi,tableau_ww1);   % utilisation de la fonction pcolor + "shading interp"
    shading interp
    colormap(blue_red_cmap)
    H=colorbar;
    ylabel(H,'          ','FontSize',12,'Rotation',0);
    grid on
    ax = gca;
    ax.Layer='Top';
    xlabel('Numéro de plongée')
    ylabel('- Pression (dbar)')
    title('Ww\_bydive\_descent')
    
    
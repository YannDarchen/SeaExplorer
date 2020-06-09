figure()
plot(explorer.time,explorer.pressure)
title('pression')
figure()
plot(explorer.time,explorer.dens)
title('densité')
figure()
plot(explorer.time,explorer.pitch_filter)
title('pitch')
figure()
plot(explorer.time,explorer.oil)
title('oil')
figure()
plot(explorer.time,explorer.temp)
title('temp')
figure()
subplot(2,1,1)
yyaxis left
plot(explorer.time,explorer.W_glider_filter)
hold on 
plot(explorer.time,W_model)
yyaxis right
plot(explorer.time,explorer.pressure)

hold off
subplot(2,1,2)
yyaxis left
plot(explorer.time,explorer.temp)
hold on 
yyaxis right
plot(explorer.time,-explorer.pressure)
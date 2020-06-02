function [cost] = cost(param,explorer)
% cost function to minimize mean water vertical velocities
[W_model] = flight_model(explorer.pressure,explorer.dens,explorer.pitch,explorer.oil,explorer.temp,param(1),param(2),param(3),explorer.M);
cost = nanmean(explorer.W_glider'.^2-W_model(1:end-5).^2);

end


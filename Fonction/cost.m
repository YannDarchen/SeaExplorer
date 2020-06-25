function [cost] = cost(param)
% cost function to minimize mean water vertical velocities
global pres dens pitch oil_vol Wglider temp mg 
[W_model] = flight_model(pres,dens,pitch,oil_vol,temp,param(1),param(2),param(3),mg);

%% Remove spikes
% ind = find(W_model<-0.3);
% W_model(ind)=NaN;
% Wglider(ind)=NaN;

%% Cost function
%cost = abs(nanmean((Wglider'.^2-W_model(1:end-5).^2)));

%cost= nanmean((Wglider'-W_model(1:end-5)).^2);

A=Wglider'.^2-W_model.^2;
cost=nansum(abs(A)); 


%A=(Wglider'-W_model(1:end-5)).^2;
%cost=sum(abs(A)); 
end


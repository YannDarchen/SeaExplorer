function [cost] = cost_wings(param)
% cost function to minimize mean water vertical velocities
global pres dens pitch oil_vol Wglider temp mg V0 eps Cd
[W_model] = flight_model2(pres,dens,pitch,oil_vol,temp,V0,eps,Cd,mg,param(1),param(2));

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


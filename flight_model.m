function [wg,U,att,Fg,Fb,Fl,Fd,att_deg] = flight_model(p,rho,pitch,dvbp,temp,Vg,eps,Cd0,mg)
%%%%%
% Glider flight model (adapted from Merckelbach et al JAOT 2010), Anthony Bosse (abolod@locean-ipsl.upmc.fr), September 2015
% INPUT variables :
% p : pressure
% rho : in situ density
% pitch : pitch angle
% dvbp : pumped oil volume
% tt : temperature
% Vg : glider volume (unknown, optimized)
% eps : glider compressibility (unknown, optimized) %Pa^-1
% Cd0 : drag coefficient (unknown, optimized)
% mg : glider mass

%OUTPUT variables : 
% U : glider velocity
% wg : glider vertical velocity
% att : attack angle
% Fg : gravity force
% Fb :  buoyancy force
% Fl : lift force
% Fd : drag force
%%%%%

% constant parameters for SLOCUM 
g = 9.81;
S = 0.1; %section hull 0.05
aw = 3.7;
ah = 2.4; % in debate
Cd1w = 0.78;
Cd1h = 2.1;

% Angles
pitch = pitch * (pi/180); % d2r

% Solve explicitely for att, angle of attack
%att = NaN(1,length(pitch));
%for i = 1:length(pitch)
%	p = pitch(i);
%if (isnan(pitch(i))==0)
%syms x1.8
%att2 = solve('(aw+ah)*x*tan(p) = Cd0+(Cd1w+Cd1h)*x^2','abs(x)<0.1'))
%	if (pitch(i)>0)
%	att(i) = att2(1);
%	end
%	if (pitch(i)<0)
%	att(i) = att2(2);
%	end
%end
%end
%att = sign(pitch)*(pi/60); % radians
% OR suppose att<<pitch and solve second order equation
%[att1,att2,delta] = solve_2deg((Cd1w+Cd1h).*ones(size(pitch))-(1./(cos(pitch).^2)).*(aw+ah),-(aw+ah)*tan(pitch),Cd0.*ones(size(pitch)));
[att1,att2] = solve_2deg((Cd1w+Cd1h).*ones(size(pitch)),-(aw+ah)*tan(pitch),Cd0.*ones(size(pitch)));
att = sign(pitch).*min(abs(att1),abs(att2)); % keep the smallest root 
att_deg = att*(180/pi);
% Vertical forces
Fg = mg*g;

alphat= 7.05e-5;% thermic compressibility
Fb = g*(rho).*(Vg*(1-eps*(p*10000)+alphat.*(temp-13.2))+dvbp/1000000);


% Glider speeds through water
U = real(sqrt( 2*(Fb-Fg)./((rho).*S.*(Cd0+(Cd1w+Cd1h).*(att.^2)).*((sin(pitch+att).^2 +...
    cos(pitch+att).^2)./sin(pitch+att))) )); % !! sign error in equation (13) from Merckelbach et al...
wg = U.*sin(att+pitch);
wg=smoothdata(wg,'movmedian','SmoothingFactor',0.09);
% Drag and lift forces
Fd = 0.5.*rho.*(Cd0+(Cd1w+Cd1h).*(att.^2)).*S.*U.*U;
Fl = 0.5.*rho.*(ah+aw).*att.*S.*U.*U; %positive upward

%%%%%%% external function to solve 2nd order equation
end 
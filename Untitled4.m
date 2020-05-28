


Vg = 0.0573;
eps = 3e-10;%Pa^-1
Cd0 = 0.4;
mg = 59;

[U,wg,att,Fg,Fb,Fl,Fd] = flight_model(explorer.pressure,explorer.dens,explorer.pitch,explorer.oil,explorer.temp,Vg,eps,Cd0,mg);
% aw = 3.7;
% ah = 2.4; % in debate
% Cd1w = 0.78;
% Cd1h = 2.1;
% pitch = explorer.pitch * (pi/180);
% syms x
% S_simple=[];
% for i=1:length(pitch)-9500
%     disp(i)
%     S_simple=[S_simple vpasolve((Cd1w+Cd1h)*x^2 - (aw+ah)*tan(pitch(i))*x + Cd0 ==0)];
% end
%S= vpasolve(eqn)
% T=[];
% for i=1:length(S_simple)
%    T=[T ;min(S_simple(2,i)]; 
% end
figure(3)
plot(att)
hold on 
plot (S_simple(2,:))
% syms x 
% a=2;
% b=3;
% c=0;
% eqn = a*x^2 + b*x +c ==0;
% S= solve(eqn)
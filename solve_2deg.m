function [v1,v2,delta] = solve_2deg(a,b,c)
delta = (b.^2-4.*a.*c);
v1 = (-b+delta.^0.5)./(2.*a);v2 = (-b-delta.^0.5)./(2.*a);
end
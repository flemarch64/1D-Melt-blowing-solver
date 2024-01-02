function [Visc] = Viscosity(A,B,G,K)
%Calculate the shear rate dependent viscosity using a Cross law
%   Input A and B as constant, G is the shear rate
Visc=(A/(1+((B*G)^0.85)))*exp(6843.83*((1/K)-(1/423.15)));
end


function [source_term_heat_exchange] = heat(alpha,n,df,Tair,Tf,rho,uf,upar,vis,lambda_air,h)
%Heat transfer with surrounding air
%   Return the source term of the heat exchange with surrounding air

%Reynolds number in each diameter
% if uf>upar
%     Red=(df*rho*((uf-upar)))/vis;
% else
%     Red=(df*rho*((upar-uf)))/vis;
% end
Red=(df*rho*(abs(uf-upar)))/vis;
%Kase-Matsuo correlation
Nud=alpha*(Red^n);

%Source term exchange
source_term_heat_exchange = h*pi*lambda_air*Nud*((Tair-Tf));
end


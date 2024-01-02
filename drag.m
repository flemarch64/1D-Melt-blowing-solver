function [Ex] = drag(alpha,n,df,uf,upar,vis,h,rho)
%Source term exchange with surrounding air
%   Return the source term of the momentum exchange with surrounding air

%Reynolds number in each diameter
% if uf>upar
    Red=(df*rho*(abs(uf-upar)))/vis;
% else
%     Red=(df*rho*(abs(upar-uf)))/vis;
% end
%Kase-Matsuo correlation
Cf=alpha/(Red^n);

%Ex=h*0.5*rho*pi*df*((Cf*(uf-upar)*(uf-upar)));
Ex=0.5*h*rho*pi*df*((Cf*abs(uf-upar)*(upar-uf)));
end


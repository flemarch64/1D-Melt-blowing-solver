%Global script 1D model
clear all;
%Inputs
Q=5.0e-5; %Fiber mass flow rate in kg/s
rho_f=860;
Qm=Q/rho_f; %volumetric flow rate in m3/s
D0=0.773e-3;%initial fiber diameter
A0=pi*D0^2*0.25;%Area of the fiber
u0=(Q/(A0*rho_f));%initial fiber velocity
Tn=363.15;%solidification temperature
zmin=0.53e-3;
zmax=0.015;
rho_air=1.225; %air density
mu_air=1.7894e-5; %air viscosity
lambda_air=0.0242; %air conductivity
vis=2.4; %Fiber zero-shear viscosity
g=9.81; %gravity
cond=0.157;%conductivity of the fiber
Cp0=750; %Specific heat of the fiber
n_drag=0.63; %Shambaugh correlation coefficient to adjust if needed

%Import data from air simulation in order to do the 1-way coupling
Tair=dlmread('temperature.csv');
Uair=dlmread('velocity.csv');

%Meshing of N points 
n=50;
h=(zmax-zmin)/(n-1);
Z=zmin:h:zmax;
Z=Z';

%Define vectors initial values of fiber area, velocity and source term
Af=ones(n,1);
Vf=ones(n,1);
Tf=ones(n,1);
Vp=ones(n,1);
Tp=ones(n,1);
B=ones(n,1);
q=((Q/n)*ones(n,1))/rho_f; %volumetric flow rate vector


%Datos temperatura fibra
T0=416.76;
A_vis=1.46e-07;
B_vis=7210.2;

%Solver beggining inputs 
iter=0;
maxiter=500;
iflag=1;

%Prescribe initial conditions for control volumes
for i=1:n
    Af(i)=A0;
    Vf(i)=u0;
    Tf(i)=T0;
    Ap(i)=Af(i);
    Vp(i)=Vf(i);
    Tp(i)=Tf(i);
    visc(i)=2.4;
end

%Solver initialization, without diameter-fiber velocity coupling
%Continuity equation
for i=1:n
    if i==1
        Af(i)=A0;
    else
        Af(i)=Qm/(Vf(i));
    end
end

%Vector initialization for TDMA
M=ones(n,1);
U=ones(n,1);
W=ones(n,1);
S=ones(n,1);
D=ones(n,1)*D0;
R1=zeros(n,1);
R2=zeros(n,1);
R1p=zeros(n,1);
R2p=zeros(n,1);
R1(1)=1;
R2(1)=1;
R1p(1)=1;
R2p(1)=1;

%Initialization
%Momentum equation

%Coefficients for TDMA M(i)*phi(i)-U(i)*phi(i+1)-W(i)*phi(i-1)=B(i)
%First order upwind
for i=1:n
    if i==1
        M(i)=1;
        U(i)=0;
        W(i)=0;
        B(i)=u0;
    elseif i==n
        M(i)=1;
        U(i)=0;
        W(i)=1;
        B(i)=0;
    else
        M(i)=(0.75*vis*(1/h))*(Af(i+1)+(2*Af(i))+Af(i-1))+(Q);
        U(i)=((0.75*vis*(1/h))*(Af(i+1)+Af(i)));
        W(i)=((0.75*vis*(1/h))*(Af(i)+Af(i-1)))+(Q);
        %Momentum source term evaluated at cell-centroid
        B(i)=(0.5*h*g*rho_f*(Af(i)));
    end
end


%TDMA algorithm
PHI=TDMA(M,U,W,B,n);
Vf=PHI;

%Continuity equation
for i=1:n
    if i==1
        Af(i)=A0;
    else
        Af(i)=Qm/(Vf(i));
    end
end

%Compute Fiber diameter
for i=2:n
    D(i)=diam(Af(i));
end

%Temperature equation

%Matrix for TDMA
O=ones(n,1);
Y=ones(n,1);
R=ones(n,1);

%Vector for shear rate
shear=ones(n,1);

%Coefficients for TDMA O(i)*Tf(i)-Y(i)*Tf(i+1)-R(i)*Tf(i-1)=S(i)
%First order upwind
for i=1:n
    if i==1
        O(i)=1;
        Y(i)=0;
        R(i)=0;
        S(i)=T0;
    elseif i==n
        O(i)=1;
        Y(i)=0;
        R(i)=0;
        S(i)=Tn;
    else 
        O(i)=(Q*Cp0)+(((cond/(2*h))*(Af(i+1)+(2*Af(i))+Af(i-1))));
        Y(i)=((cond/(2*h))*(Af(i+1)+Af(i)));
        R(i)=((cond/(2*h))*(Af(i)+Af(i-1)))+(Q*Cp0);
        S(i)=0;
    end
end

%TDMA algorithm
PHI2=TDMA(O,Y,R,S,n);
Tf=PHI2;
%end of fiber equations initialization without coupling with surrounding air


%under relaxation factors
alfa_source_temp=0.3;
alfa_source_momentum=0.5;

%Assignation of the first values used after initialization
Vp=Vf;
Tp=Tf;
Sp=S;
Bp=B;
%%End of Initialization


%Solver
while iflag==1
    
    %Continuity equation
    for i=1:n
        if i==1
            Af(i)=A0;
        else
            Af(i)=Qm/(Vp(i));
        end
    end
    
    %Non-Newtonian viscosity
    %The gradients are stored in cell centroids
    for i=1:n
        if i==1
            shear(i)=0;
        elseif i==n
            shear(i)=0;
        else
            shear(i)=(Vp(i+1)-Vp(i-1))/(2*h);
        end
    end
    %Calculate the non-newtonian viscosity
    for i=1:n
        visc(i)=Viscosity(2.4,2e-4,shear(i),Tp(i));
    end
    
    %Coefficients for TDMA M(i)*phi(i)-U(i)*phi(i+1)-W(i)*phi(i-1)=B(i)
    %1-way coupling with surrounding air
    %First order upwind
    for i=1:n
        if i==1
            M(i)=1;
            U(i)=0;
            W(i)=0;
            B(i)=u0;
        elseif i==n
            M(i)=1;
            U(i)=0;
            W(i)=1;
            B(i)=0;
        else
            M(i)=(visc(i)*(0.75/(h)))*(Af(i+1)+(2*Af(i))+Af(i-1))+(Q);
            U(i)=((visc(i)*(0.75/(h)))*(Af(i+1)+Af(i)));
            W(i)=((visc(i)*(0.75/(h)))*(Af(i)+Af(i-1)))+(Q);
            %Source term discretization with Shambaugh law Cf=alpha/(Red^n)
            %Momentum source term evaluated at cell-centroid
            B(i)=(0.5*h*g*rho_f*Af(i))+drag(1.24,n_drag,D(i),Vp(i),Uair(i),mu_air,h,rho_air);
        end
    end

    %TDMA algorithm
    PHI=TDMA(M,U,W,B,n);
    %Under relaxation
    Vf=Vp+alfa_source_momentum*(PHI-Vp);
    
    %Residuals
    R1=abs(B-(diag(M)+diag(-U(1:n-1),1)+diag(-W(2:end),-1))*Vf);
    
    %Continuity equation again to compute temperature field afterwards
    for i=1:n
        if i==1
            Af(i)=A0;
        else
            Af(i)=Qm/(Vf(i));
        end
    end

    %Compute Fiber diameter
    for i=2:n
        D(i)=diam(Af(i));
    end
    
    
    %Fiber temperature, 1-way coupling
    %First order upwind
    for i=1:n
        if i==1
            O(i)=1;
            Y(i)=0;
            R(i)=0;
            S(i)=T0;
        elseif i==n
            O(i)=1;
            Y(i)=0;
            R(i)=1;
            S(i)=0;
        else 
            O(i)=(((cond/(2*h))*(Af(i+1)+(2*Af(i))+Af(i-1))))+(Q*Cp0);
            Y(i)=((cond/(2*h))*(Af(i+1)+Af(i)));
            R(i)=((cond/(2*h))*(Af(i)+Af(i-1)))+(Q*Cp0);
            %Heat source term evaluated at cell-centroid
            S(i)=0.5*heat(4.14,0.805,D(i),Tair(i),Tp(i),rho_air,Vf(i),Uair(i),mu_air,lambda_air,h);
        end
    end
 
    %TDMA algorithm
    PHI2=TDMA(O,Y,R,S,n);
    %Under relaxation
    Tf=Tp+alfa_source_temp*(PHI2-Tp);

    %Residuals
    R2=S-((diag(O)+diag(-Y(1:n-1),1)+diag(-R(2:end),-1))*Tf);
  
    %Convergence criterion
    %Not any coma to print these values (monitoring residuals)
    for i=1:n
        error_velo=abs((sqrt(sum(R1.^2))-sqrt(sum(R1p.^2))))
        error_temp=abs((sqrt(sum(R2.^2))-sqrt(sum(R2p.^2))))
    end

    %Convergence criterion
    error=1.0e-20;
    if error_velo>error || error_temp>error
        iter=iter+1;
        Vp=Vf;
        Tp=Tf;
        Sp=S;
        Bp=B;
        R1p=R1;
        R2p=R2;
        iflag=1;
    else
        iflag=0;
    end
    if iter>maxiter
        break
    end
end


%graphics
subplot(2,1,1)
plot(Z,D,'r')
legend('EXP 1D model')
title('Diameter=f(Z) EXP')
xlabel('Z(m)')
ylabel('Diameter(m)')
grid on
subplot(2,1,2)
plot(Z,Tf,'r')
legend('EXP 1D model')
title('Fiber temperature=f(Z) EXP')
xlabel('Z(m)')
ylabel('Tf(K)')
grid on


function [phi] = TDMA(a,b,c,d,n)
%Thomas diagonalisation algorithm
%a(i)*phi(i)-b(i)*phi(i+1)-ci*phi(i-1)=d(i)
P=zeros(n,1);
Q=zeros(n,1);
%   %forward elimination
for i=1:n
    if i==1
        P(i)=b(i)/a(i);
        Q(i)=d(i)/a(i);
    else
        P(i)=b(i)/(a(i)-c(i)*P(i-1));
        Q(i)=(d(i)+c(i)*Q(i-1))/(a(i)-c(i)*P(i-1));
    end
end

%Backward substitution 
for i =n:-1:1
    if i==n
        phi(i)=Q(i);
    else
        phi(i)=P(i)*phi(i+1)+Q(i);
    end
end

%Resolution of the matrix system
L=diag(a)+diag(-b(1:n-1),1)+diag(-c(2:end),-1);
phi=inv(L)*d;
end


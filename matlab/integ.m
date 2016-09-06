function ans  = integ(alpha,Th,var)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%theta=(-Th):0.001:0.5*pi;
%plot(theta,integrand(theta));


ans=quadl(@integrand,(-Th),0.5*pi);

 function y=integrand(theta)
y=g(theta,alpha,Th).*exp((-abs(var)^(alpha/(alpha-1))).*V(theta,alpha,Th));
 end


function ansg=g(theta,alpha,Th)
ansg=(sin(alpha.*(Th+theta)-2.*theta)./sin(alpha.*(Th+theta)))-((alpha.*cos(theta).^2)./(sin(alpha.*(Th+theta)).^2));
end

function ansV=V(theta,alpha,Th)
       ansV=(clf(cos(alpha.*Th)).^(1/(alpha-1))).*((cos(theta)./sin(alpha.*(Th+theta))).^(alpha/(alpha-1))).*(cos(alpha.*Th+(alpha-1).*theta)./cos(theta)); 
    end





end

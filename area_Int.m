%This function computes the Area of the general Case of an 
%angular partition of a polygon 

function y = area_Int(omega,xd,yd,beta1,xb,yb,x0,y0,phi,theta)

%y = ( ((tan(omega).*(x0-xd)+yd-y0).^2).*(tan(theta+phi)-tan(theta)))./(2*(tan(theta+phi)-tan(omega)).*(tan(theta)-tan(omega))) + ( ((tan(beta1).*(x0-xb)+yb-y0).^2).*(tan(theta)-tan(theta+phi)))./(2*(tan(theta+phi)-tan(beta1)).*(tan(theta)-tan(beta1)));


%y = ( ((sin(omega).*(x0-xd)+cos(omega).*(yd-y0)).^2).*(tan(theta+phi)-tan(theta)))./(2*(tan(theta+phi).*cos(omega)-sin(omega)).*(tan(theta).*cos(omega)-sin(omega))) + ( ((sin(beta1).*(x0-xb)+cos(beta1).*(yb-y0)).^2).*(tan(theta)-tan(theta+phi)))./(2*(tan(theta+phi).*cos(beta1)-sin(beta1)).*(tan(theta).*cos(beta1)-sin(beta1)));

%y = ( ((sin(omega).*(x0-xd)+cos(omega).*(yd-y0)).^2).*(sin(theta+phi).*cos(theta)-sin(theta).*cos(theta+phi)))./(2*(sin(theta+phi).*cos(omega)-cos(theta+phi).*sin(omega)).*(sin(theta).*cos(omega)-cos(theta).*sin(omega))) + ( ((sin(beta1).*(x0-xb)+cos(beta1).*(yb-y0)).^2).*(sin(theta).*cos(theta+phi)-sin(theta+phi).*cos(theta)))./(2*(sin(theta+phi).*cos(beta1)-cos(theta+phi).*sin(beta1)).*(sin(theta).*cos(beta1)-cos(theta).*sin(beta1)));


y =  ((sin(omega).*(x0-xd)+cos(omega).*(yd-y0)).^2).*(sin(phi))./(2*sin(theta+phi-omega).*sin(theta-omega)) -  ((sin(beta1).*(x0-xb)+cos(beta1).*(yb-y0)).^2).*sin(phi)./(2*sin(theta+phi-beta1).*sin(theta-beta1));


for i=1:length(phi)
    if(phi(i)<10^-15)
        y(i)=0;
    end
end


%k = length(theta);

%for i=1:k

%y(i) =( ((sin(omega(i)).*(x0-xd(i))+cos(omega(i)).*(yd(i)-y0)).^2).*(tan(theta(i)+phi(i))-tan(theta(i))))./(2*(tan(theta(i)+phi(i)).*cos(omega(i))-sin(omega(i))).*(tan(theta(i)).*cos(omega(i))-sin(omega(i)))) + ( ((sin(beta1(i)).*(x0-xb(i))+cos(beta1(i)).*(yb(i)-y0)).^2).*(tan(theta(i))-tan(theta(i)+phi(i))))./(2*(tan(theta(i)+phi(i)).*cos(beta1(i))-sin(beta1(i))).*(tan(theta(i)).*cos(beta1(i))-sin(beta1(i))));

%endfor

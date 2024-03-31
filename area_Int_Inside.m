function y = area_Int_Inside(omega,xd,yd,x0,y0,phi,theta)
    y = ((tan(omega).*(x0-xd)+yd-y0).^2).*sin(phi).*(cos(omega).^2)./(2*sin(theta+phi-omega).*sin(theta-omega));
end
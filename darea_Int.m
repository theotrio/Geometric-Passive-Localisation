function darea = darea_Int(omega,xd,yd,beta1,xb,yb,x0,y0,phi,theta)
    
    darea = (((tan(omega).*(x0-xd)+yd-y0).^2).*cos(omega).*cos(omega))./(2*sin(theta+phi-omega).*sin(theta+phi-omega))- (((tan(beta1).*(x0-xb)+yb-y0).^2).*cos(beta1).*cos(beta1))./(2*sin(theta+phi-beta1).*sin(theta+phi-beta1));

end
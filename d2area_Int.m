function d2area = d2area_Int(omega,xd,yd,beta1,xb,yb,x0,y0,phi,theta)
    
    d2area = (((tan(beta1).*(x0-xb)+yb-y0).^2).*(cos(beta1).^2).*sin(2*(theta+phi-beta1)))./(4*sin(theta+phi-beta1).^4) - (((tan(omega).*(x0-xd)+yd-y0).^2).*(cos(omega).^2).*sin(2*(theta+phi-omega)))./(4*sin(theta+phi-omega).^4);
    
end
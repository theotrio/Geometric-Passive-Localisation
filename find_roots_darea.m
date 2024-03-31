function phi_roots = find_roots_darea(omega,xd,yd,beta1,xb,yb,x0,y0,theta)

    d1 = (tan(omega).*(x0-xd)+yd-y0);
    d2 = (tan(beta1).*(x0-xb)+yb-y0);
    
    phi_roots = zeros(1,2);
    
    phi_roots(1) = atan((d2*cos(beta1)*sin(theta-omega)-d1*cos(omega)*sin(theta-beta1))/(d1*cos(omega)*cos(theta-beta1)-d2*cos(beta1)*cos(theta-omega)));
    phi_roots(2) = atan(-(d2*cos(beta1)*sin(theta-omega)+d1*cos(omega)*sin(theta-beta1))/(d1*cos(omega)*cos(theta-beta1)+d2*cos(beta1)*cos(theta-omega)));

end
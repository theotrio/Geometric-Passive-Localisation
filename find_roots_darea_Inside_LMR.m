function phi_roots = find_roots_darea_Inside_LMR(omega,xd,yd,x0,y0,theta)

    m = length(omega);

    d1 = (tan(omega(1))*(x0-xd(1))+yd(1)-y0);
    d2 = (tan(omega(m))*(x0-xd(m))+yd(m)-y0);
    
    phi_roots = zeros(1,2);
    
    phi_roots(1) = atan((d2*cos(omega(m))*sin(theta(1)-omega(1))-d1*cos(omega(1))*sin(theta(m)-omega(m)))/(d1*cos(omega(1))*cos(theta(m)-omega(m))-d2*cos(omega(m))*cos(theta(1)-omega(1))));
    phi_roots(2) = atan(-(d2*cos(omega(m))*sin(theta(1)-omega(1))+d1*cos(omega(1))*sin(theta(m)-omega(m)))/(d1*cos(omega(1))*cos(theta(m)-omega(m))+d2*cos(omega(m))*cos(theta(1)-omega(1))));

end
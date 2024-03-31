function phi_roots = find_roots_darea_inside_LR(omega,xd,yd,x0,y0,theta,phi)

    m = length(omega);

    dleft = (tan(omega(1))*(x0-xd(1))+yd(1)-y0);
    dright = (tan(omega(m))*(x0-xd(m))+yd(m)-y0);
    
    %D = (dleft*cos(omega(1)))/dright*cos(omega(m));
    
    phi_roots = zeros(1,2);
    
    %phi_roots(1) = atan((D*sin(theta-omega(1))-sin(2*theta+phi-omega(m)))/(D*cos(theta-omega(1))-cos(2*theta+phi-omega(m))));
    %phi_roots(2) = atan((D*sin(theta-omega(1))+sin(2*theta+phi-omega(m)))/(D*cos(theta-omega(1))+cos(2*theta+phi-omega(m))));
    
    phi_roots(1) = atan( (dleft*cos(omega(1))*sin(theta-omega(1))-dright*cos(omega(m)))*sin(2*theta+phi-omega(m)))/(dleft*cos(omega(1)*cos(theta-omega(1))-dright*cos(omega(m))*cos(2*theta+phi-omega(m))) );
    phi_roots(2) = atan( (dleft*cos(omega(1))*sin(theta-omega(1))+dright*cos(omega(m)))*sin(2*theta+phi-omega(m)))/(dleft*cos(omega(1)*cos(theta-omega(1))+dright*cos(omega(m))*cos(2*theta+phi-omega(m))) );

end
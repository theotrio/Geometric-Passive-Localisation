function y = newton_equation_LMR(section_phi,omega,xd,yd,beta1,xb,yb,x0,y0,phi,section_theta,theta)
    m = length(section_theta);
    Al = area_Int(omega(1),xd(1),yd(1),beta1(1),xb(1),yb(1),x0,y0,theta+phi-section_theta(1),section_theta(1)); 
    Ar = area_Int(omega(m),xd(m),yd(m),beta1(m),xb(m),yb(m),x0,y0,theta-section_theta(m),section_theta(m));
    Mr = sum(area_Int(omega(2:m),xd(2:m),yd(2:m),beta1(2:m),xb(2:m),yb(2:m),x0,y0,section_phi(2:m),section_theta(2:m)));
    y = Al-Ar+Mr;
end
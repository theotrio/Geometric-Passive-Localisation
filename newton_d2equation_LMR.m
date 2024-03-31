function y = newton_d2equation_LMR(omega,xd,yd,beta1,xb,yb,x0,y0,phi,section_theta,theta)
    m = length(section_theta);
    
    Al = d2area_Int(omega(1),xd(1),yd(1),beta1(1),xb(1),yb(1),x0,y0,theta+phi-section_theta(1),section_theta(1)); 
    Ar = d2area_Int(omega(m),xd(m),yd(m),beta1(m),xb(m),yb(m),x0,y0,theta-section_theta(m),section_theta(m));

    y = Al-Ar;
end
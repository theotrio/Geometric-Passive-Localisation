% This function computes the area of 

function y = newton_equation_LR(section_phi,omega,xd,yd,beta1,xb,yb,x0,y0,phi,section_theta,theta)

    Al = area_Int(omega(1),xd(1),yd(1),beta1(1),xb(1),yb(1),x0,y0,theta+phi-section_theta(1),section_theta(1)); 
    Ar = area_Int(omega(2),xd(2),yd(2),beta1(2),xb(2),yb(2),x0,y0,theta-section_theta(2),section_theta(2));
    R = area_Int(omega(2),xd(2),yd(2),beta1(2),xb(2),yb(2),x0,y0,section_phi,section_theta(2));
    
    y = Al-Ar+R;

end
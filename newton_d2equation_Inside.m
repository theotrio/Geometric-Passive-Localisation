function y = newton_d2equation_Inside(omega,xd,yd,x0,y0,phi,theta)
    m = length(xd);
    
    ar = (((tan(omega(m))*(x0-xd(m))+yd(m)-y0)^2)*(cos(omega(m))^2)*sin(2*(theta(m)+phi(m)-omega(m))))/(4*sin(theta(m)+phi(m)-omega(m))^4); 
    al = (((tan(omega(1))*(x0-xd(1))+yd(1)-y0)^2)*(cos(omega(1))^2)*sin(2*(theta(1)+phi(1)-omega(1))))/(4*sin(theta(1)+phi(1)-omega(1))^4);
    y = ar-al;
end
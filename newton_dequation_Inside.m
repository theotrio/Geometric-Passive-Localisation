function y = newton_dequation_Inside(omega,xd,yd,x0,y0,phi,theta)
    m = length(xd);

    al = (((tan(omega(1))*(x0-xd(1))+yd(1)-y0)^2)*cos(omega(1))*cos(omega(1)))/(2*sin(theta(1)+phi(1)-omega(1))*sin(theta(1)+phi(1)-omega(1)));
    ar = (((tan(omega(m))*(x0-xd(m))+yd(m)-y0)^2)*cos(omega(m))*cos(omega(m)))/(2*sin(theta(m)+phi(m)-omega(m))*sin(theta(m)+phi(m)-omega(m)));
    y = al-ar;
end
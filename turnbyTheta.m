function y = turnbyTheta(x , theta)

  A = [cos(theta) -sin(theta);sin(theta) cos(theta)];
  
  if(size(x,1)==2)
    y = A*x;
  elseif(size(x,2)==2)
    y = A*x';
  end
  
end

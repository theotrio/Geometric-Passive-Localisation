function [A B] = turnbyTheta_Matrix(X, Y , theta)

  [n m] = size(theta);  
  A = zeros(n,m);
  B = zeros(n,m);

  for i=1:n
	for j=1:m
		x = [X(i,j) ; Y(i,j)];
		Z = turnbyTheta(x,theta(i,j));
		A(i,j) = Z(1);
		B(i,j) = Z(2);
    end
  end
  
end

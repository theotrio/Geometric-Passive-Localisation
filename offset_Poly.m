% This algorithm calculates the offset of a polygon P with vertices arranged 
% in a counter clockwise order. The algorithm enlarges or shrinks P up to a 
% distance d.
% The idea of the algorithm is to translate the edges of P to distance d
% and then recaclulate the new vertices as the intersection of the
% translated edges.
% Input: A set of the vertices of P (counterclockwise order), a distance d.
% Output: the vertices of the offset of the polygon in counter clockwise
% order.
function newVertices = offset_Poly(vertices,d,tol)

if(nargin<3)
    tol=eps;
end

 n = size(vertices,2);

a = zeros(1,n);
%Stores the edges of the polygon as lines y = ax+b where lines(1,k) = a_k
%and lines(2,k) = b_k 
lines = zeros(2,n);
%The vertices of the offset of the polygon
newVertices = zeros(2,n);

a(1) = (vertices(2,1) - vertices(2,2)) / (vertices(1,1)-vertices(1,2));
a(1) = atan(a(1));
if(abs( a(1) - pi/2 )<=tol || abs( a(1) + pi/2 )<=tol )
    if(vertices(1,1)<vertices(1,cycount(1,n,2)))
        lines(:,1) = [abs(a(1)); vertices(1,1)-d];
    else
        lines(:,1) = [abs(a(1)); vertices(1,1)+d];
    end
else
    if(tan(a(1))*(vertices(1,3)-vertices(1,1))+vertices(2,1)>vertices(2,3))
        v = turnbyTheta(vertices(:,1),-a(1));
        v(2) = v(2)+d;
        v = turnbyTheta(v,a(1));
        lines(:,1) = [a(1); v(2)-tan(a(1))*v(1)];
    else
        v = turnbyTheta(vertices(:,1),-a(1));
        v(2) = v(2)-d;
        v = turnbyTheta(v,a(1));
        lines(:,1) = [a(1); v(2)-tan(a(1))*v(1)];
    end
end

for i=2:n
    a(i) = (vertices(2,i) - vertices(2,cycount(i,n))) / (vertices(1,i)-vertices(1,cycount(i,n)));
    a(i) = atan(a(i));
    if(abs( a(i) - pi/2 )<=tol || abs( a(i) + pi/2 )<=tol )
        if(vertices(1,i)<vertices(1,cycount(i,n,2)))
            lines(:,i) = [abs(a(i)); vertices(1,i)-d];
        else
            lines(:,i) = [abs(a(i)); vertices(1,i)+d];
        end
    else
        if(tan(a(i))*(vertices(1,cycount(i,n,2))-vertices(1,i))+vertices(2,i)>vertices(2,cycount(i,n,2)))
            v = turnbyTheta(vertices(:,i),-a(i));
            v(2) = v(2)+d;
            v = turnbyTheta(v,a(i));
            lines(:,i) = [a(i); v(2)-tan(a(i))*v(1)];
        else
            v = turnbyTheta(vertices(:,i),-a(i));
            v(2) = v(2)-d;
            v = turnbyTheta(v,a(i));
            lines(:,i) = [a(i); v(2)-tan(a(i))*v(1)];
        end 
    end
    
    if(abs(lines(1,i-1) - pi/2)<=tol || abs(lines(1,i) - pi/2)<=tol)
        if(abs(lines(1,i-1) - pi/2)<=tol)
            newVertices(1,i) = lines(2,i-1);
            newVertices(2,i) = tan(lines(1,i))*newVertices(1,i)+lines(2,i);
        else
            newVertices(1,i) = lines(2,i);
            newVertices(2,i) = tan(lines(1,i-1))*lines(2,i)+lines(2,i-1);
        end
    else
        newVertices(1,i) = (lines(2,i)-lines(2,i-1))/ (tan(lines(1,i-1))-tan(lines(1,i)));
        newVertices(2,i) = tan(lines(1,i))*newVertices(1,i)+lines(2,i);
    end
end

if(abs(lines(1,n) - pi/2)<=tol || abs(lines(1,1) - pi/2)<=tol)
    if(abs(lines(1,n) - pi/2)<=tol)
        newVertices(1,1) = lines(2,n);
        newVertices(2,1) = tan(lines(1,1))*lines(2,n)+lines(2,1);
    else
        newVertices(1,1) = lines(2,1);
        newVertices(2,1) = tan(lines(1,n))*lines(2,1)+lines(2,n);
    end
else
    newVertices(1,1) = (lines(2,1)-lines(2,n))/ (tan(lines(1,n))-tan(lines(1,1)));
    newVertices(2,1) = tan(lines(1,1))*newVertices(1,1)+lines(2,1);
end

end
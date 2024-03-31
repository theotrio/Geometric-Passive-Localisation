% This functions computes the intersection between a polygon and a sector.
% Input: polygon is a set of the polygon's points, Al,Bl,Ar,Br the lines of
% the sector.
% Output: The points of the polygon's intersection.
function [npoly,pindex] = Intersect_Poly_Sector(polygon,Al,Ar,c,Bl,Br,su,tol)

if(nargin<5)
    %su shows if su=1 then the sensor has stumbled upon the emitter
    su=0;
end

if(su==1)
    %If su=1 then form a small polygon around the centre
    %cen = [Bl;Br];
    npoly = [c+[1,0]*tol,c+[0,-1]*tol,c+[-1,0]*tol,c+[0,-1]*tol];
    pindex = size(npoly,2);
    
else
    %We need to check if the centre of the sector is inside the polygon
    %We do that by checking if point c belongs to the convex hull of
    %the set polygon union c.
    p2 = [polygon c];

    p3 = ConvHull(p2);
    n = size(p3,2);
    is_outside=0;

    for i=1:n
        point_id = [abs(p3(1,i)-c(1))<=eps, abs(p3(2,i)-c(2))<=eps];
        check = and(point_id(1),point_id(2));
        if(check==1)
            is_outside=1;
        end
    end

    if(is_outside==1)
        [npoly,pindex] = Intersect_Poly_Sec1(polygon,Al,Ar,c);
    else
        [npoly,pindex] = Intersect_From_Inside2(polygon,Al,Ar,c);
    end

end

end
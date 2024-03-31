% This functions computes the intersection between a polygon and a sector.
% Input: polygon is a set of the polygon's points, Al,Bl,Ar,Br the lines of
% the sector.
% Output: The points of the polygon's intersection.
function [y,theta] = find_Maximum_Intersection(polygon,c,phi)

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
        [y,theta] = Max_Intersect_Out_Poly(polygon,c,phi);
    else
        [y,theta] = Max_Intersection_Inside_Poly(polygon,c,phi);
    end

end
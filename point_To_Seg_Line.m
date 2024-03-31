%Computes the distance between a a point and a line segment
%Input: 2x1 points, line1 and line2 are the two extreme points 
%of the segment
function [d,point] = point_To_Seg_Line(pt,line1,line2)
    
    vec1 = line2-line1;
    vec2 = pt-line1;
    
    point = line1+((vec1'*vec2)/(norm(vec1)^2))*(vec1);
    if(norm(point-line1)<=norm(vec1) && norm(point-line2)<=norm(vec1))
        d = norm(point-pt);
    else
        if(norm(line1-pt)<=norm(line2-pt))
            d = norm(line1-pt);
            point = line1;
        else
            d = norm(line2-pt);
            point = line2;
        end
    end
    
end
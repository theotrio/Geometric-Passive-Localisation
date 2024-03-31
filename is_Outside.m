function y = is_Outside(polygon,point)

    p2 = [polygon point];

    p3 = ConvHull(p2);
    n = size(p3,2);
    y=0;

    for i=1:n
        point_id = [abs(p3(1,i)-point(1))<=eps, abs(p3(2,i)-point(2))<=eps];
        check = and(point_id(1),point_id(2));
        if(check==1)
            y=1;
        end
    end
    
end
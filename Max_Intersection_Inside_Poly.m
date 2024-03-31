function [area,theta] = Max_Intersection_Inside_Poly(polygon,c,phi)

    n = length(polygon(1,:));

	SetPoints = polygon;
	x0 = c(1);
	y0 = c(2);
    
    %Compute the gradient lines between the points
    omegas = zeros(1,n);
    for i=1:n
        omegas(i) = atan2((SetPoints(2,i)-SetPoints(2,mod(i,n)+1)),(SetPoints(1,i)-SetPoints(1,mod(i,n)+1)));
    end

	%Arranging the Points of the Polygon from left to right
	xy = ones(2,n);
	xy(1,:) = x0*ones(1,n);
	xy(2,:) = y0*ones(1,n);

	points = SetPoints - xy;
    order = atan2(points(2,:),points(1,:));
    
    order = mod(order,2*pi);
    
    section_thetas = zeros(1,n+1);

	%Sort the points according to their angular position
	[section_thetas(1:n),order] = sort(order,'descend');
    %[section_thetas(1:n),order] = sort(order);
    
    section_thetas(n+1) = section_thetas(1)-2*pi;
    
    section_phis = zeros(1,length(section_thetas)-1);
    
    for i=1:length(section_thetas)-1
        section_phis(i) = section_thetas(i)-section_thetas(i+1);
    end
    
    [thetas,L,R] = just_merge(section_thetas,section_thetas-phi);
        
    xds = SetPoints(1,order);
    yds = SetPoints(2,order);
    omeg = omegas(order);

    [area,theta] = approximate_inner_direction(thetas,L,R,phi,section_thetas,section_phis,xds,yds,omeg,x0,y0);

end
% This functions computes the maximum intersection between a polygon P, and 
% a rotating sector with inner angle phi and a centre c which is outside of 
% P.
% Input: polygon is a set of the polygon's points, c is a point which the 
% centre of the rotating sector and phi is the inner angle of the secotr
% Output: The points of the polygon's intersection in counterclockwise order.

function [y, theta] = Max_Intersect_Out_Poly(polygon,c,phi)

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
    
    %Find the relative position of the polygon to the centre c
    %This will help to order the angular positions of the polygon's 
    %vertices either from -pi to pi or from 0 to 2*pi to adress the problem
    %that the angular positions are not continuous.
    dist = vecnorm(points);
    [~,idm] = min(dist);
    
    %Find the closest polygon's edge to the centre
    d1 = point_To_Seg_Line(c,SetPoints(:,idm),SetPoints(:,cycount(idm,n,1)));
    d2 = point_To_Seg_Line(c,SetPoints(:,idm),SetPoints(:,cycount(idm,n,-1)));
    
    if(d1==d2)
        if(x0<=SetPoints(1,idm))
            sl=1;
            stan=1;
        else
            sl=-1;
            stan=1;
        end
    else
        if(d1<d2)
            %idm2 = cycount(idm,n,1);
            phstan = omegas(idm);
        else
            idm2 = cycount(idm,n,-1);
            phstan = omegas(idm2);
        end
        
        stan = tan(phstan);
    
        if( abs(phstan - pi/2) <=10^-14 || abs(phstan+pi/2)<=10^-14)
            [~,idma] = max(dist);
            if(SetPoints(1,idma)>=x0)
                sl=1;
                stan = 1;
            else
                sl=1;
                stan = -1;
            end
        else
            y = stan*( x0-SetPoints(1,idm)) + SetPoints(2,idm);
            sl = (y0-y)./abs(y0-y);
            stan = stan/abs(stan);
        end
    end
    
    %If the closest edge has tangent>0 and x0 belongs in the upper plane
    %the line of the edge defines then the polygon is to the right of the
    %centre hence the oreder from -pi to pi is the correct one.
    if(sl*stan==-1)
        order = mod(order,2*pi);
    end

	%Sort the points according to their angular position
	[section_thetas,order] = sort(order,'descend');
    
    section_phis = zeros(1,length(section_thetas)-1);
    
    for i=1:length(section_thetas)-1
        section_phis(i) = section_thetas(i)-section_thetas(i+1);
    end
    
    [thetas,L,R] = just_merge(section_thetas,section_thetas-phi);

	%The sector doesn't intersect with the polygon
    if(phi-section_thetas(1)+section_thetas(n)>=10^-14)
		theta = section_thetas(n);
        y = area_Poly(SetPoints);
    else
        %Compute the line segments of the polygon 
        betaw = zeros(1,n); %Beta is the gradient of the below segments
        xbw = zeros(1,n);
        ybw = zeros(1,n);
        omegw = zeros(1,n); %Omega is the gradient of the upper segments
        xdw = zeros(1,n);
        ydw = zeros(1,n);

        betaw(1) = omegas(order(1));
        xbw(1) = SetPoints(1,order(1));
        ybw(1) = SetPoints(2,order(1));
        current_down = mod(order(1),n)+1;

        omegw(1) = omegas(cycount(order(1),n,-1));
        xdw(1) = SetPoints(1,order(1));
        ydw(1) = SetPoints(2,order(1));

        for j=2:n-1
            if(order(j)==current_down)
                    betaw(j) = omegas(order(j));
                    xbw(j) = SetPoints(1,order(j));
                    ybw(j) = SetPoints(2,order(j));
                    omegw(j) = omegw(j-1);
                    xdw(j) = xdw(j-1);
                    ydw(j) = ydw(j-1);
                    current_down = mod(order(j),n)+1;
            else
                    betaw(j) = betaw(j-1);
                    xbw(j) = xbw(j-1);
                    ybw(j) = ybw(j-1);
                    omegw(j) = omegas(cycount(order(j),n,-1));
                    xdw(j) = SetPoints(1,order(j));
                    ydw(j) = SetPoints(2,order(j));
            end
        end
        
        [y,theta] = approximate_direction1(thetas,L,R,phi,section_thetas,section_phis,xdw,ydw,omegw,xbw,ybw,betaw,x0,y0);
        
    end
    
    
end
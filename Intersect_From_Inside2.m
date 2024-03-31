%It is assumed that Ar and Al stem from atan2 so they range from -pi tp pi
function [npoly,pindex] = Intersect_From_Inside2(polygon,Al,Ar,c)

    n = length(polygon(1,:));

	SetPoints = polygon;
	x0 = c(1);
	y0 = c(2);
    tol = eps;
    
    if(mod(Al,2*pi)>=pi)
        anglesl = mod(Al,2*pi)-2*pi;
    else
        anglesl = Al;
    end
    if(mod(Ar,2*pi)<=pi)
        anglesr = mod(Ar,2*pi);
    else
        anglesr = Ar;
    end
    
    %Compute the gradient of the polygon edges
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
    
    angles = [order,anglesl,anglesr];
    k = length(angles);
    
    if(abs(anglesl-anglesr)>=pi)
        angles = mod(angles,2*pi);
    end
    
    [~,order] = sort(angles,'descend');
    
    %find the vertices inside the sector
    start=0;
    stop=0;
    
    for i=1:k
        if(order(i)==k-1)
            start = i;
        end
        if(order(i)==k)
            stop = i;
        end
    end
    
    %Compute the intersection points between the polygon edges and the sector 
    %In this case both lines intersect the same edge of the polygon
    if(stop-start==1)
        %Check if there is a vertical line
        edge_pol_check = or( abs(mod(omegas(order(cycount(stop,k,1))),2*pi)-pi/2)<=tol,abs(mod(omegas(order(cycount(stop,k,1))),2*pi)-3*pi/2)<=tol);
        Ar_check = or(abs(mod(anglesr,2*pi)-pi/2)<=tol,abs(mod(anglesr,2*pi)-3*pi/2)<=tol);
        Al_check = or(abs(mod(anglesl,2*pi)-pi/2)<=tol,abs(mod(anglesl,2*pi)-3*pi/2)<=tol);
        check_any_vertical = or(edge_pol_check,or(Ar_check,Al_check));
        if( check_any_vertical == 1 )
            if(edge_pol_check==1)
                x1 = SetPoints(1,order(cycount(stop,k,1)));
                y1 = tan(anglesl)*(x1-x0)+y0;
                y2 = tan(anglesr)*(x1-x0)+y0;
                points = [c,[x1;y1],[x1;y2]];
                npoly = ConvHull(points);
                pindex = size(npoly,2);
            end
            if(Ar_check==1)
                x1 = (tan(anglesl)*x0-tan(omegas(order(cycount(stop,k,1))))*SetPoints(1,order(cycount(stop,k,1)))+SetPoints(2,order(cycount(stop,k,1)))-y0)/(tan(anglesl)-tan(omegas(order(cycount(stop,k,1)))));
                y1 = tan(anglesl)*(x1-x0)+y0;
                x2 = x0;
                y2 = tan(omegas(order(cycount(stop,k,1))))*(x2-SetPoints(1,order(cycount(stop,k,1))))+SetPoints(2,order(cycount(stop,k,1)));
                points = [c,[x1;y1],[x2;y2]];
                npoly = ConvHull(points);
                pindex = size(npoly,2);
            end
            if(Al_check==1)
                x1 = x0;
                y1 = tan(omegas(order(cycount(stop,k,1))))*(x1-SetPoints(1,order(cycount(stop,k,1))))+SetPoints(2,order(cycount(stop,k,1)));
                x2 = (tan(anglesr)*x0-tan(omegas(order(cycount(stop,k,1))))*SetPoints(1,order(cycount(stop,k,1)))+SetPoints(2,order(cycount(stop,k,1)))-y0)/(tan(anglesr)-tan(omegas(order(cycount(stop,k,1)))));
                y2 = tan(anglesr)*(x2-x0)+y0;
                points = [c,[x1;y1],[x2;y2]];
                npoly = ConvHull(points);
                pindex = size(npoly,2);
            end
        else
            x1 = (tan(anglesl)*x0-tan(omegas(order(cycount(stop,k,1))))*SetPoints(1,order(cycount(stop,k,1)))+SetPoints(2,order(cycount(stop,k,1)))-y0)/(tan(anglesl)-tan(omegas(order(cycount(stop,k,1)))));
            y1 = tan(anglesl)*(x1-x0)+y0;
            x2 = (tan(anglesr)*x0-tan(omegas(order(cycount(stop,k,1))))*SetPoints(1,order(cycount(stop,k,1)))+SetPoints(2,order(cycount(stop,k,1)))-y0)/(tan(anglesr)-tan(omegas(order(cycount(stop,k,1)))));
            y2 = tan(anglesr)*(x2-x0)+y0;
            points = [c,[x1;y1],[x2;y2]];
            npoly = ConvHull(points);
            pindex = size(npoly,2);
        end
        if(pindex<=2)
            error("A polygon with two vertices was computed... We messed up")
        end
    %In this case, the sector's lines intersect with different edges of the
    %Polygon
    else
        %Check if there is a vertical line
        edge1_pol_check = or(abs(mod(omegas(order(cycount(start,k,1))),2*pi)-pi/2)<=tol,abs(mod(omegas(order(cycount(start,k,1))),2*pi)-3*pi/2)<=tol);
        edge2_pol_check = or(abs(mod(omegas(order(cycount(stop,k,1))),2*pi)-pi/2)<=tol,abs(mod(omegas(order(cycount(stop,k,1))),2*pi)-3*pi/2)<=tol);
        Ar_check = or(abs(mod(anglesr,2*pi)-pi/2)<=tol,abs(mod(anglesr,2*pi)-3*pi/2)<=tol);
        Al_check = or(abs(mod(anglesl,2*pi)-pi/2)<=tol,abs(mod(anglesl,2*pi)-3*pi/2)<=tol);
        check_any_vertical = or(edge1_pol_check,or(edge2_pol_check,or(Ar_check,Al_check)));
        if( check_any_vertical == 1 )
            if(edge1_pol_check==1)
                x1 = SetPoints(1,order(cycount(start,k,1)));
                y1 = tan(anglesl)*(x1-x0)+y0;
                x2 = (tan(anglesr)*x0-tan(omegas(order(cycount(stop,k,1))))*SetPoints(1,order(cycount(stop,k,1)))+SetPoints(2,order(cycount(stop,k,1)))-y0)/(tan(anglesr)-tan(omegas(order(cycount(stop,k,1)))));
                y2 = tan(anglesr)*(x2-x0)+y0;
                points = [c,[x1;y1],[x2;y2],SetPoints(:,order(cycount(start,k,1)))];
                npoly = ConvHull(points);
                pindex = size(npoly,2);
            end
            if(edge2_pol_check==1)
                x1 = (tan(anglesl)*x0-tan(omegas(order(cycount(start,k,1))))*SetPoints(1,order(cycount(start,k,1)))+SetPoints(2,order(cycount(start,k,1)))-y0)/(tan(anglesl)-tan(omegas(order(cycount(start,k,1)))));
                y1 = tan(anglesl)*(x1-x0)+y0;
                x2 = SetPoints(1,order(cycount(stop,k,1)));
                y2 = tan(anglesr)*(x2-x0)+y0;
                points = [c,[x1;y1],[x2;y2],SetPoints(:,order(cycount(start,k,1)))];
                npoly = ConvHull(points);
                pindex = size(npoly,2);
            end
            if(Ar_check==1)
                x1 = (tan(anglesl)*x0-tan(omegas(order(cycount(start,k,1))))*SetPoints(1,order(cycount(start,k,1)))+SetPoints(2,order(cycount(start,k,1)))-y0)/(tan(anglesl)-tan(omegas(order(cycount(start,k,1)))));
                y1 = tan(anglesl)*(x1-x0)+y0;
                x2 = x0;
                y2 = tan(omegas(order(cycount(stop,k,1))))*(x2-SetPoints(1,order(cycount(stop,k,1))))+SetPoints(2,order(cycount(stop,k,1)));
                points = [c,[x1;y1],[x2;y2],SetPoints(:,order(cycount(start,k,1)))];
                npoly = ConvHull(points);
                pindex = size(npoly,2);
            end
            if(Al_check==1)
                x1 = x0;
                y1 = tan(omegas(order(cycount(start,k,1))))*(x1-SetPoints(1,order(cycount(start,k,1))))+SetPoints(2,order(cycount(start,k,1)));
                x2 = (tan(anglesr)*x0-tan(omegas(order(cycount(stop,k,1))))*SetPoints(1,order(cycount(stop,k,1)))+SetPoints(2,order(cycount(stop,k,1)))-y0)/(tan(anglesr)-tan(omegas(order(cycount(stop,k,1)))));
                y2 = tan(anglesr)*(x2-x0)+y0;
                points = [c,[x1;y1],[x2;y2],SetPoints(:,order(cycount(start,k,1)))];
                npoly = ConvHull(points);
                pindex = size(npoly,2);
            end
        else
            x1 = (tan(anglesl)*x0-tan(omegas(order(cycount(start,k,1))))*SetPoints(1,order(cycount(start,k,1)))+SetPoints(2,order(cycount(start,k,1)))-y0)/(tan(anglesl)-tan(omegas(order(cycount(start,k,1)))));
            y1 = tan(anglesl)*(x1-x0)+y0;
            x2 = (tan(anglesr)*x0-tan(omegas(order(cycount(stop,k,1))))*SetPoints(1,order(cycount(stop,k,1)))+SetPoints(2,order(cycount(stop,k,1)))-y0)/(tan(anglesr)-tan(omegas(order(cycount(stop,k,1)))));
            y2 = tan(anglesr)*(x2-x0)+y0;
            points = [c,[x1;y1],[x2;y2],SetPoints(:,order(cycount(start,k,1)))];
            npoly = ConvHull(points);
            pindex = size(npoly,2);
        end
        if(pindex<=2)
            error("A polygon with two vertices was computed... We messed up")
        end
    end

end
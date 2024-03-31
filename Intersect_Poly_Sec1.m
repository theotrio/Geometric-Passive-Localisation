% This functions computes the intersection between a polygon and a sector.
% Input: polygon is a set of the polygon's points, Al,Bl,Ar,Br the lines of 
% the sector.
% Output: The points of the polygon's intersection in counterclockwise order.

function [npoly,pindex] = Intersect_Poly_Sec1(polygon,Al,Ar,c)
    
    %% Initialisation
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
    
    %Compute the gradient lines between the edges of the polygon
    omegas = zeros(1,n);
    for i=1:n
        omegas(i) = atan2((SetPoints(2,i)-SetPoints(2,mod(i,n)+1)),(SetPoints(1,i)-SetPoints(1,mod(i,n)+1)));
    end

	%% Rotational Sweep (clockwise) of the Polygon's Points 
	xy = ones(2,n);
	xy(1,:) = x0*ones(1,n);
	xy(2,:) = y0*ones(1,n);

	points = SetPoints - xy;
    order = atan2(points(2,:),points(1,:));
    
    %{
    %% Keeping the ordering of the vertices' angles consinstent
    % atan2 maps points to [-pi,pi] but if the polygon lies on the left of
    % the centre c then by calculating the angles modulo 2pi we eliminatte 
    % the discontinouity.
    
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
    
        if( abs(phstan - pi/2) <=tol || abs(phstan+pi/2)<=tol)
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
        anglesl = mod(anglesl,2*pi);
        anglesr = mod(anglesr,2*pi);
    end
    %}
    
    if(max(order)-min(order)>=pi)
        order = mod(order,2*pi);
        anglesl = mod(anglesl,2*pi);
        anglesr = mod(anglesr,2*pi);
    end

	%Sort the points according to their angular position
	[angles,order] = sort(order,'descend');

    %% Computation of the Intersection 
    
	%Compute where the sector lies relatively to the polygon
	indexl = 0; % All the points from 1 to indexl are on the left of the sector's left line
	indexr = 0; % All the points from indexr+1 to n are on the right of the sector's right line
    for j=1:n
        if(anglesl< angles(j))
            indexl = indexl+1;
        end
        if(anglesr< angles(j))
            indexr = indexr+1;
        end
    end

	%The sector doesn't intersect with the polygon
    if(indexl==n || indexr==0 || (indexl==0 && indexr==n))
		npoly = polygon;
		pindex = n;
    else
		%Compute the line segments of the polygon 
		betaw = zeros(1,n); %Beta are the gradients of the below segments
		xbw = zeros(1,n);
		ybw = zeros(1,n);
		omegw = zeros(1,n); %Omega are the gradients of the upper segments
		xdw = zeros(1,n);
		ydw = zeros(1,n);

		%Starting from the rightmost vertex compute the upper and lower edges of the polygon
		%To do so we take into consideration that the edges are in counter-clockwise order
		%Compute the edges of the lower half of the polygon 
		betaw(1) = omegas(order(1));
		xbw(1) = SetPoints(1,order(1));
		ybw(1) = SetPoints(2,order(1));
		current_down = mod(order(1),n)+1;

		%Compute the edges of the upper half of the polygon 
		omegw(1) = omegas(cycount(order(1),n,-1));
		xdw(1) = SetPoints(1,order(1));
		ydw(1) = SetPoints(2,order(1));


        for j=2:n-1
			%if the next vertex of the polygon belongs to the lower half of the polygon
			%the edges of the lower half should be updates
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
		%Either one of the edges of the sector is outside the polygon
        if(indexl==0 || indexr==n)
			% The left edge of the sector is outside of the polygon
            if(indexl==0)
				%Compute the intersection of the right edge with the two lines of the polygon
				xsol1 = (tan(anglesr)*x0-tan(omegw(indexr))*xdw(indexr)+ydw(indexr)-y0)/(tan(anglesr)-tan(omegw(indexr)));
				ysol1 = tan(anglesr)*(xsol1-x0)+y0;
				xsol2 = (tan(anglesr)*x0-tan(betaw(indexr))*xbw(indexr)+ybw(indexr)-y0)/(tan(anglesr)-tan(betaw(indexr)));
				ysol2 = tan(anglesr)*(xsol2-x0)+y0;
                if(abs(xsol1-xsol2)>tol && abs(ysol1-ysol2)>tol)
					npoly =  [SetPoints(:,order(1:indexr)) [xsol2 xsol1; ysol2 ysol1] ];
                    npoly = ConvHull(npoly);
					pindex = size(npoly,2);
				else
					npoly =  [SetPoints(:,order(1:indexr)) [xsol1 ; ysol1] ];
                    npoly = ConvHull(npoly);
					pindex = size(npoly,2);
                end
            end
			% The right edge of the sector is outside of the polygon
            if(indexr==n)
				%Compute the intersection of the left edge with the two lines of the polygon
				xsol1 = (tan(anglesl)*x0-tan(omegw(indexl))*xdw(indexl)+ydw(indexl)-y0)/(tan(anglesl)-tan(omegw(indexl)));
				ysol1 = tan(anglesl)*(xsol1-x0)+y0;
				xsol2 = (tan(anglesl)*x0-tan(betaw(indexl))*xbw(indexl)+ybw(indexl)-y0)/(tan(anglesl)-tan(betaw(indexl)));
				ysol2 = tan(anglesl)*(xsol2-x0)+y0;
                if(abs(xsol1-xsol2)>tol && abs(ysol1-ysol2)>tol)
					npoly =  [ [xsol1 xsol2; ysol1 ysol2] SetPoints(:,order(indexl+1:n)) ];
                    npoly = ConvHull(npoly);
					pindex = size(npoly,2);
				else
					npoly =  [ [xsol1 ; ysol1] SetPoints(:,order(indexl+1:n)) ];
                    npoly = ConvHull(npoly);
					pindex = size(npoly,2);
                end
            end
		%Both edges intersect with the polygon	
		else
			xsol1 = (tan(anglesr)*x0-tan(omegw(indexr))*xdw(indexr)+ydw(indexr)-y0)/(tan(anglesr)-tan(omegw(indexr)));
			ysol1 = tan(anglesr)*(xsol1-x0)+y0;
			xsol2 = (tan(anglesr)*x0-tan(betaw(indexr))*xbw(indexr)+ybw(indexr)-y0)/(tan(anglesr)-tan(betaw(indexr)));
			ysol2 = tan(anglesr)*(xsol2-x0)+y0;
            if(abs(xsol1-xsol2)>tol && abs(ysol1-ysol2)>tol)
				zr =  [xsol2 xsol1; ysol2 ysol1];
			else
				zr =  [xsol1 ; ysol1];
            end

			xsol1 = (tan(anglesl)*x0-tan(omegw(indexl))*xdw(indexl)+ydw(indexl)-y0)/(tan(anglesl)-tan(omegw(indexl)));
			ysol1 = tan(anglesl)*(xsol1-x0)+y0;
			xsol2 = (tan(anglesl)*x0-tan(betaw(indexl))*xbw(indexl)+ybw(indexl)-y0)/(tan(anglesl)-tan(betaw(indexl)));
			ysol2 = tan(anglesl)*(xsol2-x0)+y0;
            if(abs(xsol1-xsol2)>tol && abs(ysol1-ysol2)>tol)
				zl =  [xsol1 xsol2; ysol1 ysol2];
			else
				zl =  [xsol1 ; ysol1];
            end
			
			npoly = [zl zr SetPoints(:,order(indexl+1:indexr)) ];
            npoly = ConvHull(npoly);
			pindex = size(npoly,2);
            if(pindex<=2)
                error("A polygon with two vertices was computed... We messed up")
            end
        end

    end
	


end

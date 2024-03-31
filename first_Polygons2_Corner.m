function [pol_index,polygons] = first_Polygons2_Corner(xborder,yborder,Al,Ar,Bl,Br,senpos,rs,su,tol)

k = size(Al,2);
m = size(Al,1);

if(nargin<10)
    tol = 0.1;
end
if(nargin<9)
    su = zeros(k,m);
end

%Initialise the data structure to store Polygons 
%Stores the number of vertices that each polygon has
pol_index = zeros(m,k);
%Polygons are stored in a matrix that has the configuration
% sensor x emitter x 2 (x and y position) x number of polygon's vertices
polygons = zeros(m,k,2,rs);

corners = [xborder(2) xborder(1) xborder(1) xborder(2) ; yborder(2) yborder(2) yborder(1) yborder(1) ];

for i=1:m
    for j=1:k
        
        if(su(i,j)==1)
            cen = [Bl(i,j);Br(i,j)];
            polygon = [cen+[1,0]*tol,cen+[0,-1]*tol,cen+[-1,0]*tol,cen+[0,-1]*tol];
            pol_index = size(polygon,2);
        else

            if(or(Ar(i,j)>pi,Ar(i,j)<-pi))
                Ar(i,j) = mod(Ar(i,j),pi);
            end
            if(or(Al(i,j)>pi,Al(i,j)<-pi))
                Al(i,j) = mod(Al(i,j),pi)-pi;
            end

            %An 1x4 array that points to the edge of the rectangle area. There are
            %four edges North, West, South, and East represented with the numbers
            %1 throuh 4 respectively.
            Al_wall = zeros(1,4);
            Ar_wall = zeros(1,4);

            % A pointer of the left line's quadrant
            Al_quad =0;

            %Find the quadrant of the left border of the sector
            if(Al(i,j)>=pi/2 && Al(i,j)<pi)
                % Set the pointer to indicate that this is the 2nd quadrant
                Al_quad = 2;
                %Find the inteersection with the borders
                y1 = xborder(1)*tan(Al(i,j)) + Bl(i,j);
                if(y1>=yborder(2))
                    x1 = (yborder(2) - Bl(i,j))./tan(Al(i,j));
                    y1 = yborder(2);
                    %Indicate the pointer that Al intersects with wall 1
                    Al_wall(1)=1;
                else
                    x1 = xborder(1);
                    %Indicate the pointer that Al intersects with wall 1
                    Al_wall(2)=1;
                end
            end
            if(Al(i,j)>=0 && Al(i,j)<pi/2)
                Al_quad = 1;
                y1 = xborder(2)*tan(Al(i,j))+Bl(i,j);
                if(y1>=yborder(2))
                    x1 = (yborder(2) - Bl(i,j))./tan(Al(i,j));
                    y1 = yborder(2);
                    Al_wall(1)=1;
                else
                    x1 = xborder(2);
                    Al_wall(4)=1;
                end
            end
            if(Al(i,j)<-pi/2 && Al(i,j)>=-pi)
                Al_quad = 3;
                y1 = xborder(1)*tan(Al(i,j)) + Bl(i,j);
                if(y1<=yborder(1))
                    x1 = (yborder(1) - Bl(i,j))./tan(Al(i,j));
                    y1 = yborder(1);
                    Al_wall(3)=1;
                else  
                    x1 = xborder(1);
                    Al_wall(2)=1;
                end
            end
            if(Al(i,j)<0 && Al(i,j)>=-pi/2)
                Al_quad = 4;
                y1 = xborder(2)*tan(Al(i,j))+Bl(i,j);
                if(y1<=yborder(1))
                    x1 = (yborder(1) - Bl(i,j))./tan(Al(i,j));
                    y1 = yborder(1);
                    Al_wall(3)=1;
                else
                    x1 = xborder(2);
                    Al_wall(4)=1;
                end
            end

            %Find the quadrant of the left border of the sector
            if(Ar(i,j)>=pi/2 && Ar(i,j)<pi)
                %Find the inteersection with the borders
                y2 = xborder(1)*tan(Ar(i,j)) + Br(i,j);
                if(y2>=yborder(2))
                    x2 = (yborder(2) - Br(i,j))./tan(Ar(i,j));
                    y2 = yborder(2);
                    %Indicate the pointer that Ar intersects with wall 1
                    Ar_wall(1) = 1;
                else
                    x2 = xborder(1);
                    %Indicate the pointer that Ar intersects with wall 2
                    Ar_wall(2) = 1;
                end
            end
            if(Ar(i,j)>=0 && Ar(i,j)<pi/2)
                y2 = xborder(2)*tan(Ar(i,j))+Br(i,j);
                if(y2>=yborder(2))
                    x2 = (yborder(2) - Br(i,j))./tan(Ar(i,j));
                    y2 = yborder(2);
                    Ar_wall(1)=1;
                else       
                    x2 = xborder(2);
                    Ar_wall(4)=1;
                end
            end
            if(Ar(i,j)<-pi/2 && Ar(i,j)>=-pi)
                y2 = xborder(1)*tan(Ar(i,j)) + Br(i,j);
                if(y2<=yborder(1))
                    x2 = (yborder(1) - Br(i,j))./tan(Ar(i,j));
                    y2 = yborder(1);
                    Ar_wall(3)=1;
                else
                    x2 = xborder(1);
                    Ar_wall(2)=1;
                end
            end
            if(Ar(i,j)<0 && Ar(i,j)>=-pi/2)
                y2 = xborder(2)*tan(Ar(i,j))+Br(i,j);
                if(y2<=yborder(1))
                    x2 = (yborder(1) - Br(i,j))./tan(Ar(i,j));
                    y2 = yborder(1);
                    Ar_wall(3)=1;
                else
                    x2 = xborder(2);
                    Ar_wall(4)=1;
                end
            end

            %Include the appropriate corner if it is towards one
            if(Al_wall*Ar_wall'==0)
                corner = corners(:,Al_quad);
                p = [ x2 x1; y2 y1];
                points = [senpos(:,i), corner, p];
                pgon = ConvHull(points);
                polygons(i,j,:,1:4) = pgon;
                pol_index(i,j) = 4;
            else
                points = [ x2 x1; y2 y1];
                pol_index(i,j) = 3;
                pgon = [senpos(:,i), points];
                polygons(i,j,:,1:3) = pgon;
            end
        end
    end

end


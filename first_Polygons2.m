function [pol_index,polygons] = first_Polygons2(xborder,yborder,Al,Ar,Bl,Br,senpos,rs)

k = size(Al,2);
m = size(Al,1);

%Initialise the data structure to store Polygons 
%Stors the number of vertices that each polygon has
pol_index = zeros(m,k);
%Polygons are stored in a matrix that has the configuration
% sensor x emitter x 2 (x and y position) x number of polygon's vertices
polygons = zeros(m,k,2,rs);

for i=1:m
    for j=1:k

        if(or(Ar(i,j)>pi,Ar(i,j)<-pi))
            Ar(i,j) = mod(Ar(i,j),pi);
        end
        if(or(Al(i,j)>pi,Al(i,j)<-pi))
            Al(i,j) = mod(Al(i,j),pi)-pi;
        end

        %Find the quadrant of the left border of the sector
        if(Al(i,j)>=pi/2 && Al(i,j)<pi)
            %Find the inteersection with the borders
            y1 = xborder(1)*tan(Al(i,j)) + Bl(i,j);
            if(y1>=yborder(2))
                x1 = (yborder(2) - Bl(i,j))./tan(Al(i,j));
                y1 = yborder(2);
            else
                x1 = xborder(1);
            end
        end
        if(Al(i,j)>=0 && Al(i,j)<pi/2)
            y1 = xborder(2)*tan(Al(i,j))+Bl(i,j);
            if(y1>=yborder(2))
                x1 = (yborder(2) - Bl(i,j))./tan(Al(i,j));
                y1 = yborder(2);
            else
                x1 = xborder(2);
            end
        end
        if(Al(i,j)<-pi/2 && Al(i,j)>=-pi)
            y1 = xborder(1)*tan(Al(i,j)) + Bl(i,j);
            if(y1<=yborder(1))
                x1 = (yborder(1) - Bl(i,j))./tan(Al(i,j));
                y1 = yborder(1);
            else  
                x1 = xborder(1);
            end
        end
        if(Al(i,j)<0 && Al(i,j)>=-pi/2)
            y1 = xborder(2)*tan(Al(i,j))+Bl(i,j);
            if(y1<=yborder(1))
                x1 = (yborder(1) - Bl(i,j))./tan(Al(i,j));
                y1 = yborder(1);
            else
                x1 = xborder(2);
            end
        end

        %Find the quadrant of the left border of the sector
        if(Ar(i,j)>=pi/2 && Ar(i,j)<pi)
            %Find the inteersection with the borders
            y2 = xborder(1)*tan(Ar(i,j)) + Br(i,j);
            if(y2>=yborder(2))
                x2 = (yborder(2) - Br(i,j))./tan(Ar(i,j));
                y2 = yborder(2);      
            else
                x2 = xborder(1);
            end
        end
        if(Ar(i,j)>=0 && Ar(i,j)<pi/2)
            y2 = xborder(2)*tan(Ar(i,j))+Br(i,j);
            if(y2>=yborder(2))
                x2 = (yborder(2) - Br(i,j))./tan(Ar(i,j));
                y2 = yborder(2);
            else       
                x2 = xborder(2);
            end
        end
        if(Ar(i,j)<-pi/2 && Ar(i,j)>=-pi)
            y2 = xborder(1)*tan(Ar(i,j)) + Br(i,j);
            if(y2<=yborder(1))
                x2 = (yborder(1) - Br(i,j))./tan(Ar(i,j));
                y2 = yborder(1);
            else
                x2 = xborder(1);
            end
        end
        if(Ar(i,j)<0 && Ar(i,j)>=-pi/2)
            y2 = xborder(2)*tan(Ar(i,j))+Br(i,j);
            if(y2<=yborder(1))
                x2 = (yborder(1) - Br(i,j))./tan(Ar(i,j));
                y2 = yborder(1);
            else
                x2 = xborder(2);
            end
        end

        points = [ x2 x1; y2 y1];
        pol_index(i,j) = 3;
        p2 = [senpos(:,i),points];
        pgon = ConvHull(p2);
        polygons(i,j,:,1:3) = pgon;
    end

end


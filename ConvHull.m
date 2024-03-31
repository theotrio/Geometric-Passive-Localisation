% This is a function that computes the Convex Hull of points in a 2D plane.
% Input: Set Of Points in a 2D plane
% Output: The convex Hull of these points

function ch = ConvHull(SetPoints)

    n = size(SetPoints,2);

    points = sortrows(SetPoints')';
    
    tol =eps;


    upperhull(:,1:2) =  points(:,1:2);
    k=2;
    
    for i=3:n
        upperhull(:,k+1) = points(:,i);
        k = k+1;
        left_Turn = 1;	
        while(left_Turn==1 && k>2)
        %
            th=atan((upperhull(2,k)-upperhull(2,k-2))/(upperhull(1,k)-upperhull(1,k-2)));
		%
            if(abs(th-pi/2)<=tol || abs(th+pi/2)<=tol)
                upperhull(:,k-1)=[];
            else
                b = upperhull(2,k)-tan(th)*upperhull(1,k);
                if(upperhull(2,k-1)-tan(th)*upperhull(1,k-1)-b<=tol)
				upperhull(:,k-1)=[];
                else
				left_Turn=0;
                end
            end
            k = length(upperhull(1,:));
        end
    end

    lowerhull(:,1:2)=points(:,n:-1:n-1);
    k=2;

    for i=n-2:-1:1
        lowerhull(:,k+1) = points(:,i);
        k =k+ 1;
        left_Turn = 1;
        while(left_Turn==1 && length(lowerhull(1,:))>2)
            th=atan((lowerhull(2,k)-lowerhull(2,k-2))/(lowerhull(1,k)-lowerhull(1,k-2)));
            %
            if(abs(th-pi/2)<=tol || abs(th+pi/2)<=tol )
                lowerhull(:,k-1)=[];
            else
                b = lowerhull(2,k)-tan(th)*lowerhull(1,k);
                if(tan(th)*lowerhull(1,k-1)+b-lowerhull(2,k-1)<=tol)
                    lowerhull(:,k-1)=[];
                else
                    left_Turn=0;
                end
            end
            k = length(lowerhull(1,:));
        end
    end

    k1 = length(lowerhull(1,:));
    k2 = length(upperhull(1,:));

    ch = zeros(2,k1+k2-2);

    ch(:,1:k1) = lowerhull(:,k1:-1:1);
    ch(:,k1+1:k1+k2-2) = upperhull(:,k2-1:-1:2);
end
		


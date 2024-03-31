function y = cycount(value,maximum,rate)

    if (nargin<3)
        rate =1;
    end

    if(rate>0)
        y = mod(value+rate-1,maximum)+1;
    end

    if(rate<0)
        y = mod(value+rate,maximum);
        if( any (y == 0) )
            k = length(y);
            for i= 1:k
                if(y(i)==0)
                    y(i) = maximum(i);
                end
            end
        end
    end
end

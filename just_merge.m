function [y,L,R] = just_merge(vector1,vector2)

    n1 = length(vector1);
    n2 = length(vector2);
    
    %{
    i=1;
    for k=2:n1+n2
        if(abs(y(i-1)-y(i))<=10^-14)
            y(i-1) = [];
        else
            i=i+1;
        end
    end
    %}

    id1 = 1;
    id2 = 1;
    i=1;

    y = zeros(1,n1+n2);
    L = zeros(1,n1+n2-1);
    R = zeros(1,n1+n2-1);
    
    while(id1<=n1 && id2<=n2)
        if(vector1(id1)>=vector2(id2))
            y(i) = vector1(id1);
            id1 = id1+1;
            if(i-1==0)
                R(i) = 1;
            else
                R(i) = mod(R(i-1)+1,n1);
                L(i) = L(i-1);
            end
            i = i+1;
        else
            y(i) = vector2(id2);
            id2 = id2+1;
            if(i-1==0)
                L(i) = 1;
            else
                L(i) = mod(L(i-1)+1,n2);
                R(i) = R(i-1);
            end
            i= i+1;
        end
    end

    if(id1<n1)
        y(i:n1+n2) = vector1(id1:n1);
        for j=i:n1+n2-1
            R(j) = min(L(i-1),R(j-1)+1);
            L(j) = 0;
        end
    else
        y(i:n1+n2) = vector2(id2:n2);
        for j=i:n1+n2-1
            L(j) = min(R(i-1),L(j-1)+1);
            R(j) = 0;
        end
    end

end
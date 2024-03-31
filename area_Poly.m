%Given a 2D polygon as a set of its verticies;
%ordered counter-clockwise this functions returns
%the area of the polygon.

function y=area_Poly(SetVerticies)

n = length(SetVerticies);

sum1=0;
sum2=0;

for i = 1:n
    sum1 = sum1 + SetVerticies(1,i)*SetVerticies(2,mod(i,n)+1);
    sum2 = sum2 + SetVerticies(2,i)*SetVerticies(1,mod(i,n)+1);
end

y = 0.5*(sum1-sum2);


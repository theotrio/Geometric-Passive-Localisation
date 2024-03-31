function y = diameter_Poly(polygon,no_vertices)

y = inf;

for i = 1:no_vertices
   for j =i+1:no_vertices    
       dist = norm(polygon(:,i)-polygon(:,j),2);
       if (dist<=y)
           y = dist;
       end
   end
end

end
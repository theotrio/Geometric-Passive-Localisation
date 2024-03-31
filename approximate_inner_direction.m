function [area,direction] = approximate_inner_direction(thetas,L,R,phi,section_theta,section_phi,xdw,ydw,omegasw,x0,y0)

    n = length(L);
    local_sols = zeros(2,n);

    for i=1:n-1
        %%
        if( L(i)==R(i) && R(i)~=0 && L(i)~=0 )
            
            iL = cycount(L(i),length(xdw),1);
            iR = cycount(R(i),length(xdw),1);
            
            xd = [xdw(L(i)) xdw(R(i))];
            yd = [ydw(L(i)) ydw(R(i))];
            omegs = [omegasw(iL) omegasw(iR)];
            
            area_sol = zeros(1,2);
            dir_sol = zeros(1,2);
            
            area_sol(1) = area_Int_Inside(omegs(1),xd(1),yd(1),x0,y0,phi,thetas(i));
            area_sol(2) = area_Int_Inside(omegs(1),xd(1),yd(1),x0,y0,phi,thetas(i+1));
            dir_sol(1) = thetas(i);
            dir_sol(2) = thetas(i+1);
            [h1,h2] = max(area_sol);
            local_sols(:,i) = [h1;dir_sol(h2)];
        %%    
        elseif(R(i)-L(i)==1 && R(i)~=0 && L(i)~=0)
            
            iL = cycount(L(i),length(xdw),1);
            iR = cycount(R(i),length(xdw),1);

            xd = [xdw(iL) xdw(iR)];
            yd = [ydw(iL) ydw(iR)];
            omegs = [omegasw(iL) omegasw(iR)];
            
            area_sol = zeros(1,4);
            dir_sol = zeros(1,4);
            
            phi_roots = find_roots_darea_inside_LR(omegs,xd,yd,x0,y0,section_theta(R(i)),phi);
            theta_roots = section_theta(R(i))-phi_roots;
            current_theta = section_theta(R(i));
            
            for j=1:2
                if(theta_roots(j)>=thetas(i+1) && theta_roots(j)<=thetas(i)) 
                    area_sol(j) = area_Int_Inside(omegs(1),xd(1),yd(1),x0,y0,theta_roots(j)+phi-current_theta,current_theta) + area_Int_Inside(omegs(2),xd(2),yd(2),x0,y0,current_theta-theta_roots(j),theta_roots(j)); 
                    dir_sol(j) = theta_roots(j);
                else
                    area_sol(j) = -Inf;
                end
            end
            
            area_sol(3) = area_Int_Inside(omegs(1),xd(1),yd(1),x0,y0,thetas(i)+phi-current_theta,current_theta) + area_Int_Inside(omegs(2),xd(2),yd(2),x0,y0,current_theta-thetas(i),thetas(i)); 
            area_sol(4) = area_Int_Inside(omegs(1),xd(1),yd(1),x0,y0,thetas(i+1)+phi-current_theta,current_theta) + area_Int_Inside(omegs(2),xd(2),yd(2),x0,y0,current_theta-thetas(i+1),thetas(i+1)); 
            dir_sol(3) = thetas(i);
            dir_sol(4) = thetas(i+1);
            [h1,h2] = max(area_sol);
            local_sols(:,i) = [h1;dir_sol(h2)];
        
        %%    
        elseif(R(i)-L(i)>=2 && L(i)~=0 && R(i)~=0 )
            
            iL = cycount(L(i),length(xdw),1);
            iR = cycount(R(i),length(xdw),1);
            
            if(iR==1)
                omegs = [omegasw(iL:R(i)),omegasw(1)];
                stheta = [section_theta(iL:R(i)),section_theta(1)];
            else
                omegs = omegasw(iL:iR);
                stheta = section_theta(iL:iR);
            end
            
            xd = xdw(L(i):R(i));
            yd = ydw(L(i):R(i));
            
            phis = section_phi(L(i):R(i));
            
            m = length(xd);
            
            area_sol = zeros(1,4);
            dir_sol = zeros(1,4);
            
            theta_roots = find_roots_darea_Inside_LMR(omegs,xd,yd,x0,y0,stheta);
            
            middle = sum(area_Int_Inside(omegs(2:m-1),xd(2:m-1),yd(2:m-1),x0,y0,phis(2:m-1),stheta(2:m-1)));
            
            for j=1:2
                if(theta_roots(j)>=thetas(i+1) && theta_roots(j)<=thetas(i)) 
                    left_right =  area_Int_Inside(omegs(1),xd(1),yd(1),x0,y0,theta_roots(j)+phi-stheta(1),stheta(1)) + area_Int_Inside(omegs(m),xd(m),yd(m),x0,y0,stheta(m-1)-theta_roots(j),theta_roots(j));
                    area_sol(j) = left_right+middle; 
                    dir_sol(j) = theta_roots(j);
                else
                    area_sol(j) = -Inf;
                end
            end
            
            area_sol(3) = area_Int_Inside(omegs(1),xd(1),yd(1),x0,y0,thetas(i)+phi-stheta(1),stheta(1)) + area_Int_Inside(omegs(m),xd(m),yd(m),x0,y0,stheta(m-1)-thetas(i),thetas(i)) + middle;
            area_sol(4) = area_Int_Inside(omegs(1),xd(1),yd(1),x0,y0,thetas(i+1)+phi-stheta(1),stheta(1)) + area_Int_Inside(omegs(m),xd(m),yd(m),x0,y0,stheta(m-1)-thetas(i+1),thetas(i+1)) + middle;
            dir_sol(3) = thetas(i);
            dir_sol(4) = thetas(i+1);
            [h1,h2] = max(area_sol);
            local_sols(:,i) = [h1;dir_sol(h2)];
             
        end
    
    end
    %%
    [area,index] = max(local_sols(1,:));
    direction = local_sols(2,index);

end
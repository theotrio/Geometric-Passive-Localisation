function [area,direction] = approximate_direction1(thetas,L,R,phi,section_theta,section_phi,xdw,ydw,omegasw,xbw,ybw,betasw,x0,y0,iter_nr)

    if(nargin<15)
        iter_nr = 10^2;
    end

    n = length(L);
    local_sols = zeros(2,n);
    
    acc = 10^-3;
    step = 5;

    for i=1:n-1
        %%
        if( L(i)==R(i) && R(i)~=0 && L(i)~=0 )
            
            xd = [xdw(L(i)) xdw(R(i))];
            yd = [ydw(L(i)) ydw(R(i))];
            omegs = [omegasw(L(i)) omegasw(R(i))];
            xb = [xbw(L(i)) xbw(R(i))];
            yb = [ybw(L(i)) ybw(R(i))];
            betas = [betasw(L(i)) betasw(R(i))];
            
            domain = thetas(i+1):phi:thetas(i);
            if(domain(length(domain))~=thetas(i))
                domain =[domain thetas(i)];
            end
            m = length(domain)+1;
            dir_sol = zeros(1,m);
            area_sol = zeros(1,m);
            
            for j=1:m-2
                
                current_theta = [domain(j)+phi domain(j)];
                
                phir = find_roots_darea(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,domain(j));
                phil = find_roots_darea(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,domain(j)+phi);
                
                nr_dom_unsrt = [phir phil current_theta];
                [nr_dom,ids] = sort(nr_dom_unsrt);
                
                for k=1:length(nr_dom)
                    if(ids(k)==6)
                        start = k;
                    end
                    if(ids(k)==5)
                        stop = k;
                    end
                end
                
                k0=length(nr_dom(start:stop));
                max_area_dom=0;
                d_best = 0;
                area_sol_depth =0;
                dir_sol_depth = 0;
                
                for k=1:k0-1
                    test_point = (nr_dom(start+k-1)+nr_dom(start+k))/2;
                    dal = darea_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point-domain(j),domain(j)+phi);
                    dar = darea_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point-domain(j),domain(j));
                    d2al = d2area_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point-domain(j),domain(j)+phi);
                    d2ar = d2area_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point-domain(j),domain(j));
                    if(dal*dar>=0 && (d2al-d2ar)*(dal-dar)<=0 )
                        initial_point=nr_dom(start+k-1):(nr_dom(start+k)-nr_dom(start+k-1))/step:nr_dom(start+k);
                        dir_sol_depth = zeros(1,length(initial_point));
                        area_sol_depth = zeros(1,length(initial_point));
                        for nri=1:length(initial_point)
                            dir_sol_depth(nri) = Newton_Raphson_Mod("newton_dequation_LMR","newton_d2equation_LMR",initial_point(nri),acc,iter_nr,omegs,xd,yd,betas,xb,yb,x0,y0,phi,current_theta);
                            if(dir_sol_depth(nri)<=domain(j) || dir_sol_depth(nri)>=domain(j+1) )
                                area_sol_depth(nri) = -Inf;
                            else
                                area_sol_depth(nri) = area_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,phi,dir_sol_depth(nri));
                            end 
                        end
                    end
                        
                    if(max_area_dom <= max(area_sol_depth))
                        [max_area_dom,idnri] = max(area_sol_depth);
                        d_best = dir_sol_depth(idnri);
                    end
                    
                end
                
                area_sol(j) = max_area_dom;
                dir_sol(j) = d_best;

            end
            area_sol(m-1) = area_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,phi,thetas(i));
            area_sol(m) = area_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,phi,thetas(i+1));
            dir_sol(m-1) = thetas(i);
            dir_sol(m) = thetas(i+1);
            [h1,h2] = max(area_sol);
            local_sols(:,i) = [h1;dir_sol(h2)];
        %%    
        elseif(R(i)-L(i)==1 && R(i)~=0 && L(i)~=0)
            
            domain = thetas(i+1):phi:thetas(i);
            if(domain(length(domain))~=thetas(i))
                domain =[domain thetas(i)];
            end
            m = length(domain)+1;
            dir_sol = zeros(1,m);
            area_sol = zeros(1,m);
            xd = [xdw(L(i)) xdw(R(i))];
            yd = [ydw(L(i)) ydw(R(i))];
            omegs = [omegasw(L(i)) omegasw(R(i))];
            xb = [xbw(L(i)) xbw(R(i))];
            yb = [ybw(L(i)) ybw(R(i))];
            betas = [betasw(L(i)) betasw(R(i))];
            phis = [section_phi(L(i)) section_phi(R(i))];
            stheta = [section_theta(L(i)+1) section_theta(R(i)+1)];
            
            for j=1:m-2
                current_theta = [domain(j)+phi domain(j)];
                
                phir = find_roots_darea(omegs(2),xd(2),yd(2),betas(2),xb(2),yb(2),x0,y0,domain(j));
                phil = find_roots_darea(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,domain(j)+phi);
                
                nr_dom_unsrt = [phir phil current_theta];
                [nr_dom,ids] = sort(nr_dom_unsrt);
                
                for k=1:length(nr_dom)
                    if(ids(k)==6)
                        start = k;
                    end
                    if(ids(k)==5)
                        stop = k;
                    end
                end
                
                k0=length(nr_dom(start:stop));
                max_area_dom=0;
                d_best = 0;
                area_sol_depth =0;
                dir_sol_depth = 0;
                
                for k=1:k0-1
                    test_point = (nr_dom(start+k-1)+nr_dom(start+k))/2;
                    dal = darea_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point+phi-stheta(1),stheta(1));
                    dar = darea_Int(omegs(2),xd(2),yd(2),betas(2),xb(2),yb(2),x0,y0,test_point-stheta(2),stheta(2));
                    d2al = d2area_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point-domain(j),domain(j)+phi);
                    d2ar = d2area_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point-domain(j),domain(j));
                    if(dal*dar>=0 && (d2al-d2ar)*(dal-dar)<=0 )
                        initial_point=nr_dom(start+k-1):(nr_dom(start+k)-nr_dom(start+k-1))/step:nr_dom(start+k);
                        dir_sol_depth = zeros(1,length(initial_point));
                        area_sol_depth = zeros(1,length(initial_point));
                        for nri=1:length(initial_point)
                            dir_sol_depth(nri) = Newton_Raphson_Mod("newton_dequation_LMR","newton_d2equation_LMR",initial_point(nri),acc,iter_nr,omegs,xd,yd,betas,xb,yb,x0,y0,phi,stheta);
                            if(dir_sol_depth(nri)<=thetas(i+1) || dir_sol_depth(nri)>=thetas(i) )
                                area_sol_depth(nri) = -Inf;
                            else
                                area_sol_depth(nri) = newton_equation_LMR(phis,omegs,xd,yd,betas,xb,yb,x0,y0,phi,stheta,dir_sol_depth(nri));
                            end  
                        end
                    end
                    
                    if(max_area_dom <= max(area_sol_depth))
                        [max_area_dom,idnri] = max(area_sol_depth);
                        d_best = dir_sol_depth(idnri);
                    end
                    
                end
                
                area_sol(j) = max_area_dom;
                dir_sol(j) = d_best;

            end
            
            area_sol(m-1) = newton_equation_LMR(phis,omegs,xd,yd,betas,xb,yb,x0,y0,phi,stheta,thetas(i));
            area_sol(m) = newton_equation_LMR(phis,omegs,xd,yd,betas,xb,yb,x0,y0,phi,stheta,thetas(i+1));
            dir_sol(m-1) = thetas(i);
            dir_sol(m) = thetas(i+1);
            [h1,h2] = max(area_sol);
            local_sols(:,i) = [h1;dir_sol(h2)];
        
        %%    
        elseif(R(i)-L(i)>=2 && L(i)~=0 && R(i)~=0 )
            domain = thetas(i+1):phi:thetas(i);
            m = length(domain);
            if(domain(m)~=thetas(i))
                domain =[domain thetas(i)];
            end
            m = length(domain)+1;
            dir_sol = zeros(1,m);
            area_sol = zeros(1,m);
            xd = xdw(L(i):R(i));
            yd = ydw(L(i):R(i));
            omegs = omegasw(L(i):R(i));
            xb = xbw(L(i):R(i));
            yb = ybw(L(i):R(i));
            betas = betasw(L(i):R(i));
            phis = section_phi(L(i):R(i));
            stheta = section_theta(L(i)+1:R(i)+1);
            
            for j=1:m-2
                current_theta = [domain(j)+phi domain(j)];
                kf = length(stheta);
                
                phir = find_roots_darea(omegs(kf),xd(kf),yd(kf),betas(kf),xb(kf),yb(kf),x0,y0,domain(j));
                phil = find_roots_darea(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,domain(j)+phi);
                
                nr_dom_unsrt = [phir phil current_theta];
                [nr_dom,ids] = sort(nr_dom_unsrt);
                
                for k=1:length(nr_dom)
                    if(ids(k)==6)
                        start = k;
                    end
                    if(ids(k)==5)
                        stop = k;
                    end
                end
                
                k0=length(nr_dom(start:stop));
                max_area_dom=0;
                d_best=0;
                area_sol_depth =0;
                dir_sol_depth = 0;
                
                for k=1:k0-1
                    test_point = (nr_dom(start+k-1)+nr_dom(start+k))/2;
                    dal = darea_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point+phi-stheta(1),stheta(1));
                    dar = darea_Int(omegs(kf),xd(kf),yd(kf),betas(kf),xb(kf),yb(kf),x0,y0,test_point-stheta(kf),stheta(kf));
                    d2al = d2area_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point-domain(j),domain(j)+phi);
                    d2ar = d2area_Int(omegs(1),xd(1),yd(1),betas(1),xb(1),yb(1),x0,y0,test_point-domain(j),domain(j));
                    if(dal*dar>=0 && (d2al-d2ar)*(dal-dar)<=0 )
                        initial_point=nr_dom(start+k-1):(nr_dom(start+k)-nr_dom(start+k-1))/step:nr_dom(start+k);
                        dir_sol_depth = zeros(1,length(initial_point));
                        area_sol_depth = zeros(1,length(initial_point));
                        for nri=1:length(initial_point)
                            dir_sol_depth(nri) = Newton_Raphson_Mod("newton_dequation_LMR","newton_d2equation_LMR",initial_point(nri),acc,iter_nr,omegs,xd,yd,betas,xb,yb,x0,y0,phi,stheta);
                            if(dir_sol_depth(nri)<=thetas(i+1) || dir_sol_depth(nri)>=thetas(i) )
                                area_sol_depth(nri) = -Inf;
                            else
                                area_sol_depth(nri) = newton_equation_LMR(phis,omegs,xd,yd,betas,xb,yb,x0,y0,phi,stheta,dir_sol_depth(nri));
                            end  
                        end
                    end
                    
                    if(max_area_dom <= max(area_sol_depth))
                        [max_area_dom,idnri] = max(area_sol_depth);
                        d_best = dir_sol_depth(idnri);
                    end
                    
                end
                
                area_sol(j) = max_area_dom;
                dir_sol(j) = d_best;
            end
            area_sol(m-1) = newton_equation_LMR(phis,omegs,xd,yd,betas,xb,yb,x0,y0,phi,stheta,thetas(i));
            dir_sol(m-1) = thetas(i);
            area_sol(m) = newton_equation_LMR(phis,omegs,xd,yd,betas,xb,yb,x0,y0,phi,stheta,thetas(i+1));
            dir_sol(m) = thetas(i+1);
            [h1,h2] = max(area_sol);
            local_sols(:,i) = [h1;dir_sol(h2)];
        end
    
    end
    %%
    [area,index] = max(local_sols(1,:));
    direction = local_sols(2,index);

end
% Load the necessary packages for Octave
%pkg load statistics matgeom;
clear;
close;

tic

no_exp =10*26;

run_ratio = floor(no_exp/26);

start = 20231128001;

maxstep = 100;
%The observed area at every step
area_obs = zeros(no_exp,maxstep);

% Set the number of sensors
m = 3;

velocities = zeros(2,no_exp);
velocities(1,:) = 2;
%count = 0;
count = 0;
for i = 1:26
    velocities(2,(i-1)*run_ratio+1:i*run_ratio) = count;
    count = count+0.25;
end
clear i count

ratiov = velocities(2,:)./velocities(1,:);

tol = 35*eps;

diameter_polygons = zeros(3,maxstep,no_exp);
distance_from_centre_storg = zeros(m,maxstep,no_exp);
starting_pos_storg = zeros(2,m,no_exp);
sen_trajectories = zeros(2,m,no_exp,maxstep);

empos_save = zeros(2,4,no_exp,maxstep);

sen_dist = zeros(m*(m-1)/2,maxstep,no_exp);

for i11=start+1:(start+no_exp)

% Initialise the starting positions

if(i11==start+1)
    ndirectory = strcat('Ratio_',num2str(ratiov(i11-start),'%.3f')); 
    mkdir(ndirectory);
    cd(strcat('./',ndirectory));
elseif(ratiov(i11-start)~=ratiov(i11-start-1))
    cd('../');
    ndirectory = strcat('Ratio_',num2str(ratiov(i11-start),'%.3f')); 
    mkdir(ndirectory);
    cd(strcat('./',ndirectory));
end

senpos = 60*rand(2,m)-30;

starting_pos_storg(:,:,i11-start) = senpos;

empos = [15 15 -15 -15; 15 -15 15 -15];

%The maximum velocity of the meitter
v_emitter = velocities(2,(i11-start));

v_sensor = velocities(1,(i11-start));

% Get the number of emitters
k = size(empos,2);

centroid = mean(empos,2);

%The angle of the sector phi
%phi = [pi/12; pi/12; pi/12; pi/12];
phi  = pi/12*ones(m,1);

%Specify the noise generated during transmissions
mnoise = 0.01;

%The number of maximum vertices stored in the data structure
rs = 10;
crs = rs;

print_gif=1;
color = ['b','r','g','m','o'];

step = zeros(2,m);

%Store the previous steps to display the trajectory
previous_steps = zeros(m,2,maxstep+1);

% We define a rectangle area for exploration with
% left-bottom point (xborder(1),yborder(1)) and 
% right-top point (xborder(2),yborder(2))
xborder = 1.5*[-30,30];
yborder = 1.5*[-30,30];
ax = [xborder, yborder];

%---Initialise the set of possible movements that the sensor can make------
no_directions = 8;
step_no = 0:no_directions-1;
steps = (2*pi/no_directions)*step_no;
D = zeros(2,no_directions);
D(1,1:no_directions) = cos(steps);
D(2,1:no_directions) = sin(steps);
clear steps step_no;
%--------------------------------------------------------------------------

%---------------------------  Initial Step  -------------------------------
% Get the sector's left and right border
% after listening to a transmission from the sectors
[Al,Bl,Ar,Br] = listen(senpos,empos,phi,mnoise);

%Initialise the first polygons with the first measurements
[pol_index,polygons]=first_Polygons2(xborder,yborder,Al,Ar,Bl,Br,senpos,rs);

cpolygons = zeros(k,2,rs);
cp_index = zeros(1,k);

for j=1:k
    pol = zeros(2,pol_index(1,j));
    pol(:,1:pol_index(1,j))=polygons(1,j,:,1:pol_index(1,j));
    for seni = 2:m
        [pol,ip] = Intersect_Poly_Sector(pol,Al(seni,j),Ar(seni,j),senpos(:,seni));
    end
    %If the new polygon has more vertices than the storage matrix
    while(ip>crs)
        rpoly = zeros(k,2,2*crs);
        rpoly(:,:,1:crs) = cpolygons; 
        cpolygons = rpoly;
        crs = 2*crs;
    end
    cp_index(j)=ip;
    cpolygons(j,:,1:ip)=pol;
end


%Simulation Presentation

if(print_gif==1)
    %Create a gif
    h = figure;
    %filename = 'G_avg_area.gif';
    filename = strcat(num2str(i11),"_","G_Sym_avg_area_Central_vs",num2str(v_sensor),"_ve",num2str(v_emitter),".gif");
%end

plot(senpos(1,:),senpos(2,:),'o')
axis(ax)
hold("on");
plot(empos(1,:),empos(2,:),'x')
%{
for j=1:m
    for i=1:k
        x = ones(2,3);
        x(:,1:3)=polygons(j,i,:,1:pol_index(j,i));
        pgon = polyshape(x');
        plot(pgon)
    end
end
%}

%end

%if(print_gif==1)
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File 
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    %fnam =strcat(num2str(i11),"_sim1_Central_step_0.png");
    %saveas(gcf,fnam)
end

%pause(1.2)

%--------------------------------------------------------------------------

for i=1:maxstep

    %Use the algorithm to figure out where to move next
    max_area = zeros(k,size(D,2));

    for seni = 1:m
        for j = 1:k
            pol = zeros(2,cp_index(j));
            pol(:,1:cp_index(j))=cpolygons(j,:,1:cp_index(j));
            pol = offset_Poly(pol,v_emitter,tol);
            for j2=1:size(D,2)
                y = find_Maximum_Intersection(pol,senpos(:,seni)+v_sensor*D(:,j2),phi(seni));
                max_area(j,j2) = max(y);
            end
        end
        
        centroid = mean(empos,2);
        distance_from_centre_storg(:,i,i11-start) = vecnorm(senpos - centroid*ones(1,m));

        %Store the sensor's position
        sen_trajectories(:,seni,i11-start,i) = senpos(:,seni);
        
        %
        %Move to the direction that minimises the avg area
        avg_area = sum(max_area,1);
        [~,direction] = min(avg_area);
        step(:,seni) = D(:,direction);
        %}
        
        %{
        %Move to the direction that minimises the max observed area
        [search1,idline] = max(max_area);
        [~,idmax] = max(search1);
        [~,direction] = min(max_area(idline(idmax),:));
        step(:,seni) = D(:,direction);
        %}
        
    end
    
    %Store the previous step
    for seni = 1:m
        previous_steps(seni,:,i) = senpos(:,seni);
    end
    %Move to the calculated position and get a measurement
    senpos = senpos+ v_sensor.*step;
    
    distcount = 1;
    for disti=1:m
        for distj=disti+1:m
            sen_dist(distcount,i,i11-start) = norm(senpos(:,disti)-senpos(:,distj));
            distcount = distcount + 1;
        end
    end
    
    %The emitters move to a random direction
    em_directions = int16(round((no_directions-1)*rand(1,size(empos,2))+1));
    for emi=1:size(empos,2)
        empos(:,emi) = rand(1)*v_emitter*D(:,em_directions(emi))+empos(:,emi);
        empos_save(:,emi,i11-start,i) = empos(:,emi);
    end
     
    
    %Add noise to the current position
    %ran = normrnd(0,0.05);
    %dir_pos = 2*pi*rand(1);
    %pos_noise = ran*[cos(dir_pos),sin(dir_pos)];
    %senpos(:,2) = senpos(:,2)+ pos_noise';
    %clear ran dir_pos pos_noise
    
    [Al,Bl,Ar,Br,meas] = listen(senpos,empos,phi,mnoise);
    %-----------------------------------------------------
    
    %Compute for every emitter the new polygons
    
    for seni=1:m
        for j=1:k
            pol = zeros(2,pol_index(seni,j));
            pol(:,1:pol_index(seni,j))=polygons(seni,j,:,1:pol_index(seni,j));
            pol = offset_Poly(pol,v_emitter,tol);
            [pol,ip] = Intersect_Poly_Sector(pol,Al(seni,j),Ar(seni,j),senpos(:,seni));
            %If the new polygon has more vertices than the storage matrix
            while(ip>rs)
                rpoly = zeros(m,k,2,2*rs);
                rpoly(:,:,:,1:rs) = polygons; 
                polygons = rpoly;
                rs = 2*rs;
            end
            pol_index(seni,j)=ip;
            polygons(seni,j,:,1:ip)=pol;
            %sum_area(seni,i) = sum_area(seni,i) + area_Poly(pol);
        end
    end
    
    
    for j=1:k
        pol1 = zeros(2,cp_index(j));
        pol1(:,1:cp_index(j))=cpolygons(j,:,1:cp_index(j));
        pol1 = offset_Poly(pol1,v_emitter,tol);
        for seni = 1:m
            [pol,ip] = Intersect_Poly_Sector(pol1,Al(seni,j),Ar(seni,j),senpos(:,seni));
            if( is_Outside(pol,empos(:,j)) )
                warning("The emitters position is outside of the polygon") 
            end
            pol1 =pol;
        end
        %If the new polygon has more vertices than the storage matrix
        while(ip>crs)
            rpoly = zeros(k,2,2*crs);
            rpoly(:,:,1:crs) = cpolygons; 
            cpolygons = rpoly;
            crs = 2*crs;
        end
        cp_index(j)=ip;
        cpolygons(j,:,1:ip)=pol;
        area_obs(i11-start,i) = area_obs(i11-start,i) + area_Poly(pol);
        diameter_polygons(j,i,i11-start) = diameter_Poly(pol,ip);
    end
    
    %Simulation Presentation
    if(print_gif==1)
    clf;
    plot(senpos(1,:),senpos(2,:),'o')
    axis(ax)
    hold("on");
    plot(empos(1,:),empos(2,:),'x')
    for seni=1:m
        xpos = zeros(1,i);
        ypos = zeros(1,i);
        xp = previous_steps(seni,1,1:i);
        xpos(1:i) = xp(1:i);
        yp = previous_steps(seni,2,1:i); 
        ypos(1:i) = yp(1:i);
        plot(xpos,ypos,'--gs',...
        'LineWidth',1,...
        'MarkerSize',3,...
        'MarkerEdgeColor',color(seni));
        
    end
    
    for j=1:k
        pol = zeros(2,cp_index(j));
        pol(:,1:cp_index(j))=cpolygons(j,:,1:cp_index(j));
        pgon = polyshape(pol');
        plot(pgon)
    end
    
    end
    if(print_gif==1)
        drawnow
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File 
        imwrite(imind,cm,filename,'gif','WriteMode','append');
        %fnam =strcat(num2str(i11),"_sim1_step_",num2str(i),".png");
        %saveas(gcf,fnam)
    end
    
    
end

close(h);

strcat("The simulation is at ",num2str((i11-start)*100/(no_exp))," %")

end

cd('../');

total_running_time = toc;
strcat("The simulation needed ",num2str(floor(total_running_time/3600)), " hour(s), ", num2str(floor(mod(total_running_time/60,60))), " minutes, and ", num2str(round(mod(total_running_time,60)))," seconds")
save(strcat("experiment_",num2str(start),".mat"),"distance_from_centre_storg","starting_pos_storg","area_obs","ratiov","sen_dist","diameter_polygons","sen_trajectories","empos","empos_save")
%save(strcat("experiment_",num2str(start),"_MM.mat"),"distance_from_centre_storg","starting_pos_storg","area_obs","ratiov","sen_dist")




%ps1 is the positon of the sensors 2xm in 2D
%ps2 is the positon of the emitter 2xk in 2D for k emitters

%Returns a matrix Y that contains the lines of the sector and the measurements of the
% emitters' position. 
%Output y has the following sructure to store lines ax+b
% A1(i,j) = The tangets' amgle of the left lines a of the i-th sector for the j-th emitter
% B1(i,j) = The constant b to translate the j-th left lines, to pass through the i-th sensor
% A2(i,j) = The tangets' angle of the right lines a of the i-th sector for the j-th emitter
% B2(i,j) = The constant b to translate the j-th right lines, to pass through the i-th sensor
% measurements(i,:) is the array that stores the measurement of each emitter as the j-th sensor percieves them.

function [A1,B1,A2,B2,measurements] = listen(ps1,ps2,phi,mnoise)

  %pkg load statistics
  
  if(nargin<4)
      mnoise=1;
  end
  
  %Get the number of the emitters  
  k = size(ps2,2);
  %Get the number of the sensors
  m = size(ps1,2);
  
  if(size(mnoise,2)==1)
      mnoise = mnoise*ones(1,m);
  end
  
  for i=1:m
  if(mnoise(i)*4*(pi/180)>=phi(i)/2)
     mnoise(i)= phi(i)/9;
  end
  end

  %Find the line between the sensros and the emitters

  %Compute the emitters' relative position for x position
  Xs = ps1(1,:)'.*ones(m,k);
  Xe = ps2(1,:).*ones(m,k);
  Xd = Xe-Xs;

  %The same for y position
  Ys = ps1(2,:)'.*ones(m,k);
  Ye = ps2(2,:).*ones(m,k);
  Yd = Ye-Ys;


  %For every sensor; create the measurement using some gaussian white noise 
  sigma = ones(m,k).*mnoise'.*pi/180;
  noise = normrnd(0,sigma);

  % Add noise to the relative position of the emitters
  [Xd,Yd] = turnbyTheta_Matrix(Xd,Yd,noise);

  %measurements = zeros(m,k);
  measurements = atan2(Yd,Xd);

  % The tanget of the left line of the sector
  A1 = measurements+phi/2;
  % The constant to translate the left line to pass through the sensor
  B1 = -tan(A1).*Xs+Ys;
  % The tanget of the right line of the sector
  A2 = measurements-phi/2;
  % The constant to translate the right line to pass through the sensor
  B2 = -tan(A2).*Xs+Ys;

  
end

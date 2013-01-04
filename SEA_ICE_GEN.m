% Generates all the exf stuff I need to simulate sea ice

% According to the packages exf and seaice, I need:
% 10m winds
% relative humidity
% SW and LW radiation
% Evaporation, Precipitation, Runoff
%
% This instantiates a front in the system. 
clear 

out_file = 1;

%% Gridding 
Nt = 24; %Hours. Can do minutes/seconds if need be.

%Grid Dimensions
nx = 120;
ny = 40;
nz = 37;

%Length Scales
Lx = 120e3; dx = Lx/nx;
Ly = 40e3; dy = Ly/ny;
Lz = 800; dz = 5; %Specified for top 200 meters



%Grid %Commented out is a centered grid
x = (1:nx)*dx; %x = x - mean(x);
y = (1:ny)*dy; %y = y - mean(y); 
z = -(1:21)*dz; %z = z - mean(z); 


%FROM BFK: Grid size increases after first 100m of z

for i = 21:nz
    dz = dz*(1.2); 
    z(i) = z(i-1) - dz;
end

[X,Y,Z]=ndgrid(x,y,z);
X0=floor(Lx/2);
Y0=floor(Ly/2);
%% Setting Initial Stratification

prof_z = 0;
p0 = 1000; %kg/m^3
alpha = 2e-4; %Thermal Expansion
T_0 = -1; %Bottom Temperature
delTn = 3;
delTh = 4;
delta = 200;
D = 800;
H = 100;

ML = Z >= -H; %Mixed Layer
%Cold Halocline Layer
DL = Z <= -200; %Below Mixed Layer
CHL = ones(size(Z)) - DL - ML;

TCHL = T_0 + (-100 - Z)/(50) + .5*(delTh/2)*tanh((Y)/(2*Ly)); 
T_1 = 1;

TDL = (T_1/2)*tanh((Z + 500)/400) + .5*(delTh/2)*tanh((Y)/(2*Ly));

%This generates a single front corresponding to a .25 degree increase in
%temperature in the mixed layer across a width of a few grid boxes

INFRONT = X > Lx/2;
BEHIND = X <= Lx/2;

Front = .25*exp(-((X-Lx/2).^2/(400*Lx))).*INFRONT + .25*BEHIND; 

TDL = exp(1).^((Z+H)/delta).*(delTn + (delTh/2)*tanh((Y - mean(y))/Ly));%.*(exp(1).^(-((X-mean(x)).^2/(2*dx))))));% + .01*(rand(size(X))- 1))); %Deep Layer Strat

TML = bsxfun(@plus,TCHL(:,:,6),0*X); %Mixed Layer Strat is taken to be constant from the first below-H deep layer value through the mixed layer

T = TDL.*DL + TML.*ML + TCHL.*CHL;% +Front.*(ML);

%This has corners. Lets smooth it so that we can smoothly move from the CHL
%to the Deep Layer. 
Tinit = T;

for i = 1:nx
    for j = 1:ny
        T(i,j,:) = smooth(squeeze(T(i,j,:)));
    end
end

TFront = T + Front;

%Display the stratification
plot(squeeze(mean(mean(T,2),1)),z)
plot(squeeze(mean(mean(Tinit,2),1)),z,'b.',squeeze(mean(mean(T,2),1)),z,'r-')
legend('Original Data','Smoothed Data Using ''M/smooth''',...
       'Temperature','Vertical')

% Invert for T via linear EOS. 

fprintf(out_file,' creating Arctic_Temp files ... ');
fid = fopen('Arctic_Temp.bin','w','b'); fwrite(fid,T,'real*8'); fclose(fid);
fid = fopen('Front_Arctic_Temp.bin','w','b'); fwrite(fid,TFront,'real*8'); fclose(fid);
fprintf(out_file,'  done.\n');
 
%% Salinity
MLS = Z >= -50;
DLS = Z <= - 200;
CHLS = ones(size(Z)) - MLS - DLS; 

SDL = 34.75 + 0*X;
SML = 32.9 + 0*X;
SCHL = 34.75 + (Z + 200)/(-81);

S = SDL.*DLS + SML.*MLS + SCHL.*CHLS;


fprintf(out_file,' creating Arctic_Salt.bin ... ');
fid = fopen('Arctic_Salt.bin','w','b'); fwrite(fid,S,'real*8'); fclose(fid);
fprintf(out_file,'  done.\n');

%% 10-m winds
w0 = 5; %m/s
Ytop = Y(:,:,1);
ymean = mean(mean(mean(Y)));
%Wind field is zero in y direction, but essentially a gaussian @5 m/s
%centered on the middle in x. 
for i = 1:Nt
    AirU(:,:,i) = w0*cos((pi/2)*((Ytop - ymean)/ymean)) + .1*w0*(rand(size(Ytop))-1).*cos((pi/2)*((Ytop - ymean)/ymean));
    AirV(:,:,i) = 0*X(:,:,i);
end

fprintf(out_file,' creating AirU.bin, AirV.bin ... ');
fid = fopen('AirU.bin','w','b'); fwrite(fid,AirU,'real*8'); fclose(fid);
fid = fopen('AirV.bin','w','b'); fwrite(fid,AirV,'real*8'); fclose(fid);
fprintf(out_file,'   done.\n');

%% 2-m air temperature

Tave = -2; 

for i = 1:7

    AirT(:,:,i) = 0*X(:,:,i) + Tave; %+ (i-1)*40/6; % + (2*rand(size(X(:,:,i))) - 1);
    %Air Temperature has seasonal cycle if you include the third term. We don't
    %want our cooling to come from the Air Temperature.
end
%    AirT(:,:,6) = AirT(:,:,5);
%   AirT(:,:,7) = AirT(:,:,6);
for i = 8:12
    AirT(:,:,i) = AirT(:,:,7);
end
for i = 13:20
    AirT(:,:,i) = AirT(:,:,13-i+8);
end
for i = 21:26
    AirT(:,:,i) = AirT(:,:,20);
end

fprintf(out_file,' creating AirT.bin ... ');
fid = fopen('AirT.bin','w','b'); fwrite(fid,AirT,'real*8'); fclose(fid);
fprintf(out_file,'  done.\n');

%% 2-m specific humidity

AQave = .01
for i = 1:Nt
    AirAQ(:,:,i) = 0*X(:,:,i) + AQave; %.005*(rand(size(X(:,:,i)))-1); 
end


fprintf(out_file,' creating AirAQ.bin ... ');
fid = fopen('AirAQ.bin','w','b'); fwrite(fid,AirAQ,'real*8'); fclose(fid);
fprintf(out_file,'  done.\n');


%% Precipitation

%% Evaporation

%% River and Glacial Runoff

%% LW/SW Radiative fluxes
q0 = 200; %W/m^2
qd = -1834; %W/m^2

for i = 1:Nt
    q(:,:,i) = q0 + 0*X(:,:,1); % + qd*(max(cos(2*pi*i/24),cos(pi*(1/3))) - cos(pi*(1/3))) + 0*X(:,:,1);
end

for i = 1:Nt %Diurnal Cycle
     if q(1,1,i) > 0
     LWflux(:,:,i) = q(:,:,i);
     swflux(:,:,i) = 0 + 0*X(:,:,1);
     else
        LWflux(:,:,i) = 0*X(:,:,1);
        swflux(:,:,i) = q(:,:,i);
     end
end


fprintf(out_file,' creating LWdiurn.bin SWdiurn.bin ... ');
fid=fopen('LWdiurn.bin','w','b');
fid2=fopen('SWdiurn.bin','w','b');

fwrite(fid,LWflux,'real*8');
fwrite(fid2,swflux,'real*8');
fclose(fid);
fclose(fid2);
fprintf(out_file,'  done.\n');

%% Sea Ice Files
% The following four IC files:
% AREA: Fractional Area which is sea ice (0-1)
% Hsnow: Amount of snowcover (meters)
% Heff: Initial thickness in meters
% Hsalt: Salinity of ice in g/m^2
% IceAge: Ice Age

%% Ice Area
%Want to establish that we have a hole in the middle
Midx = Lx/2;
Midy = Ly/2;
Radius = Ly/5;
Xtop = X(:,:,1);
Ytop = Y(:,:,1);
Inhole = (Xtop - Midx).^2 + (Ytop - Midy).^2 < Radius^2;

Area = ones(size(Xtop)) - Inhole;

%% Hsnow
% Zero to begin with
Hsnow = zeros(size(X(:,:,1)));

%% Hsalt
% Wont set

%% Heff
% Random between 0 and 3 meters
Heff = 3*Area;

%% SEA ICE OBCS STUFF
fprintf(out_file, ' creating North and South OBC files ... ');
fid=fopen('Southtemp.OBC','w','b');
fid2=fopen('Northtemp.OBC','w','b');
fid3 = fopen('Southvvel.OBC','w','b');
fid4 = fopen('Northvvel.OBC','w','b');
for t=1:Nt
	Southtemp(:,:,t) = squeeze(T(:,1,:));
    Northtemp(:,:,t) = squeeze(T(:,ny,:));
    Southvvel(:,:,t) = zeros(size(Southtemp(:,:,t)));
    Northvvel(:,:,t) = zeros(size(Northtemp(:,:,t)));
end
%Enforcing the same Temp on the boundaries as were the initial conditions.
fwrite(fid,Southtemp,'real*8');
fwrite(fid2,Northtemp,'real*8');
fwrite(fid3,Southvvel,'real*8');
fwrite(fid4,Northvvel,'real*8');
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fprintf(out_file,'  done.\n');

%% IceAge
% Random age, either a year or fresh (1 day)
%IceAge = 86400*300*(randi(2) - 1) + 86400;


fprintf(out_file,' creating Sea Ice initial conditions ... ');
fid = fopen('Ice_Area.bin','w','b'); fwrite(fid,Area,'real*8'); fclose(fid);
fid = fopen('Ice_Snow.bin','w','b'); fwrite(fid,Hsnow,'real*8'); fclose(fid);
%fid = fopen('Ice_Salt.bin','w','b'); fwrite(fid,Hsalt,'real*8'); fclose(fid);
fid = fopen('Ice_Eff.bin','w','b'); fwrite(fid,Heff,'real*8'); fclose(fid);
%fid = fopen('Ice_Age.bin','w','b'); fwrite(fid,IceAge,'real*8'); fclose(fid);
fprintf(out_file,'  done.\n');




%% Bathymetry
slope = .1;
h0 = 400;
h = X*0 - Lz; %+ h0*(slope*Y)/max(max(Y));
fid = fopen('topog.box','w','b'); fwrite(fid,h,'real*8'); fclose(fid);




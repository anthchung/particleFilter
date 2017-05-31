fileNum = 8;
num = int2str(fileNum);
networkFile = strcat('network',num,'.csv');
gpsFile = strcat('gps',num,'.csv');
dataNetwork = csvread(networkFile, 1);
dataGPS = csvread(gpsFile,1);
networkNo = dataNetwork(:,1);
networkTime = dataNetwork(:,2);
networkLat = dataNetwork(:,3);
networkLong = dataNetwork(:,4);
networkAcc = dataNetwork(:,5);
networkSpeed = dataNetwork(:,6);

gpsNo = dataGPS(:,1);
gpsTime = dataGPS(:,2);
gpsLat = dataGPS(:,3);
gpsLong = dataGPS(:,4);
gpsAcc = dataGPS(:,5);
gpsSpeed = dataGPS(:,6);

%sizeno = size(no);

B615Lat= -33.91783464;
B615Long = 151.2304754;
B616ALat=-33.91798754;
B616ALong=151.231596;


figure
scatter(networkLat, networkLong,'x');
[dx,dy,d2] = diffLatLong(B616ALat,networkLat,B616ALong,networkLong);
%plot(dx,dy,'*');
%legend('Network location', 'GPS location', 'B616A (start)', 'B615 (end)');
%legend('
title('Plot of coordinates');
xlabel('Latitude (°)');
ylabel('Longitude (°)');
hold on
scatter(gpsLat, gpsLong,'*');
%hold on
scatter(B615Lat,B615Long,'o');
scatter(B616ALat,B616ALong,'o');
legend('Network location', 'GPS location', 'B616A (start)', 'B615 (end)');

%new_latitude  = latitude  + (dy / r_earth) * (180 / pi);
%new_longitude = longitude + (dx / r_earth) * (180 / pi) / cos(latitude * pi/180);
%r = 6371e3;
%dy = (B615Lat - B616ALat)*r*pi/180;
%dx = (B615Long - B616ALong)*r*pi/180*cos(B615Lat*pi/180);
%d2 = sqrt(dy^2+dx^2);


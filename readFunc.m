function [networkArray,gpsArray] = readFunc(fileNum)
%reads a network and gps 
%fileNum = 8;
num = int2str(fileNum);
networkFile = strcat('network',num,'.csv');
gpsFile = strcat('gps',num,'.csv');
dataNetwork = csvread(networkFile, 1);
dataGPS = csvread(gpsFile,1);
networkNo = dataNetwork(:,1);
networkTime = dataNetwork(:,2);
if fileNum <=11
    networkTime = (networkTime-networkTime(1))./1e9;
end
if fileNum >=12
    networkTime = (networkTime-networkTime(1));
end
networkLat = dataNetwork(:,3);
networkLong = dataNetwork(:,4);
networkAcc = dataNetwork(:,5);
networkSpeed = dataNetwork(:,6);

networkArray = [networkNo,networkTime,networkLat,networkLong,networkAcc,networkSpeed];

gpsNo = dataGPS(:,1);
gpsTime = dataGPS(:,2);
gpsTime = (gpsTime-gpsTime(1))./1e9;
gpsLat = dataGPS(:,3);
gpsLong = dataGPS(:,4);
gpsAcc = dataGPS(:,5);
gpsSpeed = dataGPS(:,6);

gpsArray = [gpsNo,gpsTime,gpsLat,gpsLong,gpsAcc,gpsSpeed];

end
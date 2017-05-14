function [dx,dy,d] = diffLatLong(lat1,lat2,long1,long2)
%calculates the difference in lat,long and overall distance, usage: [dx,dy,d] = diffLatLong(lat1,lat2,long1,long2)
    r = 6371e3;
    dy = (lat1 - lat2).*r.*pi/180;
    dx = (long1 - long2).*r.*pi/180.*cos(lat1.*pi/180);
    d = sqrt(dy.^2+dx.^2);
end
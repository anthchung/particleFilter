function [network_error, particle_error] = ParticleFilterFunc(fileNum)
%summary
%initalise particles
%add noise to the particles
%move particles according to motion model
%update observation model with particles
%generate probability of the particles based using Gaussian distribution
%with mean = ??? and variance of observations
%normalise the probabilities, dividing by the sum of the probabilities
%resample, selecting variables with higher probabilities
%mean of the particles is the final estimate

%w = warning('query','last');
%id = w.identifier;
id='MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame';
warning('off',id);
close all

B615Lat= -33.91783464;
B615Long = 151.2304754;
B616ALat=-33.91798754;
B616ALong=151.231596;
[x0,y0,d0] = diffLatLong(B616ALat,B616ALat, B616ALong, B616ALong);
[x1,y1,d1] = diffLatLong(B616ALat,B615Lat, B616ALong, B615Long);
[x,y,d] = diffLatLong(B616ALat,B615Lat, B616ALong, B615Long);
theta = atan(y/x);
velocity = 1;
vx = velocity*cos(theta);
vy = velocity*sin(theta);

%initialise variables
t = (0:1:100); %time steps
x = vx*t;
y = vy*t;
var_system = 2; %noise variance in system
var_measure = 100; %noise variance in measurement

%read in the collected data
[networkArray,gpsArray] = readFunc(fileNum);
networkLat = networkArray(:,3);
networkLong = networkArray(:,4);
[x_networkArray,y_networkArray,d_network] = diffLatLong(B616ALat,networkLat,B616ALong,networkLong);
gpsLat = gpsArray(:,3);
gpsLong = gpsArray(:,4);
[x_gps,y_gps,d_gps] = diffLatLong(B616ALat,gpsLat,B616ALong,gpsLong);

%initialise particles
n = 300; %no of particles
prob = 1/n; %initially the particles are equiprobable
x_particles = zeros(1,n);
y_particles = zeros(1,n);
for j = 1:n
    x_particles(j) = x_networkArray(1) + sqrt(var_system)*randn; %randn draws random numbers from a normal/Gaussian distribution
    y_particles(j) = y_networkArray(1) + sqrt(var_system)*randn; 
end
initialAcc = networkArray(1,5);
%particles_moved = particles;
end_time = size(t);
prev_time = 1;
prevNetworkTime = -1;
x_position_out = [sum(prob.*x_particles)];
y_position_out = [sum(prob.*y_particles)];
x_actual_out = [0];
y_actual_out = [0];
x_particleArray = [x_particles];
y_particleArray = [y_particles];
x_network_error = [];
y_network_error = [];




for time = 2:end_time(2)
    %move particles
    for i = 1:n
        x_particles(i) = x_particles(i) + vx*(t(time) - t(prev_time)) + sqrt(var_system)*randn;
        y_particles(i) = y_particles(i) + vy*(t(time) - t(prev_time)) + sqrt(var_system)*randn;
    end
    
    networkTime = find(networkArray(:,2)<t(time));
    [m,o]=size(networkTime);
    lastNetworkTime = networkArray(networkTime(m),2);
    lastNetworkLat = networkArray(networkTime(m),3);
    lastNetworkLong = networkArray(networkTime(m),4);
    [x_network,y_network,d] = diffLatLong(B616ALat,lastNetworkLat, B616ALong, lastNetworkLong);
    
    x_actual = x(time);
    y_actual = y(time);
    accuracy = networkArray(networkTime(m),5);
    if prevNetworkTime ~=lastNetworkTime && initialAcc > 100 && initialAcc > accuracy
        for j = 1:n
            x_particles(j) = x_network(1) + sqrt(var_system)*randn; %randn draws random numbers from a normal dist
            y_particles(j) = y_network(1) + sqrt(var_system)*randn;
        end
    end
    %calculate probability only if there is a new network location
    if prevNetworkTime ~=lastNetworkTime && accuracy < 100
        %create normal distribution 
        prob = (1/sqrt(2*pi*var_measure))*exp(-(sqrt((x_particles-x_network).^2+(y_particles-y_network).^2)).^2/(2*var_measure));
        prob = prob./sum(prob);
        %resample
        probThreshold = 1/(2*n);
        lowProbIndex = find(prob<probThreshold);
        highProb = find(prob>=probThreshold);
        lowSize = size(lowProbIndex);
        highProbIndex = datasample(highProb,lowSize(2));

        %now replace low_prob indexes with newIndex
        x_particles(lowProbIndex) = x_particles(highProbIndex);
        y_particles(lowProbIndex) = y_particles(highProbIndex);
        prob = (1/sqrt(2*pi*var_measure))*exp(-(sqrt((x_particles-x_network).^2+(y_particles-y_network).^2)).^2/(2*var_measure));
        prob = prob./sum(prob);
    else
        prob = prob;
    end

    %assign weights to particle, expected position is p(particle)*particle
    x_position = sum(prob.*x_particles);
    y_position = sum(prob.*y_particles);
    %fprintf('time is: %f, %f, position is %f,%f\n',t(time), lastNetworkTime,x_position,y_position)
    if prevNetworkTime ~=lastNetworkTime
        x_network_error = [x_network_error, x_network - x_actual];
        y_network_error = [y_network_error, y_network - y_actual];
    end
    
    prev_time = time;
    prevNetworkTime = lastNetworkTime;
    %change to preallocate array to speed up, use size(t)
    x_position_out = [x_position_out x_position];
    y_position_out = [y_position_out y_position];
    x_actual_out = [x_actual_out x_actual];
    y_actual_out = [y_actual_out y_actual];
    x_particleArray = [x_particleArray; x_particles];
    y_particleArray = [y_particleArray; y_particles];
    

end

x_error = x_position_out - x_actual_out;
y_error = y_position_out - y_actual_out;
error = (x_error.^2+y_error.^2);
particle_error = sqrt(mean(error));
network_error_array=(x_network_error.^2+y_network_error.^2);
network_error = sqrt(mean(network_error_array));
%create figure
%figure;

%subplot(1,3,1)
%plot(x_position_out, y_position_out, 'X',x_networkArray,y_networkArray,'*', x_gps,y_gps,'o');
%title('Estimated position from starting point');
%legend('Estimated position with particle filter', 'Raw network location', 'GPS location');
%xlabel('x metres from starting point'); ylabel('y metres from starting point');
%{
subplot(1,3,2)
%plot(t,error);
plot(t, x_particleArray);
title('Movement of Particles in the x direction over time');
xlabel('Time (s)');
ylabel('Distance from starting point in the x direction (m)');

subplot(1,3,3)
plot(t, y_particleArray)
title('Movement of Particles in the y direction over time');
ylabel('Distance from starting point in the y direction (m)');
%}
%fprintf('mean error is: %fm\n',mean(error))

%maximise window, code from: http://stackoverflow.com/questions/15286458/automatically-maximize-figure-in-matlab
%drawnow;
%set(get(handle(gcf),'JavaFrame'),'Maximized',1);

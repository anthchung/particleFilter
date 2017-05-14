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
%clear all
close all
%clc

B615Lat= -33.91783464;
B615Long = 151.2304754;
B616ALat=-33.91798754;
B616ALong=151.231596;
[x,y,d] = diffLatLong(B616ALat,B615Lat, B616ALong, B615Long);
theta = atan(y/x);
velocity = 1;
vx = velocity*cos(theta);
vy = velocity*sin(theta);

%initialise variables
%vx = 2; %velocity in x direction
%vy = 2; %velocity in y direction
t = (0:0.5:100); %time steps
%t2 = (100.5:0.5:200);
x = vx*t;
y = vy*t;
var_system = 1; %noise variance in system
var_measure = 1; %noise variance in measurement

%initialise particles
n = 100; %no of particles
x_particles = zeros(1,n);
y_particles = zeros(1,n);
for j = 1:n
    x_particles(j) = x(1) + sqrt(var_system)*randn; %randn draws random numbers from a normal dist
    y_particles(j) = y(1) + sqrt(var_system)*randn;
end
%particles_moved = particles;
end_time = size(t);
prev_time = 1;
x_position_out = [0];
y_position_out = [0];
x_actual_out = [0];
y_actual_out = [0];
x_particleArray = [x_particles];
y_particleArray = [y_particles];

%read in the collected data
[networkArray,gpsArray] = readFunc(8);

for time = 2:end_time(2)
    %move particles
    for i = 1:n
        x_particles(i) = x_particles(i) + vx*(t(time) - t(prev_time)) + sqrt(var_system)*randn;
        y_particles(i) = y_particles(i) + vy*(t(time) - t(prev_time)) + sqrt(var_system)*randn;
    end
    
    networkTime = find(networkArray(:,2)<t(time));
    [m,n]=size(networkTime);
    lastNetworkTime = networkArray(networkTime(m),2);
    x_actual = x(time);
    y_actual = y(time);
    x_observe = x_actual + sqrt(var_measure)*randn;
    y_observe = y_actual + sqrt(var_measure)*randn;
    %calculate probability only if there is a new network location
    if prevNetworkTime ~=lastNetworkTime
        %create normal distribution 
        prob = (1/sqrt(2*pi*var_measure))*exp(-(abs(x_particles - x_observe)+ abs(y_particles - y_observe)).^2/(2*var_measure));
        prob = prob/sum(prob);
    end
    %prob = 1/n;
    %{
    %resample
    probThreshold = 1/(2*n);
    lowProbIndex = find(prob<probThreshold);
    highProb = find(prob>=probThreshold);
    lowSize = size(lowProbIndex);
    highProbIndex = datasample(highProb,lowSize(2));
    
    %now replace low_prob indexes with newIndex
    x_particles(lowProbIndex) = x_particles(highProbIndex);
    y_particles(lowProbIndex) = y_particles(highProbIndex);
    prob = (1/sqrt(2*pi*var_measure))*exp(-(abs(x_particles - x_observe)+ abs(y_particles - y_observe)).^2/(2*var_measure));
    prob = prob/sum(prob);
    %assign weights to particle, expected position is p(particle)*particle
    x_position = sum(prob.*x_particles);
    y_position = sum(prob.*y_particles);
    %}
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
error = sqrt(x_error.^2+y_error.^2);
%{
%create figure
figure;

subplot(1,3,1)
plot(x_position_out, y_position_out, 'X', x_actual_out, y_actual_out)
title('Estimated position from starting point');
legend('Estimated position', 'Actual position');
xlabel('x metres from starting point'); ylabel('y metres from starting point');

subplot(1,3,2)
%plot(t,error);
plot(t, x_particleArray);
title('Movement of Particles in the x direction over time');
xlabel('Metres from starting point in the x direction');

subplot(1,3,3)
plot(t, y_particleArray)
title('Movement of Particles in the y direction over time');
xlabel('Metres from starting point in the y direction');

fprintf('mean error is: %fm\n',mean(error))

%maximise window, code from: http://stackoverflow.com/questions/15286458/automatically-maximize-figure-in-matlab
drawnow;
set(get(handle(gcf),'JavaFrame'),'Maximized',1);
%}
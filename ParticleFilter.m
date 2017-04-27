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

%initialise variables
vx = 2; %velocity in x direction
vy = 2; %velocity in y direction
t = (0:0.5:100); %time steps
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
%x_error = [0];
%y_error = [0];
for time = 2:end_time(2)
    %move particles
    %particles = particles + v*(t(time) - t(prev_time)) + sqrt(var_system)*randn;
    for i = 1:n
        x_particles(i) = x_particles(i) + vx*(t(time) - t(prev_time)) + sqrt(var_system)*randn;
        y_particles(i) = y_particles(i) + vy*(t(time) - t(prev_time)) + sqrt(var_system)*randn;
    end
    x_actual = x(time);
    y_actual = y(time);
    x_observe = x_actual + sqrt(var_measure)*randn;
    y_observe = y_actual + sqrt(var_measure)*randn;
    %create normal distribution 
    %prob = (1/sqrt(2*pi*var_measure))*exp(-(particles - mean(particles)).^2/(2*var_measure));
    %prob = (1/sqrt(2*pi*var_measure))*exp(-(particles - observe).^2/(2*var_measure));
    prob = (1/sqrt(2*pi*var_measure))*exp(-(x_particles - x_observe+ y_particles - y_observe).^2/(2*var_measure));
    prob = prob/sum(prob);
    %prob = 1/n;
   
    %resample
    probThreshold = 1/(2*n);
    lowProbIndex = find(prob<probThreshold);
    highProb = find(prob>=probThreshold*2);
    lowSize = size(lowProbIndex);
    highProbIndex = datasample(highProb,lowSize(2));
    %now replace low_prob indexes with newIndex
    x_particles(lowProbIndex) = x_particles(highProbIndex);
    y_particles(lowProbIndex) = y_particles(highProbIndex);
    prob = (1/sqrt(2*pi*var_measure))*exp(-(x_particles - x_observe+ y_particles - y_observe).^2/(2*var_measure));
    prob = prob/sum(prob);
    %assign weights to particle, expected position is p(particle)*particle
    x_position = sum(prob.*x_particles);
    y_position = sum(prob.*y_particles);
    
    %replace the particles with low prob index with high prob ones
    %for probCounter = 0:size(low_prob);
    %    index = low_prob(probcounter);
    %    newIndex = randsample(high_prob,1);
    %end
    
    prev_time = time;
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
%create figure
figure;

%create plot with positions over time
%subplot(1,2,1)
%plot(t, position_out,'X', t, actual_out)
%title('Estimated and actual positions over time')
%legend('Estimated position', 'Actual position')
%xlabel('Time (s)'); ylabel('Metres from starting position (m)');
subplot(1,2,1)
plot(x_position_out, y_position_out, 'X', x_actual_out, y_actual_out)
title('Estimated position from starting point');
legend('Estimated position', 'Actual position');
xlabel('x metres from starting point'); ylabel('y metres from starting point');

subplot(1,2,2)
plot(t,error);
fprintf('mean error is: %fm\n',mean(error))
%create plot with particles over time
%subplot(1,2,2)
%plot(t, particleArray)
%title('Graph of particles over time')
%xlabel('Time (s)'); ylabel('Metres from starting position (m)');

%maximise window, code from: http://stackoverflow.com/questions/15286458/automatically-maximize-figure-in-matlab
drawnow;
set(get(handle(gcf),'JavaFrame'),'Maximized',1);
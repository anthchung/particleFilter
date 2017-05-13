
%Pr = 1e-7;
%Pt = 50e-3;

%Gt = 1;
%Gr = 1;
%c = 3e8;
%f = 2.4e9;

%Pr = (Pt*Gt*Gr*c^2)/(16*pi^2*d0^2*f^2);

d0 = 5;
d = 10:200; %can specify more d's if needed for more APs
n = 2; %free space

Pt = 18; %dBm, estimated
Pr = -40; %dBm, measured about 5m away from home router
PLd0 = Pt-Pr;

PLd = PLd0 + 10*n*log(d/d0);
Prd = Pt - PLd;
Prd = Prd - 30; %convert to dBW for awgn function
Prdnoise = awgn(Prd, 10);
%plot(Prdnoise)

v = 1;
t = (0:0.001:5);
x = v*t;
%xnoise = awgn(x, 10);
%plot(xnoise)

%initialise particles
n = 10; %no of particles
%upper = 1;
%lower = -1;
%avg = x(1);
%var = 1;
%spread = (-1:(upper-lower)/n:1);
%particle = normpdf(spread, avg, var);
particles = zeros(1,n);
for i = 1:n
    particles(i) = x(1) + randn;
end
%move particles
%particle = particle + v*t(2);

%create normal distribution 
%f = (1/sqrt(2*pi*sqrt(var)))*exp(-(particle-mean(particle)).^2/(2*var)); %unsure of what the mean and var should be

%assign weights to particle, expected position is p(particle)*particle

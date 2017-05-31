fileNum = 16;
[network_error, particle_error_array] = ParticleFilterFunc(fileNum);
for i = 1:100
    [~,particle_error] = ParticleFilterFunc(fileNum);
    [particle_error_array] = [particle_error_array, particle_error];
end

network_error
mean(particle_error_array)
(network_error-mean(particle_error_array))/network_error*100
function [spikelist] = createspikelist(simulation,neuron, IIneuron)
%createspikelist create a single list of spikes for plotting using
%spikeraster
% create spikelist from neuron(i).spikes and IIneuron(i).spikes
% Started LSS 9 Dec 2024
%
% get length of spike raster
n_spikes = 0 ;
for tpno = 1:simulation.N_TPNs
    n_spikes = n_spikes + neuron(tpno).spikecount ;
end
for iino = 1:simulation.N_IIs
    n_spikes = n_spikes + IIneuron(iino).spikecount ;
end
spikelist = zeros([n_spikes 2]) ; % predeclare and then fill
spikeindex = 1 ;
for tpno = 1:simulation.N_TPNs
    spikelist(spikeindex: spikeindex + neuron(tpno).spikecount - 1, 1) = tpno ;
    spikelist(spikeindex: spikeindex + neuron(tpno).spikecount - 1, 2) = neuron(tpno).spikes(1:neuron(tpno).spikecount) ;
    spikeindex = spikeindex + neuron(tpno).spikecount ;
end

for iino = 1:simulation.N_IIs
    spikelist(spikeindex: spikeindex + IIneuron(iino).spikecount - 1, 1) = tpno + iino ;
    spikelist(spikeindex: spikeindex + IIneuron(iino).spikecount - 1, 2) = IIneuron(iino).spikes(1:IIneuron(iino).spikecount) ;
    spikeindex = spikeindex + IIneuron(iino).spikecount ;
end

% make times in seconds, not time steps
spikelist(:, 2) = spikelist(:, 2) * simulation.timestep ;
end % createspikelist
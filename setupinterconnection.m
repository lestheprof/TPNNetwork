function [neuron, IIneuron] = setupinterconnection(simulation, neuron, IIneuron, connectionfile)
%% setupinterconnection: sets up the connections between the neurons
%
% read in table of connections from connectionfile
[n_arcs, nettable] = readnetwork(connectionfile) ;
% format is 1st line is from_ntype	from_nno	to_ntype	to_nno	to_syntype	to_synno	delay
% subsequent lines are
% TPN|II n TPN|II n A|B|AS|BS|S n dt
% tab-separated
%
% Last modified LSS 4 Dec 2024: added delaysamples
%
% each target is a structure with form <neuron_type neuron_number
% synapse_type>
targetno_tp = ones([1 simulation.N_TPNs]) ;
targetno_ii = ones([1 simulation.N_IIs]) ;
for arcs = 1:n_arcs % test each arc
    if strcmp(nettable(arcs,1).from_ntype,'TPN')
        % add this to the neuron targets for TPN neuron
        % nettable(arcs,2).from_nno
        from = nettable(arcs,2).from_nno ;
        tgt.to_ntype = nettable(arcs,3).to_ntype ;
        tgt.to_nno = nettable(arcs,4).to_nno ;
        tgt.to_syntype = nettable(arcs,5).to_syntype ;
        tgt.to_synno = nettable(arcs,6).to_synno ;
        tgt.delay = nettable(arcs,7).delay ;
        tgt.delaysamps = ceil(tgt.delay/simulation.timestep) ; % delay in samples
        neuron(from).targets(targetno_tp(from)) = tgt ;
        targetno_tp(from) = targetno_tp(from) + 1 ;
    end
    if strcmp(nettable(arcs,1).from_ntype,'II')
        % add this to the neuron targets for II neuron
        % nettable(arcs,2).from_nno
        from = nettable(arcs,2).from_nno ;
        tgt.to_ntype = nettable(arcs,3).to_ntype ;
        tgt.to_nno = nettable(arcs,4).to_nno ;
        tgt.to_syntype = nettable(arcs,5).to_syntype ;
        tgt.to_synno = nettable(arcs,6).to_synno ;
        tgt.delay = nettable(arcs,7).delay ;
        tgt.delaysamps = ceil(tgt.delay/simulation.timestep) ; % delay in samples
        IIneuron(from).targets(targetno_ii(from)) = tgt ;
        targetno_ii(from) = targetno_ii(from) + 1 ;
    end
end


end

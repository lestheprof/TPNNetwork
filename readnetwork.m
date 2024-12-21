function [n_arcs, nettable] = readnetwork(filename)
%testreadnetwork reads a network (or set of weights) from the file filename
%  Either:
%% reads network into a table. Ntwork is expected to have 
% from_ntype	from_nno	to_ntype	to_nno	to_syntype	to_synno	delay
% as first line, where 
% (arc source,axon)
% from_ntype is TPN or II, from_nno is number of this neuron
% (arc target, synapse)
% to_ntype is TPN or II
% to_nno is number of target neuron
% to_syntype is synapse type: s for simple (II), a for apical, b for basal,
% as for apical shunt, bs for basal shunt
% to_synno is synapse number of this type
% delay is delay of this arc in seconds

% or
% reads weight table in, returning the number of weights and the weight
% table
% 1st line is neuron_type	neuron_number	syn_type	syn_number	weight
%
% LSS 21 October 2024.
    opts = detectImportOptions(filename);
    opts.DataLines = 2; % start data at line 2
    nettable = readtable(filename, opts) ;
    n_arcs = height(nettable) ; % number of arcs
end
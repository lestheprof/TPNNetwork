function [apicalactivation, basalactivation,ahactiv,  spikelist] = Network(varargin)
% Network: encapsulates a whole network.
% Includes TPNs and II neurons
% 
% All parameters are  <name,value> pairs using varargin
%% simulation parameters (in network.xxx)
% duration: time in seconds for running the neuron
% timestep:  timestep length in seconds
% N_TPNs: no of TPNs
% N_IIs: number of inhibitory interneurons
% Needs definition of external inputs
% Needs definition of internalconnectivity
%
%% TPN generic values (in TPNs.xxx)
% tau_apical: time constant for apical synapses (alpha function)
% tau_basal: time constant for basal synapses (alpha function)

%% parameters describing synapses for TPNs (in TPNs.xxx, also as an arrau,1 per TPN)
% concept is that dendrite has a local capacitance and membrane
% resistance. In addition there is a resistance describing the maximal
% conductance (=minimum resistance) of the dendrite that governs the amount
% of charge transferred. Keep one set of values for apical and one for
% basal.
% assume these are same or all TPNs
% c_apical: capacitance at apical dendrite
% r_apical: apical leakage resistance
% c_basal: capacitance at basal dendrite
% r_basal: basal leakage resistance
% r_synap_dendrite: resistance (1/conductance) of apical dendrite.
% r_synba_dendrite: resistance of basal dendrite:
% apicalshuntduration: time that apical shunt is active for (in seconds)
% basalshuntduration: time that basal shunt is active for (in seconds)

%% values for each TPN (also in TPNs.xxx, but as an array, 1 per TPN.)
% noapicalinputs:  number of apical inputs for each TPN
% nobasalinputs: number of basal inputs for each TPN
% noapicalshunts: number of apical shunting synapses for each TPN
% nobasalshunts: number of basal shunting synapses for each TPN

%% actual weights: 1-d array for each TPN
% apicalsynapseweights: weights on apical synapses, must be same length as noapicalinputs
% basalsynapsewights: weights on basal synapses, must be same length as nobasalinputs

%% parameters for shunting synapses: 1-d array for each TPN
% apicalshuntweights: weights for apical shunt: 0 leq value leq 1
% basalshuntweights: weights for basal shunt: 0 leq value leq 1

%% parameters for spike generation: 1 poer TPN
% thresh_value: initial threshold for dtetermining spike (applied to ahactiv
% value)
% refractoryperiod: refractory period during which no spikes will be
% generated.
% relrefperiod: relative refractory period after which threshold returns to
%  initial value
% thresh_leap : leap in threshold value after Refractory period
% thresh_decay : exponent of decay (exp(-thresh_decay * time))
%% parameters for returning spike list
% neuronid: identity of this neuron (integer)
% maxnospikes: maximum number of spikes this neuron can return

%% omit (these are either externalor internal)
% apicalinputs: actual apical spike inputs in <time synapse_number> format
% basalinputs: same format
% apicalshuntinputs: actual apical shunt inputs, in <time synapse_number> format
% basalshuntinputs: actual basal shunt inputs, in <time synapse_number> format

% returns apicalactivation, basalactivation
%  axon hillock activity (ahactiv) level and spike times array <neuronid time>
%
% TPN.m started 16 July 2024 by LSS, this new function is based on
% the version 12 8 24 LSS
% TPN_spines started 16 August, LSS.
% TPN_spines_shun started 19 August 2024 LSS
% added apical and basal shunting synapses (working) 3 Sept 2024

% get the parameters using  varargin name value system
%
% default neuronid
neuronid = 1 ;
% default max number of spikes
maxnospikes = 1000 ;
% default timestep
timestep= 0.0001; % 100 microseconds
% lest there be no shunting information, ensure program doesn't crash
noapicalshunts = 0 ;
nobasalshunts = 0 ;

%
i=1 ;
while(i<=size(varargin,2))
    switch lower(varargin{i})
        case 'duration' % simulation
            network.duration=varargin{i+1};
            i=i+1;
        case 'timestep' % simulation
            network.timestep=varargin{i+1};
            i=i+1;
        case 'noapicalinputs' % TPN no of apical inputs
            noapicalinputs =varargin{i+1};
            i=i+1;
        case 'nobasalinputs' % TPN no of basal inputs
            nobasalinputs=varargin{i+1};
            i=i+1;
        case 'apicalinputs' % (this simulation only) actual apical spike inputs
            apicalinputs = varargin{i+1};
            i=i+1;
        case 'basalinputs' % (this simulation only) actual basal spike inputs
            basalinputs = varargin{i+1};
            i=i+1;
        case 'apicalsynapseweights' % TPN weights on apical synapses
            % Must be noapicalinputs long
            apicalsynapseweights = varargin{i+1};
            if (length(apicalsynapseweights) ~= noapicalinputs)
                error("Apical synapse weight vector not same as number of apical inputs") ;
            end
            i=i+1;
        case 'basalsynapsewights' % TPN weights on basal synapses
            % Must be nobasalinputs long
            basalsynapseweights = varargin{i+1};
            if (length(basalsynapseweights) ~= nobasalinputs)
                error("Basal synapse weight vector not same as number of basal inputs") ;
            end
            i=i+1;
        case 'tau_apical' % all TPNs used in alpha function for apical synapses
            TPNs.tau_apical = varargin{i+1};
            i=i+1;
        case 'tau_basal' % all TPNs used in alpha function for basal synapses
            TPNs.tau_basal = varargin{i+1};
            i=i+1;
        case 'c_apical' %  all TPNs capacitance of apical dendrite (lumped: used after summing currents)
            c_apical = varargin{i+1};
            i=i+1;
        case 'c_basal' %  all TPNs capacitance of basal dendrite (lumped: used after summing currents)
            c_basal = varargin{i+1};
            i=i+1;
        case 'r_apical' % all TPNs resistance across c_apical: inverse of leak (lumped: used after summing currents)
            r_apical = varargin{i+1};
            i=i+1;
        case 'r_basal' %  all TPNs resisitance across c_basal: inverse of leak (lumped: used after summing currents)
            r_basal = varargin{i+1};
            i=i+1;
        case 'r_synap_dendrite' % all TPNs  resistance of apical synapses
            r_synap_dendrite = varargin{i+1};
            i=i+1;
        case 'r_synba_dendrite' %  all TPNs resistance of basal synapses
            r_synba_dendrite = varargin{i+1};
            i=i+1;
        case 'c_apical_spine' % all TPNs  capacitance of apical dendrite (lumped: used after summing currents)
            c_apical_spine = varargin{i+1};
            i=i+1;
        case 'c_basal_spine' %  all TPNs capacitance of basal dendrite (lumped: used after summing currents)
            c_basal_spine = varargin{i+1};
            i=i+1;
        case 'r_apical_spine' % all TPNs  resistance across c_apical: inverse of leak (lumped: used after summing currents)
            r_apical_spine = varargin{i+1};
            i=i+1;
        case 'r_basal_spine' % all TPNs  resisitance across c_basal: inverse of leak (lumped: used after summing currents)
            r_basal_spine = varargin{i+1};
            i=i+1;
        case 'r_synap_spine' %  all TPNs resistance of apical synapses
            r_synap_spine = varargin{i+1};
            i=i+1;
        case 'r_synba_spine' % all TPNs  resistance of basal synapses
            r_synba_spine = varargin{i+1};
            i=i+1;
        case 'thresh_value' %  all TPNs initial threshold for spiking
            thresh_value = varargin{i+1};
            i=i+1;
        case 'thresh_leap' %  all TPNs amount by which threshold jumps
            thresh_leap = varargin{i+1};
            i=i+1;
        case 'thresh_decay' %  all TPNs decay rate for threshold
            thresh_decay = varargin{i+1};
            i=i+1;
        case 'refractoryperiod' %  all TPNs period during which no spikes can be generated after a spike
            refractoryperiod = varargin{i+1};
            i=i+1;
        case 'relrefperiod' %  all TPNs relative refractory period
            relrefperiod = varargin{i+1};
            i=i+1;
        case 'neuronid' % Neuron  identity of neuron: integer
            neuronid = varargin{i+1};
            i=i+1;
        case 'maxnospikes' % Neuron max number of spikes for this neuron
            maxnospikes = varargin{i+1};
            i=i+1;
        case 'noapicalshunts' % TPN number of apical shunting synapses
            noapicalshunts = varargin{i+1};
            i=i+1;
        case 'nobasalshunts' % TPN number of basal shunting synapses
            nobasalshunts = varargin{i+1};
            i=i+1;
        case 'apicalshuntinputs' % (this simulation only) actual apical shunt inputs, in <time synapse_number> format'
            apicalshuntinputs = varargin{i+1};
            i=i+1;
        case 'basalshuntinputs' % (this simulation only) actual basal shunt inputs, in <time synapse_number> format'
            basalshuntinputs = varargin{i+1};
            i=i+1;
        case 'apicalshuntweights' % TPN weights for apical shunt: 0 leq value leq 1'
            apicalshuntweights = varargin{i+1};
            if (length(apicalshuntweights) ~= noapicalshunts)
                error("Apical shunting synapse weight vector not same as number of apical shunts") ;
            end
            i=i+1;
        case 'basalshuntweights' % TPN weights for basal shunt: 0 leq value leq 1'
            basalshuntweights = varargin{i+1};
            if (length(basalshuntweights) ~= nobasalshunts)
                error("Basal shunting synapse weight vector not same as number of basal shunts") ;
            end
            i=i+1;
        case 'apicalshuntduration' % all TPNs  ime that apical shunt is active for (in seconds)'
            apicalshuntduration = varargin{i+1};
            i=i+1;
        case 'basalshuntduration' % All TPNs time that basal shunt is active for (in seconds)
            basalshuntduration = varargin{i+1};
            i=i+1;
        otherwise
            error('TPN.m: Unknown argument %s given',varargin{i});
    end
    i=i+1;
end

% length in timesteps of simulation
simlength = ceil(duration/timestep) ;

% calculate alpha functions for apical and basal dendrites (tau need not be
% the same: we might later want a variety of these for excitatory and inhibitory synapses)
alpha_apical = setupAlphaFunction(timestep, tau_apical) ;
alpha_basal = setupAlphaFunction(timestep, tau_basal );

% array for apical synapse activations (extra length to stop late spikes
% falling off the end)
% the per-synapse post-tynaptic values are used in this version
apicalpostsynapse = zeros([noapicalinputs simlength  + length(alpha_apical)]);
% array for basal activations
basalpostsynapse = zeros([nobasalinputs simlength + length(alpha_basal)]);
% vector for apical current
apicalcurrent = zeros([1 simlength + length(alpha_apical)]) ;
% vector for basal current
basalcurrent = zeros([1 simlength + length(alpha_basal)]) ;
% vector for apical activation (charged by apicalcurrent, capacitor and
% parallel resistor)
apicalactivation = zeros([1 simlength + length(alpha_apical)]) ;
% vector for basal activation (charged by basalcurrent, capacitor and
% parallel resistor)
basalactivation = zeros([1 simlength + length(alpha_basal)]) ;
% vector for axon hillock activation
ahactiv = zeros([1 simlength]) ;
% vector for threshold: useful for seeing what has been going on
% threshold = ones([1 simlength]) * thresh_value ; NO: cannot declare till
% calc-thresh_increment has been calculated
% vector for recording spikes
% spikevector = zeros([1 simlength]) ; no longer required


% sort apical and basal inputs, additive and shunting, into time order
% both apical and basal inputs are a 2d array with
% N_spikes rows 2 columns <time neuron>
if  (exist('apicalinputs','var') && ~isempty(apicalinputs) )% allow empty apical input
    apicalinputs = sortrows(apicalinputs, 1) ;
    % replace times with timestep numbers
    apicalinputs(:,1) = round(apicalinputs(:,1)/timestep);
    % check that the apical spike inputs are within range
    if (max(apicalinputs(:,2)) > noapicalinputs)
        error("TPN.m: apical input out of range");
    end
end
if  (exist('basalinputs','var') && ~isempty(basalinputs)) % allow empty basal input
    basalinputs = sortrows(basalinputs, 1) ;
    % replace times with timestep numbers
    basalinputs(:,1) = round(basalinputs(:,1)/timestep);
    % check that the basal spike inputs are within range
    if (max(basalinputs(:,2)) > nobasalinputs)
        error("TPN.m: basal input out of range");
    end
end
if (exist('apicalshuntinputs', 'var') && ~isempty(apicalshuntinputs))
    apicalshuntinputs = sortrows(apicalshuntinputs, 1);
    % replace times with timestep numbers
    apicalshuntinputs(:,1) = round(apicalshuntinputs(:,1)/timestep);
    % check that the basal spike inputs are within range
    if (max(apicalshuntinputs(:,2)) > noapicalshunts)
        error("TPN.m: apical shunt input out of range");
    end
end
if (exist('basalshuntinputs', 'var') && ~isempty(basalshuntinputs))
    basalshuntinputs = sortrows(basalshuntinputs, 1);
    % replace times with timestep numbers
    basalshuntinputs(:,1) = round(basalshuntinputs(:,1)/timestep);
    % check that the basal spike inputs are within range
    if (max(basalshuntinputs(:,2)) > nobasalshunts)
        error("TPN.m: basal shunt input out of range");
    end
end

% set up shunting synapses but only if there are shunting synapses
if (noapicalshunts > 0)
    asduration = floor(apicalshuntduration/timestep) ; % get apical shunt duration in timesteps
    if (asduration < 1)
        asduration = 1 ;
    end % 0 is possible, but needs reset to 1
    % calculate shunt value per timestep for each shunting synapse weight
    ap_sh_wt_pt = zeros([1 length(apicalshuntweights)]) ; % preallocate
    for sno = 1: length(apicalshuntweights)
        ap_sh_wt_pt(sno) = (1-apicalshuntweights(sno))^(1/asduration); % weight to apply per timestep
    end
end
if (nobasalshunts  > 0)
    bsduration = floor(basalshuntduration/timestep) ; % get basal shunt duration in timesteps
    if (bsduration < 1)
        bsduration = 1 ;
    end % 0 is possible, but needs reset to 1
    ba_sh_wt_pt = zeros([1 length(basalshuntweights)]) ; % preallocate
    for sno = 1: length(basalshuntweights)
        ba_sh_wt_pt(sno) = (1 - basalshuntweights(sno))^(1/bsduration); % weight to apply per timestep
    end
end


% calculate apical leakiness per timestep:how much leaks away from
% apicalactivation voltage in 1 timestep
% C * dV/dt = I = V/R => dV/dt = V/(R*C)
% so fraction that leaks away is timestep/(R*C)
apicalfracleak = timestep / (r_apical * c_apical) ;
basalfracleak = timestep / (r_basal * c_basal) ;
apicalspinefracleak = timestep/(r_apical_spine  * c_apical_spine) ;
basalspinefracleak = timestep/(r_basal_spine  * c_basal_spine) ;


% calculate the amount tio be added to the threshold whne a spike occurs.
thresh_increment = calc_thresh_increment(thresh_leap, thresh_decay, ...
    refractoryperiod, relrefperiod, timestep) ;
th_inc_length = length(thresh_increment) ;
% vector for threshold: useful for seeing what has been going on
threshold = ones([1 (simlength + th_inc_length)]) * thresh_value ;



% simulate by running step by step
apicalspikeno = 1;  % where are we in the list of apical spikes
basalspikeno = 1;  % where are we in the list of basal spikes
apicalshuntspikeno = 1; % where we are in list of apical shunting spikes
if (exist('apicalshuntinputs', 'var') && ~isempty(apicalshuntinputs))
    inapicaltimeinterval = zeros([1 size(apicalshuntinputs, 1)]) ;
    countdown_ap = zeros([1 size(apicalshuntinputs, 1)]) ;
end
if (exist('basalshuntinputs', 'var') && ~isempty(basalshuntinputs))
    basalshuntspikeno = 1 ; % where we are in list of basal shunting spikes
    inbasaltimeinterval = zeros([1 size(basalshuntinputs, 1)]) ;
    countdown_bs = zeros([1 size(basalshuntinputs, 1)]) ;
end
no_of_spikes = 0 ;
spikelist = zeros([maxnospikes 2]) ;
spikelist(: , 1) = neuronid ;
for tstep = 1:simlength
    % has a spike been generated? If so, store it
    spiked = runstep(tstep) ;
    if (spiked)
        no_of_spikes = no_of_spikes + 1 ;
        if (no_of_spikes == maxnospikes)
            error('TPN.m: max number of spikes for neuron %i exceeded.', neuronid);
        end
        spikelist(no_of_spikes, 2) = tstep * timestep ;
    end

end

spikelist = spikelist(1:no_of_spikes, :) ; % retrun only the used part


    function isspike = runstep(ts)
        % run simulation for one timestep, returning true if a spike is
        % generated
        % calculate apical synaptic depolarisation: fix, there can be > 1
        % spike in a single timestep so jse while
        while ((apicalspikeno <= size(apicalinputs,1)) && ...
                (apicalinputs(apicalspikeno,1) == ts)) % calculate for each synapse: is there a new apical spike input?
            % apical spike detected: add alpha function * weight
            % calculate for each synapse used in this version
            apicalpostsynapse(apicalinputs(apicalspikeno,2), ts:(ts + length(alpha_apical) -1)) = ...
                (apicalpostsynapse(apicalinputs(apicalspikeno,2), ts:(ts + length(alpha_apical) -1)) * (1-apicalspinefracleak)) + ...
                (alpha_apical * apicalsynapseweights(apicalinputs(apicalspikeno,2)) / (r_synap_spine + r_synap_dendrite) );
            % adjust using spine capacitance leakage and resistance

            % add this change to the total depolarisation (currently no
            % geometry on apical dendrite)
            apicalcurrent(ts:(ts + length(alpha_apical) -1)) = ...
                apicalcurrent(ts:(ts + length(alpha_apical) -1)) + ...
                apicalpostsynapse(apicalinputs(apicalspikeno,2), ts:(ts + length(alpha_apical) -1) ) ;
            % (alpha_apical * apicalsynapseweights(apicalinputs(apicalspikeno,2)) / r_synap_dendrite )  ;
            apicalspikeno = apicalspikeno + 1 ;
        end

        % calculate apical activation:
        % total charge is weight coulombs!!! That's why the voltages are
        % so high: fixed.
        % is the apical current a current (I) or a charge (I * Delta T = Q)?
        if (ts > 1)
            apicalactivation(ts) = apicalactivation(ts-1) * (1 - apicalfracleak) + ...
                (apicalcurrent(ts)/c_apical) ;
            % ((timestep * apicalcurrent(ts))/c_apical) ;
            % (timestep * (apicalcurrent(ts)/c_apical)) ;
        end
        % apply apical shunt, if any.
        if ((noapicalshunts > 0) && (exist('apicalshuntinputs','var')))
            %   apply shunting synaptic inhibition here.
            %   for each apical shunting synapse, add up the inhibition and then apply.
            % is there a shunting synapse at this time, or has there been one
            % within the last asduration timesteps?
            while (apicalshuntspikeno <= size(apicalshuntinputs, 1) && ...
                    (apicalshuntinputs(apicalshuntspikeno,1) == ts))
                if (apicalshuntinputs(apicalshuntspikeno,1) == ts) % if there's a new shunting input, initialise countdown
                    inapicaltimeinterval(apicalshuntspikeno)  = true;
                    % keep the id of the spike synapse number, apicalshuntinputs(apicalshuntspikeno,2)
                    % and use it when indexing ap_sh_wt_pt
                    countdown_ap(apicalshuntspikeno) = asduration ;
                end
                apicalshuntspikeno = apicalshuntspikeno + 1;
            end % while
            % calculate total shunting effect
            shunteffect = 1 ;
            for apsspikeindex = 1:length(inapicaltimeinterval)
                if (inapicaltimeinterval(apsspikeindex)  == true)
                     shunteffect = shunteffect * ap_sh_wt_pt(apicalshuntinputs(apicalshuntspikeno - 1,2)) ;
                    % shunteffect = shunteffect * ap_sh_wt_pt(apsspikeindex) ; % issue: ap_sh_wt_pt see above.
                    countdown_ap(apsspikeindex) = countdown_ap(apsspikeindex) - 1;
                    if (countdown_ap(apsspikeindex) == 0)
                        inapicaltimeinterval(apsspikeindex) = false ;
                    end % if
                end% if
            end  % for
            % apply shunteffect
            apicalactivation(ts) = apicalactivation(ts) * shunteffect ;
        end % if



        % calculate basal synaptic depolarisation
        while ((basalspikeno <= size(basalinputs,1)) ...
                && (basalinputs(basalspikeno,1) == ts))
            % basal spike detected: add alpha function * weight
            % calculate for ecach synapse (not used in this version)
            basalpostsynapse(basalinputs(basalspikeno,2), ts:(ts + length(alpha_basal) -1)) = ...
                (basalpostsynapse(basalinputs(basalspikeno,2), ts:(ts + length(alpha_basal) -1)) *  (1-basalspinefracleak)) + ...
                (alpha_basal * basalsynapseweights(basalinputs(basalspikeno,2)) / (r_synba_spine + r_synba_dendrite)) ;
            % add this change to the total depolarisation
            basalcurrent(ts:(ts + length(alpha_basal) -1)) = ...
                basalcurrent(ts:(ts + length(alpha_basal) -1)) + ...
                basalpostsynapse(basalinputs(basalspikeno,2), ts:(ts + length(alpha_basal) -1)) ;
            % (alpha_basal * basalsynapseweights(basalinputs(basalspikeno,2))  / r_synba_dendrite) ;
            basalspikeno = basalspikeno + 1 ;
        end

        % calculate basal activation:
        % is the basal current a current (I) or a charge (I * Delta T = Q)?
        if (ts > 1)
            basalactivation(ts) = basalactivation(ts-1) * (1 - basalfracleak) + ...
                (basalcurrent(ts)/c_basal) ;
            % ((timestep * basalcurrent(ts))/c_basal) ;
            % (timestep * (basalcurrent(timestep)/c_basal)) ;
        end
        if ((nobasalshunts  > 0)  && (exist('basalshuntinputs', 'var')))
            %   apply shunting synaptic inhibition here.
            %   for each basal shunting synapse, add up the inhibition and then apply            
            %   apply shunting synaptic inhibition here.
            % is there a shunting synapse at this time, or has there been one
            % within the last asduration timesteps?
            while (basalshuntspikeno <= size(basalshuntinputs, 1) && ...
                    (basalshuntinputs(basalshuntspikeno,1) == ts))
                if (basalshuntinputs(basalshuntspikeno,1) == ts) % if there's a new shunting input, initialise countdown
                    inbasaltimeinterval(basalshuntspikeno)  = true;
                    % keep the id of the spike synapse number, basalshuntinputs(apicalshuntspikeno,2)
                    % and use it when indexing ap_sh_wt_pt
                    countdown_bs(basalshuntspikeno) = bsduration ;
                end
                basalshuntspikeno = basalshuntspikeno + 1;
            end % while
            % calculate total shunting effect
            shunteffect = 1 ;
            for bssspikeindex = 1:length(inbasaltimeinterval)
                if (inbasaltimeinterval(bssspikeindex)  == true)
                     shunteffect = shunteffect * ba_sh_wt_pt(basalshuntinputs(basalshuntspikeno - 1,2)) ;
                    countdown_bs(bssspikeindex) = countdown_bs(bssspikeindex) - 1;
                    if (countdown_bs(bssspikeindex) == 0)
                        inbasaltimeinterval(bssspikeindex) = false ;
                    end % if
                end% if
            end  % for
            % apply shunteffect
            basalactivation(ts) = basalactivation(ts) * shunteffect ;


        end

        % Try equation from Reza & Ahsan (R^2 + 2*R +2*C *
        % (1+mod(R)))for axon hillock activation
        ahactiv(ts) = basalactivation(ts)^2  +  2 *  basalactivation(ts) +...
            2 * apicalactivation(ts) * (1 + abs(basalactivation(ts))) ;

        % decide whether to spike at this timestep, and store spike if so.
        if (ahactiv(ts) > threshold(ts)) % spike!
            isspike = 1 ;
            % update threshold
            threshold(ts:ts + th_inc_length -1) = threshold(ts:ts + th_inc_length -1) + ...
                thresh_increment ;
        else
            isspike = 0 ;
        end
    end % runstep


end
function [isspike, neuron,   apicalcurrent, basalcurrent, apicalactivation, basalactivation, ...
    ahactiv,  apical, basal, shunts] = TPN_runstep(ts, tpnno, simulation, neuron,  apical, basal, shunts, ... % parameters
    apicalcurrent, basalcurrent, apicalactivation, basalactivation, ahactiv)
% run TPN simulation for one neuron for onetimestep, returning true if a spike is
% generated
% new version: standalone function, instead of being inside an
% enclosing function: started 25 Sept 2024, appears working 1 Oct 2024.
% new version (3) for use with a network. Note that apicalpostsynapse now
% part of apical, and basalpostsynapse part of basal structures. Will also
% need to be integrated into the interconnection of the neurons (network).
% started 18 Oct 2024.
% modded 29 Oct:
% currently sending wholeparameter set to and from this function, then choosing the
% specific neuron one using tpnno
% similarly, the whole of the arrays are sent and returned, but the
% specific slice for a single neuron  is selectred using the tpnno
% parameter
% Modified Nov 19 so thatsynaptic spines use the alpha function  * weight
% to provide the charge coming in which charges the capacitor
% C_apical(basal)_spine, with leakage R_synapse_spine, calculated as
% apical(basal)fracleak. Started. 
% Modified 27 Nov to consider synaptic input as charge, charging spine end
% up, and transferring a current to the dendrite, causing activity to alter
% as capacitor charged up.
%
% what to do after firing is next to be updated
%
%% parameters
% ts: timestep,
% tpnno: TPneuron number
% neuron, apical, basal, shunts parameters
% apicalpostsynapse, basalpostsynapse, apicalcurrent, basalcurrent, apicalactivation, basalactivation,
% threshold, ahactiv: supplied (and return) values for each timestep describing values
% inside neuron.
%% returns
% isspike: true if there is a spike at this timestep
% neuron:some parameters that may be altered are returned
% apicalpostsynapse, basalpostsynapse, apicalcurrent, basalcurrent,
%   apicalactivation, basalactivation, threshold,...
%   ahactiv:  returned values for variables inside neuron
% apical, basal, shunts: return parameters that may be changed

while ((apical(tpnno).apicalspikeno <= size(apical(tpnno).apicalinputs,1)) && ...
        (apical(tpnno).apicalinputs(apical(tpnno).apicalspikeno,1) == ts)) % calculate for each synapse: is there a new apical spike input?    
    % a presynaptic spike occurred at time ts on neuron tpnno at apical
    % synapse apical(tpnno).apicalinputs(apical(tpnno).apicalspikeno,2)
    %
    % v4: apical spike detected: integrate alpha function effect * weight on capacitor
    % apical(tpnno).C_apical_spine, and use
    % apical(tpnno).apicalspinefracleak to gradually discharge
    % calculate for each synapse used in this version

    % use alphafunctioneffect to calculate the spine voltage, then feed this through R_synap_dendrite to get a current into the apical dendrite
    [alphafn, alphalen] = alphafunctioneffect(apical(tpnno).alpha_apical, apical(tpnno).C_apical_spine,apical(tpnno).apicalsynapseweights(apical(tpnno).apicalinputs(apical(tpnno).apicalspikeno,2)), ...
        apical(tpnno).apicalspinefracleak, simulation.timestep, 3 * length(apical(tpnno).alpha_apical)) ;
    % add the voltage change
    apical(tpnno).apicalpostsynapse(apical(tpnno).apicalinputs(apical(tpnno).apicalspikeno,2), ts:(ts + alphalen) -1) = ...
        apical(tpnno).apicalpostsynapse(apical(tpnno).apicalinputs(apical(tpnno).apicalspikeno,2), ts:(ts + alphalen) -1) + (alphafn * ...
        apical(tpnno).synapsemultiplier) ;



    % add this change to the total depolarisation (currently no
    % geometry on apical dendrite) by feeding voltage through
    % apical(tpnno).R_synap_dendrite to get a curent
    apicalcurrent(tpnno, ts:(ts + alphalen) -1) = ...
        apicalcurrent(tpnno,ts:(ts + alphalen) -1) + ...
        apical(tpnno).apicalpostsynapse(apical(tpnno).apicalinputs(apical(tpnno).apicalspikeno,2), ts:(ts + alphalen) -1) / apical(tpnno).R_synap_dendrite ;
    %
    apical(tpnno).apicalspikeno = apical(tpnno).apicalspikeno + 1 ;
end

% calculate apical activation:
% total charge is weight coulombs!!! That's why the voltages are
% so high: fixed.
% apical current is a current (I)
% apicalactivation is a voltage: voltage increment is (i * delta t)/C
if (ts > 1)
    apicalactivation(tpnno,ts) = apicalactivation(tpnno,ts-1) * (1 - apical(tpnno).apicalfracleak) + ... % decrement due to leakage
        ((apicalcurrent(tpnno,ts) * simulation.timestep)/apical(tpnno).C_apical) ; % increment due to incoming current
end
% apply apical shunt, if any.
if ((shunts(tpnno).noapicalshunts > 0) && (isfield(shunts, 'apicalshuntinputs')))
    %   apply shunting synaptic inhibition here.
    %   for each apical shunting synapse, add up the inhibition and then apply.
    % is there a shunting synapse at this time, or has there been one
    % within the last asduration timesteps?
    while (shunts(tpnno).apicalshuntspikeno <= size(shunts(tpnno).apicalshuntinputs, 1) && ...
            (shunts(tpnno).apicalshuntinputs(shunts(tpnno).apicalshuntspikeno,1) == ts))
        if (shunts(tpnno).apicalshuntinputs(shunts(tpnno).apicalshuntspikeno,1) == ts) % if there's a new shunting input, initialise countdown
            shunts(tpnno).inapicaltimeinterval(shunts(tpnno).apicalshuntspikeno)  = true;
            % keep the id of the spike synapse number, apicalshuntinputs(apicalshuntspikeno,2)
            % and use it when indexing ap_sh_wt_pt
            shunts(tpnno).countdown_ap(shunts(tpnno).apicalshuntspikeno) = shunts(tpnno).apicalshuntduration ;
        end
        shunts(tpnno).apicalshuntspikeno = shunts(tpnno).apicalshuntspikeno + 1;
    end % while
    % calculate total shunting effect
    shunteffect = 1 ;
    for apsspikeindex = 1:length(shunts(tpnno).inapicaltimeinterval)
        if (shunts(tpnno).inapicaltimeinterval(apsspikeindex)  == true)
            shunteffect = shunteffect * shunts(tpnno).ap_sh_wt_pt(shunts(tpnno).apicalshuntinputs(shunts(tpnno).apicalshuntspikeno - 1,2)) ;
            % shunteffect = shunteffect * ap_sh_wt_pt(apsspikeindex) ; % issue: ap_sh_wt_pt see above.
            shunts(tpnno).countdown_ap(apsspikeindex) = shunts(tpnno).countdown_ap(apsspikeindex) - 1;
            if (shunts(tpnno).countdown_ap(apsspikeindex) == 0)
                shunts(tpnno).inapicaltimeinterval(apsspikeindex) = false ;
            end % if
        end% if
    end  % for
    % apply shunteffect
    apicalactivation(tpnno,ts) = apicalactivation(tpnno,ts) * shunteffect ;
end % if

% calculate basal synaptic depolarisation
while ((basal(tpnno).basalspikeno <= size(basal(tpnno).basalinputs,1)) ...
        && (basal(tpnno).basalinputs(basal(tpnno).basalspikeno,1) == ts))
    % a presynaptic spike occurred at time ts on neuron tpnno at basal
    % synapse basal(tpnno).basalinputs(basal(tpnno).basalspikeno,2)
    %
    % v4: basal spike detected: integrate alpha function effect * weight on capacitor
    % basal(tpnno).C_basal_spine, and use
    % basal(tpnno).basalspinefracleak to gradually discharge
    % calculate for each synapse used in this version
    %
    % use alphafunctioneffect to calculate the spine voltage, then feed this through R_synba_dendrite to get a current into the basal dendrite
    [alphafn, alphalen] = alphafunctioneffect(basal(tpnno).alpha_basal, basal(tpnno).C_basal_spine,basal(tpnno).basalsynapseweights(basal(tpnno).basalinputs(basal(tpnno).basalspikeno,2)), ...
        basal(tpnno).basalspinefracleak, simulation.timestep, 3 * length(basal(tpnno).alpha_basal)) ;

    % add the voltage change
    basal(tpnno).basalpostsynapse(basal(tpnno).basalinputs(basal(tpnno).basalspikeno,2), ts:(ts + alphalen) -1) = ...
        basal(tpnno).basalpostsynapse(basal(tpnno).basalinputs(basal(tpnno).basalspikeno,2), ts:(ts + alphalen) -1) + (alphafn * ...
        basal(tpnno).synapsemultiplier) ;

  
  
    % add this change to the total depolarisation (currently no
    % geometry on basal dendrite) by feeding voltage through
    % basal(tpnno).R_synba_dendrite to get a curent
    basalcurrent(tpnno, ts:(ts + alphalen) -1) = ...
        basalcurrent(tpnno,ts:(ts + alphalen) -1) + ...
        basal(tpnno).basalpostsynapse(basal(tpnno).basalinputs(basal(tpnno).basalspikeno,2), ts:(ts + alphalen) -1) / basal(tpnno).R_synba_dendrite ;

    basal(tpnno).basalspikeno = basal(tpnno).basalspikeno + 1 ;
end

% calculate basal activation:
% the basal current a current (I)
if (ts > 1)
    basalactivation(tpnno,ts) = basalactivation(tpnno,ts-1) * (1 - basal(tpnno).basalfracleak) + ...
        ((basalcurrent(tpnno,ts)*simulation.timestep)/basal(tpnno).C_basal) ;
end
if ((shunts(tpnno).nobasalshunts  > 0)  && (~isempty(shunts(tpnno).basalshuntinputs)))
    %   apply shunting synaptic inhibition here.
    %   for each basal shunting synapse, add up the inhibition and then apply
    %   apply shunting synaptic inhibition here.
    % is there a shunting synapse at this time, or has there been one
    % within the last asduration timesteps?
    while (shunts(tpnno).basalshuntspikeno <= size(shunts(tpnno).basalshuntinputs, 1) && ...
            (shunts(tpnno).basalshuntinputs(shunts(tpnno).basalshuntspikeno,1) == ts))
        if (shunts(tpnno).basalshuntinputs(shunts(tpnno).basalshuntspikeno,1) == ts) % if there's a new shunting input, initialise countdown
            shunts(tpnno).inbasaltimeinterval(shunts(tpnno).basalshuntspikeno)  = true;
            % keep the id of the spike synapse number, basalshuntinputs(apicalshuntspikeno,2)
            % and use it when indexing ap_sh_wt_pt
            shunts(tpnno).countdown_bs(shunts(tpnno).basalshuntspikeno) = shunts(tpnno).bsduration ;
        end
        shunts(tpnno).basalshuntspikeno = shunts(tpnno).basalshuntspikeno + 1;
    end % while
    % calculate total shunting effect
    shunteffect = 1 ;
    for bssspikeindex = 1:length(shunts(tpnno).inbasaltimeinterval)
        if (shunts(tpnno).inbasaltimeinterval(bssspikeindex)  == true)
            shunteffect = shunteffect * shunts(tpnno).ba_sh_wt_pt(shunts(tpnno).basalshuntinputs(shunts(tpnno).basalshuntspikeno - 1,2)) ;
            shunts(tpnno).countdown_bs(bssspikeindex) = shunts(tpnno).countdown_bs(bssspikeindex) - 1;
            if (shunts(tpnno).countdown_bs(bssspikeindex) == 0)
                shunts(tpnno).inbasaltimeinterval(bssspikeindex) = false ;
            end % if
        end% if
    end  % for
    % apply shunteffect
    basalactivation(tpnno,ts) = basalactivation(tpnno,ts) * shunteffect ;


end

% Try equation from Reza & Ahsan (R^2 + 2*R +2*C *
% (1+mod(R)))for axon hillock activation
ahactiv(tpnno,ts) = basalactivation(tpnno,ts)^2  +  2 *  basalactivation(tpnno,ts) +...
    2 * apicalactivation(tpnno,ts) * (1 + abs(basalactivation(tpnno, ts))) ;

% decide whether to spike at this timestep, and store spike if so.
if (ahactiv(tpnno,ts) > neuron(tpnno).TP_threshold(ts)) % spike!
    isspike = 1 ;
    neuron(tpnno).spikecount = neuron(tpnno).spikecount + 1 ;
    if (neuron(tpnno).spikecount) > neuron(tpnno).maxnospikes
        error("TPN_runstep_v3: maximum spike count in TPN %d exceeded at timestep %d", tpnno, ts) ;
    end
    neuron(tpnno).spikes(neuron(tpnno).spikecount) = ts ;
    % update threshold
    neuron(tpnno).TP_threshold(ts:ts + neuron(tpnno).th_inc_length -1) = neuron(tpnno).TP_threshold(ts:ts + neuron(tpnno).th_inc_length -1) + ...
        neuron(tpnno).thresh_increment ;

    % what to do after firing gets inserted here
else
    isspike = 0 ;
end
end % runstep

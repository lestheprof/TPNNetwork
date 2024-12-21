function [Cvoltage, endsample] = alphafunctioneffect(alpha,C,weight,fracleak, timestep, endsample)
%alphafunctioneffect: consider alpha function * weight as charge input, the
%integrate this on capacitor C, while leaking using fracleak.
%   The charge incoming is the weight times the alpha function: calculated timestep
%   by timestep. This is integrated on C, whilst leaking away through R
%   (calculated as fracleak earlier)
% alpha is alphafunction (array)
% weight is the weight at this spine synapse
% endsample is the length of this alphaeffect (which is really a voltage
% trace), necessarily truncated at endsample
% fracleak precomputed from R and C at the synapse in setupnetwork.m
% started LSS 19 11 2024.
Cvoltage = zeros([1 endsample]) ; % preallocate
alphalength = length(alpha) ;
Cvoltage(1) = (alpha(1)*weight*timestep)/C ;
for t = 2:endsample
    if (t <= alphalength) % add charge from alphafunction
        Cvoltage(t) = Cvoltage(t-1) + (alpha(t)*weight*timestep)/C ; % V=(Q(delta T))/C
    else
        Cvoltage(t) = Cvoltage(t-1) ;
    end
    % and let some leak away
    Cvoltage(t) = Cvoltage(t)* (1-fracleak);
end
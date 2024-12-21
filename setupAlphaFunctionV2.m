function alphafunction = setupAlphaFunctionV2(timestep,tau)
            % how long does the alpha function need to be stored for?
            % according to Koch (p100), alpha function reaches 1 percent of
            % the peak at 7.64 * tpeak, and tpeak = tau. So we'll use 10 *
            % tau for a little safety.
            % SetUpAlphaFunctionV2 is used in v4, where the alpha function
            % describes the transfer of charge through  the synapse, so
            % instead of summing to 1, the integral (sum) of the function
            % over time is 1
            %
            % LSS last updated 16 July 2024
            %
            alphalength = round(10 * tau / timestep) ;
            alphafunction = zeros([1  alphalength]) ;
            for i=1:alphalength
                alphafunction(i) = (i*timestep/tau) * exp(1 - (i* timestep)/tau) ;
            end
            % normalise sum of alphafunction elements to 1
            % alphafunction = alphafunction/sum(alphafunction) ; % normalise
            % but is this right? if we use the alpha function to charge a
            % capacitor, we would want the sum \sum_{lenght of alpha
            % function} alpha(i) * timestep to be constant, independent of
            % timestep. Normalising as above gives this sum to be equal to
            % timestep.
            % This can be achieved by multiplying the denominator by the
            % timestep
            alphafunction = alphafunction/(sum(alphafunction) * timestep) ;
        end
function cohcs = calc_guardrail_cohc(cfs, species, guardrail_mode)
% cohcs = CALC_GUARDRAIL_COHC(cfs, species, guardrail_mode)
%
% Calculates a "guardrail" COHC value at specified CF values.
% 
% Based on a series of simulations reported in [[cite]], guardrail COHC
% values are values of COHC for a specific CF and species that prevent ΔL
% and ΔTH values due to efferent activity from exceeding what is plausible
% given MOC electrical stimulation experiments. These simualations were
% parametrized in terms of normalized cochlear distance according to:
%     Greenwood, D. D. (1990). A cochlear frequency-position function for
%     several species—29 years later. The Journal of the Acoustical Society
%     of America, 87(6), 2592–2605. https://doi.org/10.1121/1.399052
%
% Hence, below we convert CF to a cochlear distance based on species and
% then determine an appropriate COHC value based on a pre-determined linear
% mapping between cochlear distance and log(COHC).

% Convert from CF to cochlear distance based on species
if species == 1  % cat
	cd = log10(cfs ./ 456.0 + 0.8) ./ 2.1;
elseif species == 2  % human
	cd = log10(cfs ./ 165.4 + 0.88) ./ 2.1;
else
	error("Species settings other than 1 (cat) or 2 (Human Shera) are not supported.");
end

% Return COHC based on mode and cochlear distance
% If "none", we have no guardrails, so lowest allowed COHC is 0.0
if guardrail_mode == "none"
	cohcs = zeros(size(cfs));
% if "standard", use lm between cochlear distance and log(cohc)
elseif guardrail_mode == "standard"
	cohcs = exp(-2.003942 + 1.753976 .* cd);
elseif guardrail_mode == "standard+5"
	cohcs = exp(-2.452641 + 1.942155 .* cd);
elseif guardrail_mode == "standard+10"
	cohcs = exp(-3.010872 + 2.252671 .* cd);
elseif guardrail_mode == "standard+15"
	cohcs = exp(-3.665091 + 2.659225 .* cd);
elseif guardrail_mode == "standard+20"
	cohcs = exp(-4.587940 + 3.358810 .* cd);
else
	error("guradrail_mode not recognized!")
end



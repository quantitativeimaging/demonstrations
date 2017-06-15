% % EJR 2017
% 
% % 1. Generate one (or more) decay times
% %    Assume 1 excited fluorophore
% % Problems: uniform prior is weighted to too-high numbers if few data
% %           Likelihood hits small-number range over about 200 fluorophores - use log
% % likelihood.
% 
% NOTES - CHECK ADDITION OF LIKELIHOODS! (15-6-2017)
%       - I think the addition is sometimes log when it might be linear
%       - Probably under-weights terms with missing observations. 
%       - 
%       - Probably fixed (15/6)
%       - But am I now over-counting miss chance
%       - second one should not be (m-1)/N_D
%       - But a correction to (m-1)/N_D from (L/N_D). Poss (m-1)/L


number_of_simulations = 20;
EST_TAUS = zeros(number_of_simulations, 1);
for lpGrand = 1:number_of_simulations
	
% To run this: 
%  1. Generate some tObs by running section 3.

% Detection times - given tau ~ 4.1, Q = 1, n_D = 32.
% The list below is one such series of 9 detection times (10 arrivals)

% tObs  = [0.3820; 0.4403; 1.7521; 2.1195; 2.2572; 2.4480; 3.1544; 3.7733; 10.3995]; % From tau=4.1, 10 arrivals
% tObs  = [3.6689; 5.4470; % From tau = 4.1, 5 arrivals
% tObs  = [0.8508; 1.0155; 1.1545; 1.5293; 3.3626]; % 5 arrivals, tau = 4.1
% tObs  = [1.5441; 1.8589; 2.8371; 4.1173] % from 5 arrivals, tau=4.1. Gives nice result...
%   tObs    = [1.1568; 1.6567; 2.0638; 3.3979; 4.1945]; % from 5 arrivals, tau=4.1. Gives nice result...

% 1. PHOTON DYNAMICS SIMULATION:
% OR SIMULATE SOME DETECTION TIMES:
tObsList  = [];   % Detected photon times
n_D       = 32;   % Number of detector elements
T         = 25;   % Upper time limit of detection
Q         = 1.00; % Constant quantum efficiency term
tau_known = 4.1;  % Defined Lifetime, ns
lambda_known = 1/tau_known; % in case 1/tau is useful
Nf        = 10; % number of excited fluorophores
tEms = exprnd(tau_known, Nf,1); 
tEmsSorted = sort(tEms);
listDetRands = rand(size(tEmsSorted));
m = 0; % Considering what happens when the m'th photon has been detected
listDetected = zeros(size(tEmsSorted)); % Is this photon detected?
for lp = 1:length(tEmsSorted)
	detector_efficiciency = Q * (1 - m/n_D);
	if(listDetRands(lp) <= detector_efficiciency)
		listDetected(lp) = 1;
		m = m+1;
	end
end
tObs = tEmsSorted(listDetected==1);
tObs(tObs > T) = [];


% 2. INFERENCE METHOD
L   = length(tObs); % Number of detected photons
t_L = tObs(end);    % Final detected time (followed by blank up to T)
			 
Q = 1.00; % Constant quantum efficiency term (should be a known value)
T = 25;   % nanoseconds. Upper time-limit of detector
N_D = 32; % Number of detector elements. Known.

poss_ns = [L:(ceil(2*L))];
poss_taus = 1:0.2:20;
[array_ns, array_taus] = meshgrid(poss_ns, poss_taus); 

list_ns = array_ns(:);
list_taus = array_taus(:);
list_log_likelihood_total = zeros(size(list_ns));
% list_log_fmds = zeros([ length(tObs), length(list_ns)] );

for lpPrior = 1:length(list_ns)

% Get value of priors for this likelihood calculation
n   = list_ns(lpPrior);
tau = list_taus(lpPrior);
lam = 1/tau;

% Initialise log-likelihood value 
% ... Start with likelihood = 1 for this value of (n, lam)
log_likelihood_total = 0; % 
likelihood_jLs = 0;

% WORK OUT PHOTON DETECTION 'SLICING' PROBABILITY
% Arrivals: (L+jL), jT, (n-L-jL-jT) in times (0,t_L), (t_L, T), (T,inf) 
% Not-detected: (jL, jT, n-L-jL-jT) in times   "         "         "
% For long lifetime, probability that <<L arrive in mid-range is low
% ?Issue - Since we now have jL, we may need to redo j math in t_m_d's
%  Work with jL misses during detection phase, within each slice...
for jL = 0:(n-L)
	log_likelihood_total = 0; % For this value of jL...
	f_slicing = 0 ; % Probability for a particular value of jL. Initialise
	for jT = 0:(n-L-jL)
	% Combinations / arrival probs / necessary miss chances
	% Combinations: n! / n1! n2! n3!
	% Before t_L; between tL and T, after T
	f_gap_inc = ( factorial(n) / ( factorial(L+jL)*factorial(jT)*factorial(n-L-jL-jT) ) ) * ... 
		(1-exp(-lam*t_L))^(L+jL) * (exp(-lam*t_L) - exp(-lam*T))^(jT) * ...
		(exp(-lam*T))^(n-L-jL-jT)* ...
		(1-(N_D-L)/N_D)^(jL+jT);
	f_slicing = f_slicing + f_gap_inc;
	end
	% Have added up all p( jT and miss) values that produce this jL
	% NOW THERE ARE L + jL photons arriving in the period (0,t_L)
	% ? We have already accounted for the probability that the jL are missed
	log_likelihood_total = log_likelihood_total + log( f_slicing );

	% WORK OUT THE Probability of the m'th photon detection time, 
	% Given that there are n arrivals in t=(0,inf) and 
	% ... a maximum of jL photons can be missed during the (0,tL) range
	for m = 1:(length(tObs) - 0);
		t_m_d = tObs(m); % The m'th (present) detection time
		
		if( m==1 )   % FOR THE FIRST DETECTION. If Q=1, no arrivals can be missed before t = t_m_d
			f_m_d = 0; % Probability density for m'th detection. Initialise
			for k = 1:(1+jL) % skipping k-1 arrivals, detecting k'th
				% Probability of k'th arrival at time t_m_d, given n arrivals
				% f_n_a = lam * (L+jL) * nchoosek(L+jL-1,k-1) * (1-exp(-lam*DT))^(k-1) * ...
				% 	exp( -(L+jL-k+1)*lam*DT );
				f_n_a = lam * (n) * nchoosek(n-1,k-1) * (1-exp(-lam*t_m_d))^(k-1) * ...
					exp( -(n-k+1)*lam*t_m_d );

				% Increment for f_m_d. The (1) is for (perfect) spatial_efficiency
				f_inc = Q*((1-Q)^(k-1)) * (1) * f_n_a;
				f_m_d = f_m_d + f_inc;
			end
		log_fmd = log(f_m_d);
		log_likelihood_total = log_likelihood_total + log_fmd;
		% list_log_fmds(1, lpPrior) = log_fmd;
		end

		if(m>1) % FOR THE SECOND DETECTION AND ONWARDS
		f_m_d = 0;
		for j = 0:jL % j = # misses up to the m'th detection. Photons arriving in (0,t_m) inclusive is m + up to jL
		  % Sum up: Prob miss j photons * Prob detect k'th photon, k=m+j 
		  f_inc = nchoosek((n-1),(j+m-1)) * ... 
				      lam*(n)*(1-exp(-lam*t_m_d))^(j+m-1)*exp(-(n-j-m+1)*lam*t_m_d) * ...  
						  Q* ( 1 - Q*(N_D+1-m)/N_D )^j;
			f_m_d = f_m_d + f_inc;
		end
		log_fmd = log(f_m_d);
		log_likelihood_total = log_likelihood_total + log_fmd;
		end
	end % Loop contains calculation of m'th photon detection time probability
  likelihood_jLs = likelihood_jLs + exp(log_likelihood_total);
	
end % End of loop containing this jL value (jL photons missed before t_L)
list_log_likelihood_total(lpPrior) = log(likelihood_jLs);
% list_log_likelihood_total(lpPrior) = log_likelihood_total;

end % End of loop through this value of the prior (n, lam)

array_log_likelihood = reshape(list_log_likelihood_total, size(array_ns));

figure(5)
mesh(array_ns, array_taus, exp(array_log_likelihood));
xlabel('n');
ylabel('tau')
zlabel('likelihood')
set(gca, 'fontSize', 14)
title('Imperfect inference of (n,\Tau) given detection times')
% % PUT PRIOR END STATEMENTS HERE

p_marginal = sum(sum(exp(array_log_likelihood)));
est_tau = sum(sum(exp(array_log_likelihood).*array_taus)) / p_marginal

figure(6)
plot(poss_taus, (p_post_taus_log) );
xlabel('lifetime / ns', 'fontSize', 14)
ylabel('likelihood / arb', 'fontSize', 14)


p_post_taus = sum(exp(array_log_likelihood), 2) / p_marginal;
% p_post_taus_total = p_post_taus_total + p_post_taus;
p_post_taus_log   = p_post_taus_log + log(p_post_taus);
% p_post_taus_log = log(p_post_taus);

figure(7)
plot(poss_taus, p_post_taus / (sum(p_post_taus) ));
xlabel('lifetime (discretised) / ns', 'fontSize', 14)
ylabel('likelihood / normalised', 'fontSize', 14)
title('Likelihood - if photons may miss', 'fontSize', 14)
set(gca, 'fontSize', 14)
% ylim([0 0.12])
drawnow

figure(8)
hist(tObs, [0.5:1:24.5])
xlabel('time / ns', 'fontSize', 14)
ylabel('count', 'fontSize', 14)
title('Detected times', 'fontSize', 14)
set(gca, 'fontSize', 14)

% figure(9)
% mesh(list_log_fmds);

% % % %  REPEAT FOR MANY PULSES:
% % Repeat simulation and inference:
% % p_post_taus = sum(exp(array_log_likelihood), 2) / p_marginal;
% p_post_grand = p_post_taus;
% p_post_grand = p_post_grand + p_post_taus;
% 


EST_TAUS(lpGrand) = est_tau;
end
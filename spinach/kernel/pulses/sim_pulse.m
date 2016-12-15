% Simultaneous multi-channel pulse function. Resamples the grids of the input
% waveforms onto a common time grid, generates a unified Liouvillian stack on
% the resulting common grid and applies the resulting propagators to the state
% vector. The pulse centers are aligned in the resulting combined waveform.
% 
% Each nucleus can only appear once on the nucleus list: using e.g. two diffe-
% rent carbon channels is not allowed.
%
% Input variables:
%
%    L:             the drift Liouvillian
%
%    offsets:       [off1 off2 ...]           (Hz, row vector)
%
%    phases:       {[phi1 phi2 ... phiN],...
%                   [phi1 phi2 ... phiK]}     (degrees, cell array of row vectors)
%
%    amplitudes:   {[a1 a2 ... aN],...
%                   [b1 b2 ... bK]}           (rad/s, cell array of row vectors)
%
%    durations:     [dur1 dur2 ...]           (seconds, row vector)
%
%    spins:        {'1H', '13C', ...}         (cell array of strings)
%
% Note: the function is only applicable in situations where the drift
%       Liouvillian commutes with offset operators. 
%
% i.kuprov@soton.ac.uk

function rho=sim_pulse(spin_system,L,rho,spins,offsets,phases,amplitudes,durations)

% Check consistency
grumble(spin_system,spins,offsets,phases,amplitudes,durations)

% Find out which pulse is the longest in the stack
longest_pulse_duration=max(durations);

% Find out which pulse has got the largest frame density
nframes=zeros(size(durations));
for n=1:length(amplitudes)
    nframes(n)=numel(amplitudes{n});
end
densest_pulse_fps=max(nframes./durations);

% Generate a time grid that is at least as long as the longest pulse and
% at least as dense as the pulse with the largest frame density.
new_grid_duration=longest_pulse_duration;
new_grid_fps=densest_pulse_fps;
new_grid_nframes=ceil(new_grid_duration*new_grid_fps);
new_grid_step=new_grid_duration/new_grid_nframes;
new_grid=new_grid_duration*((0:(new_grid_nframes-1))-(new_grid_nframes-1)/2)/(new_grid_nframes-1);

% Preallocate the new amplitude and phase arrays
new_amplitudes=cell(size(amplitudes));
new_phases=cell(size(phases));

% Resample all pulses onto the new time grid.
for n=1:length(offsets)
    old_grid=durations(n)*((0:(nframes(n)-1))-(nframes(n)-1)/2)/(nframes(n)-1);
    new_amplitudes{n}=interp1(old_grid,amplitudes{n},new_grid,'pchip',0);
    new_phases{n}=interp1(old_grid,phases{n},new_grid,'pchip',0);
end

% Pre-compute primitive operators for all channels
Lx=cell(size(spins)); Ly=cell(size(spins)); Lz=cell(size(spins));
for n=1:length(spins)
    Lm=operator(spin_system,'L-',spins{n});
    Lx{n}=(Lm'+Lm)/2; Ly{n}=(Lm'-Lm)/2i; Lz{n}=(Lm'*Lm-Lm*Lm')/2;
end

% Apply the pulse stack
for n=1:new_grid_nframes
    
    % Concatenate operators across RF channels at the current grid step
    total_operator=L;
    for k=1:length(spins)
        pulse_operator=new_amplitudes{k}(n)*(Lx{k}*cosd(new_phases{k}(n))+Ly{k}*sind(new_phases{k}(n)));
        total_operator=total_operator+pulse_operator-2*pi*offsets(k)*Lz{k};
    end
   
    % Apply the pulse chunk
    rho=step(spin_system,total_operator,rho,new_grid_step);
    
end

end

% Consistency enforcement
function grumble(spin_system,spins,offsets,phases,amplitudes,durations)
if ~ischar(spins)
    error('spins parameter must be a character string.');
end
if ~ismember(spins,spin_system.comp.isotopes)
    error('the requested spins are not present in the spin system.');
end
if (~isnumeric(offsets))||(~isreal(offsets))
    error('offsets parameter must be a vector of real numbers.');
end
if (~isnumeric(durations))||(~isreal(durations))||any(durations<=0)
    error('durations parameter must be a vector of positive real numbers.');
end
if (~isnumeric(phases))||(~isreal(phases))
    error('phases parameter must be an array of real numbers.');
end
if (~isnumeric(amplitudes))||(~isreal(amplitudes))
    error('amplitudes parameter must be an array of real numbers.');
end
if all(size(phases)~=size(amplitudes))
    error('phases and amplitudes arrays must have the same number of elements.');
end
end

% Had I the heavens' embroidered cloths,
% Enwrought with golden and silver light,
% The blue and the dim and the dark cloths
% Of night and light and the half light,
% I would spread the cloths under your feet:
% But I, being poor, have only my dreams;
% I have spread my dreams under your feet;
% Tread softly because you tread on my dreams.
%
% W.B. Yeats (1865–1939)


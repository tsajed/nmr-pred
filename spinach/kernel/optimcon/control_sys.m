% Options and structure creation for a control system. Sets the various
% options used throughout optimal control algorithms.
%
% Modifications to this function are discouraged -- the accuracy settings
% should be modified by setting the sub-fields of the ctrl_param structure.
%
% <http://spindynamics.org/wiki/index.php?title=Control_sys.m>

function ctrl_param=control_sys(spin_system,control_system)

% list of CONTROL algorithms
% >>>>-add to lists as appropriate-<<<<
optimcon_algos={'grape','tannor','zhu-rabitz'};

% ========================== STRUCT = OPTIM.SYS ==========================
% store the cost function handle; throw error is not supplied - should
% always be supplied if used with fminnewton: inherited from fminnewton
% inputs

% control banner
if isfield(control_system,'method') && ismember(control_system.method,optimcon_algos)
    ctrl_param.method=control_system.method;
    control_system=rmfield(control_system,'method'); % parsed -> remove
else
    error('must provide a control method e.g. grape, tannor, or zhu-rabitz')
end

optim_report(spin_system.sys.output,' ');
optim_report(spin_system.sys.output,'===========================================');
optim_report(spin_system.sys.output,'=                                         =');
optim_report(spin_system.sys.output,'=                 CONTROLS                =');
optim_report(spin_system.sys.output,'=                                         =');
optim_report(spin_system.sys.output,'===========================================');
optim_report(spin_system.sys.output,' ');
optim_report(spin_system.sys.output,[pad('Control algorithm',80) pad(ctrl_param.method,20) ' ']);

controls=control_system.control_ops; num_ctrls=numel(controls);
ctrl_param.control_ops=controls;
ctrl_param.comm_rel=zeros(num_ctrls);
for n=1:num_ctrls
    for m=1:num_ctrls
        ctrl_param.comm_rel(n,m)=1-any(any(controls{n}*controls{m}-controls{m}*controls{n}));
    end
end
optim_report(spin_system.sys.output,[pad('Number of control operators',80) pad(int2str(num_ctrls),20) ' (calculated)']);
control_system=rmfield(control_system,'control_ops'); % parsed -> remove
optim_report(spin_system.sys.output,[pad('Pairwise non-commuting control operators',80) pad(int2str((sum(ctrl_param.comm_rel(:)==0))./2),20) ' (calculated)']);
% test for hermitian control_ops
if ~all(cellfun(@ishermitian,ctrl_param.control_ops))
    optim_report(spin_system.sys.output,pad('WARNING: not all control operators are Hermitian',80));
end
if isfield(control_system,'target')
    ctrl_param.target=control_system.target;
    control_system=rmfield(control_system,'target'); % parsed -> remove
end
if isfield(control_system,'rho')
    ctrl_param.rho=control_system.rho;
    control_system=rmfield(control_system,'rho'); % parsed -> remove
end
if isfield(control_system,'ref_waveform') && ~ismember(ctrl_param.method,{'grape'})
    ctrl_param.ref_waveform=control_system.ref_waveform;
    control_system=rmfield(control_system,'ref_waveform'); % parsed -> remove
    optim_report(spin_system.sys.output,[pad('Reference waveform provided',80) pad('true',20) ' (user-specified)']);
elseif ~ismember(ctrl_param.method,{'grape'})
    optim_report(spin_system.sys.output,[pad('Reference waveform provided',80) pad('false',20) ' (default)']);
end
if isfield(control_system,'power_level') && control_system.power_level>0
    ctrl_param.power_level=control_system.power_level;
    control_system=rmfield(control_system,'power_level'); % parsed -> remove
    optim_report(spin_system.sys.output,[pad('Control pulse control amplitude (rad/s)',80) pad(num2str(ctrl_param.power_level,'%0.6g'),20) ' (user-specified)']);
end
if ismember(ctrl_param.method,{'tannor','zhu-rabitz'})
    if isfield(control_system,'krotov_microiteration')
        ctrl_param.krotov_microiteration=control_system.krotov_microiteration;
        control_system=rmfield(control_system,'krotov_microiteration'); % parsed -> remove
        optim_report(spin_system.sys.output,[pad('Krotov microiteration tolerance',80) pad(num2str(ctrl_param.krotov_microiteration,'%0.8e'),20) ' (user-specified)']);
    else
        ctrl_param.krotov_microiteration=1e-6;
        optim_report(spin_system.sys.output,[pad('Krotov microiteration tolerance',80) pad(num2str(ctrl_param.krotov_microiteration,'%0.8e'),20) ' (safe default)']);
    end
    
    if isfield(control_system,'lambda')
        ctrl_param.lambda=control_system.lambda;
        control_system=rmfield(control_system,'lambda'); % parsed -> remove
        optim_report(spin_system.sys.output,[pad('krotov lagrange multiplier, lambda',80) pad(num2str(ctrl_param.lambda,'%0.8e'),20) ' (user-specified)']);
    else
        ctrl_param.lambda=1e-3;
        optim_report(spin_system.sys.output,[pad('krotov lagrange multiplier, lambda',80) pad(num2str(ctrl_param.lambda,'%0.8e'),20) ' (safe default)']);
    end
end
if isfield(control_system,'pulse_duration') && control_system.pulse_duration>0
    ctrl_param.pulse_duration=control_system.pulse_duration;
    control_system=rmfield(control_system,'pulse_duration'); % parsed -> remove
    optim_report(spin_system.sys.output,[pad('Total duration of control pulses (s)',80) pad(num2str(ctrl_param.pulse_duration,'%0.6g'),20) ' (user-specified)']);
    if isfield(control_system,'nsteps') && mod(1,control_system.nsteps+1) && control_system.nsteps>0
        ctrl_param.nsteps=control_system.nsteps;
        control_system=rmfield(control_system,'nsteps'); % parsed -> remove
        optim_report(spin_system.sys.output,[pad('Number of control pulse',80) pad(int2str(ctrl_param.nsteps),20) ' (user-specified)']);
        if isfield(control_system,'time_step') && control_system.time_step==ctrl_param.pulse_duration/ctrl_param.nsteps
            ctrl_param.time_step=control_system.time_step;
            control_system=rmfield(control_system,'time_step'); % parsed -> remove
            optim_report(spin_system.sys.output,[pad('Length of single control pulse (s)',80) pad(num2str(ctrl_param.time_step,'%0.6g'),20) ' (user-specified)']);
        end
    end
end

%optim_report(optim_param.sys.output,[pad('Penalty function type',80) pad(int2str(1),20) 'TODO']);


unprocessed=fieldnames(control_system);
if ~isempty(unprocessed)
    optim_report(spin_system.sys.output,'----------------------------------------------------------------------------------------------');
    for n=1:numel(unprocessed)
        optim_report(spin_system.sys.output,[pad('WARNING: unprocessed option',80) pad(unprocessed{n},20)]);
    end
    optim_report(spin_system.sys.output,'----------------------------------------------------------------------------------------------');
end



end


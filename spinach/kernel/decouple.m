% Obliterates all interactions and populations in the subspace of states
% that involve user-specified spins in any way. The specified spins would
% not contribute to the system dynamics until the Liouvillian is rebuilt
% from scratch. Syntax:
%
%                 [L,rho]=decouple(spin_system,L,rho,spins)
%
%     spins: spins to be wiped, specified either by name, such as '13C',
%            or by a list of numbers, such as [1 2 3].
%
%     L:     Liouvillian superoperator
%
%     rho:   state vector or a horizontal stack thereof
%
% Note: this function is an analytical equivalent of a running decoupling
%       pulse sequence on the specified spins.
%
% Note: this function requires sphten-liouv formalism.
%
% i.kuprov@soton.ac.uk

function [L,rho]=decouple(spin_system,L,rho,spins)

% Return if the spin list is empty
if isempty(spins), return; end

% Check consistency
grumble(spin_system,L,rho,spins);

% Find the nuclei to be decoupled
if isnumeric(spins)
    dec_mask=false(1,spin_system.comp.nspins);
    dec_mask(spins)=true;
else
    dec_mask=ismember(spin_system.comp.isotopes,spins);
end

% Inform the user
report(spin_system,[num2str(nnz(dec_mask)) ' spins to be frozen and depopulated.']);

% Get the list of states to be wiped
zero_mask=(sum(spin_system.bas.basis(:,dec_mask),2)~=0);

% Inform the user
report(spin_system,['zeroing ' num2str(nnz(zero_mask)) ' rows and columns in the Liouvillian.']);

% Zero the corresponding rows and columns of the Liouvillian
L(zero_mask,:)=0; L(:,zero_mask)=0;

% Zero the corresponding rows of the state vector stack
if nargout==2
    report(spin_system,['zeroing ' num2str(nnz(zero_mask)) ' rows in the state vector.']);
    rho(zero_mask,:)=0;
end

end

% Consistency enforcement
function grumble(spin_system,L,rho,spins)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('analytical decoupling is only available for sphten-liouv formalism.');
end
if (~isempty(rho))&&(size(L,2)~=size(rho,1))
    error('matrix dimensions of L and rho must agree.');
end
if size(L,1)~=size(L,2)
    error('Liouvillian must be a square matrix.');
end
if (~isnumeric(spins))&&(~iscell(spins))
    error('spins parameter must either be a list of numbers or a cell array of strings.');
end
if iscell(spins)&&any(~ismember(spins,spin_system.comp.isotopes))
    error('the system does not contain the spins specified.');
end
if isnumeric(spins)&&(size(spins,1)~=1)
    error('if spins are specified by number, a row vector of numbers must be used.');
end
if isnumeric(spins)&&(max(spins)>spin_system.comp.nspins)
    error('the spin number specified is greater than the number of spins in the system.');
end
if isnumeric(spins)&&(any(~isreal(spins))||any(spins<1))
    error('spin numbers must be real positive integers.');
end
end

% It's not worth doing something unless you were doing something that
% someone, somewere, would much rather you weren't doing.
%
% Terry Pratchett 


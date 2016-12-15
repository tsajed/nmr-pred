% Splits the relaxation superoperator into longitudinal, transverse and
% mixed components.
%
% i.kuprov@soton.ac.uk

function [R1,R2,Rm]=relax_split(spin_system,R)

% Interpret the basis
[L,M]=lin2lm(spin_system.bas.basis);

% Index single- and multi-spin orders (sso and mso)
sso_mask=(sum(logical(spin_system.bas.basis),2)==1);
mso_mask=(sum(logical(spin_system.bas.basis),2)>1 );

% Index longlitudinal and transverse states
long_sso_mask=any((L>0)&(M==0),2)&sso_mask;
tran_sso_mask=any((L>0)&(M~=0),2)&sso_mask;

% Split the relaxation superoperator
R1=R; R1(~long_sso_mask,~long_sso_mask)=0;
R2=R; R2(~tran_sso_mask,~tran_sso_mask)=0;
Rm=R; R2(~mso_mask,~mso_mask)=0;

end

% He that reproveth a scorner getteth to himself shame: and he that
% rebuketh a wicked man getteth himself a blot. Reprove not a scorn-
% er, lest he hate thee: rebuke a wise man, and he will love thee.
% Give instruction to a wise man, and he will be yet wiser: teach a
% just man, and he will increase in learning.
%
% Proverbs 9:7


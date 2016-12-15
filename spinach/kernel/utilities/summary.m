% Prints various summaries on behalf of the spin system setup modules
% of Spinach kernel. Edits to this function are discouraged.
%
% i.kuprov@soton.ac.uk

function summary(spin_system,topic,header)

switch topic
    
    case 'zeeman'
        report(spin_system,header);
        report(spin_system,'==============================================================================================');
        report(spin_system,'#    Spin  2S+1                 Matrix                  norm(rank0)  norm(rank1)  norm(rank2) ');
        report(spin_system,'----------------------------------------------------------------------------------------------');
        for n=1:spin_system.comp.nspins
            if significant(spin_system.inter.zeeman.matrix{n},0)
                [rank0,rank1,rank2]=mat2sphten(spin_system.inter.zeeman.matrix{n});
                report(spin_system,[pad(num2str(n),6) pad(spin_system.comp.isotopes{n},6) pad(num2str(spin_system.comp.mults(n)),6)...
                                    num2str(spin_system.inter.zeeman.matrix{n}(1,:),'%+10.3e  %+10.3e  %+10.3e')]);
                report(spin_system,['                  ' num2str(spin_system.inter.zeeman.matrix{n}(2,:),'%+10.3e  %+10.3e  %+10.3e')...
                                    '    '  pad(num2str(norm(rank0),'%10.4e'),13) pad(num2str(norm(rank1),'%10.4e'),13) num2str(norm(rank2),'%10.4e')]);
                report(spin_system,['                  ' num2str(spin_system.inter.zeeman.matrix{n}(3,:),'%+10.3e  %+10.3e  %+10.3e')]);
                if n<spin_system.comp.nspins, report(spin_system,''); end
            end
        end
        report(spin_system,'==============================================================================================');
    
    case 'coordinates'
        report(spin_system,header);
        report(spin_system,'======================================');
        report(spin_system,'N    Spin     X         Y         Z   ');
        report(spin_system,'--------------------------------------');
        for n=1:spin_system.comp.nspins
            report(spin_system,[strjust([num2str(n) blanks(3-length(num2str(n)))],'left') ' '...
                                strjust([spin_system.comp.isotopes{n} blanks(5-length(spin_system.comp.isotopes{n}))],'center') '  '...
                                num2str(spin_system.inter.coordinates{n},'%+5.3f    ') '  ' spin_system.comp.labels{n}]);
        end
        report(spin_system,'======================================');
        
    case 'pbc'
        report(spin_system,header);
        report(spin_system,'===============================');
        report(spin_system,'     X         Y         Z     ');
        report(spin_system,'-------------------------------');
        for n=1:numel(spin_system.inter.pbc)
            report(spin_system,['   ' pad(num2str(spin_system.inter.pbc{n}(1),'%+5.3f   '),10)...
                                      pad(num2str(spin_system.inter.pbc{n}(2),'%+5.3f   '),10)...
                                      pad(num2str(spin_system.inter.pbc{n}(3),'%+5.3f   '),10)]);
        end
        report(spin_system,'===============================');
        
    case 'couplings'
        report(spin_system,header);
        report(spin_system,'=============================================================================================');
        report(spin_system,'Spin A  Spin B                 Matrix                 norm(rank0)  norm(rank1)  norm(rank2)  ');
        report(spin_system,'---------------------------------------------------------------------------------------------');
        [rows,cols,~]=find(cellfun(@norm,spin_system.inter.coupling.matrix)>spin_system.tols.inter_cutoff);
        for n=1:numel(rows)
            [rank0,rank1,rank2]=mat2sphten(spin_system.inter.coupling.matrix{rows(n),cols(n)});
            report(spin_system,['  ' pad(num2str(rows(n)),6) '  ' pad(num2str(cols(n)),6) ...
                                num2str(spin_system.inter.coupling.matrix{rows(n),cols(n)}(1,:),'%+10.3e  %+10.3e  %+10.3e')]);
            report(spin_system,['                ' num2str(spin_system.inter.coupling.matrix{rows(n),cols(n)}(2,:),'%+10.3e  %+10.3e  %+10.3e')...
                                '    '  pad(num2str(norm(rank0),'%10.4e'),13) pad(num2str(norm(rank1),'%10.4e'),13) num2str(norm(rank2),'%10.4e')]);
            report(spin_system,['                ' num2str(spin_system.inter.coupling.matrix{rows(n),cols(n)}(3,:),'%+10.3e  %+10.3e  %+10.3e')]);
            if n<numel(rows), report(spin_system,''); end
        end
        report(spin_system,'=============================================================================================');
         
    case 'chemistry'
        if numel(spin_system.chem.parts)>1
            for n=1:numel(spin_system.chem.parts)
                report(spin_system,['chemical subsystem ' num2str(n) ' contains spins: ' num2str(spin_system.chem.parts{n})]);
            end
            report(spin_system,'inter-subsystem reaction rates:');
            report(spin_system,'===============================');
            report(spin_system,' N(from)   N(to)    Rate(Hz)   ');
            report(spin_system,'-------------------------------');
            [rows,cols,vals]=find(spin_system.chem.rates);
            for n=1:length(vals)
                report(spin_system,[' ' strjust([num2str(rows(n)) blanks(3-length(num2str(rows(n))))],'left') '       '...
                                        strjust([num2str(cols(n)) blanks(3-length(num2str(cols(n))))],'left') '      '...
                                                 num2str(vals(n),'%+0.3e')]);
            end
            report(spin_system,'===============================');
        end
        [rows,cols,vals]=find(spin_system.chem.flux_rate);
        if numel(vals)>0
            report(spin_system,'point-to-point flux rates:');
            report(spin_system,'===============================');
            report(spin_system,' N(from)   N(to)    Rate(Hz)   ');
            report(spin_system,'-------------------------------');
            for n=1:length(vals)
                report(spin_system,[' ' strjust([num2str(rows(n)) blanks(3-length(num2str(rows(n))))],'left') '       '...
                                        strjust([num2str(cols(n)) blanks(3-length(num2str(cols(n))))],'left') '      '...
                                                 num2str(vals(n),'%+0.3e')]);
            end
        end

    case 'rlx_rates_t1_t2'
        report(spin_system,header);
        report(spin_system,'========================================');
        report(spin_system,'N    Spin        R1             R2      ');
        report(spin_system,'----------------------------------------');
        for n=1:spin_system.comp.nspins
            report(spin_system,[strjust([num2str(n) blanks(3-length(num2str(n)))],'left') ' '...
                                strjust([spin_system.comp.isotopes{n} blanks(5-length(spin_system.comp.isotopes{n}))],'center') '  '...
                                num2str(spin_system.rlx.r1_rates(n),'%+0.5e   ') '  '...
                                num2str(spin_system.rlx.r2_rates(n),'%+0.5e   ') '  '...
                                spin_system.comp.labels{n}]);
        end
        report(spin_system,'========================================');
        
    case 'rlx_rates_lindblad'
        report(spin_system,header);
        report(spin_system,'========================================');
        report(spin_system,'N    Spin        R1             R2      ');
        report(spin_system,'----------------------------------------');
        for n=1:spin_system.comp.nspins
            report(spin_system,[strjust([num2str(n) blanks(3-length(num2str(n)))],'left') ' '...
                                strjust([spin_system.comp.isotopes{n} blanks(5-length(spin_system.comp.isotopes{n}))],'center') '  '...
                                num2str(spin_system.rlx.lind_r1_rates(n),'%+0.5e   ') '  '...
                                num2str(spin_system.rlx.lind_r2_rates(n),'%+0.5e   ') '  '...
                                spin_system.comp.labels{n}]);
        end
        report(spin_system,'========================================');
        
    case 'rlx_rates_nott'
        report(spin_system,' ');
        report(spin_system,'=== Nottingham DNP relaxation theory ===');
        report(spin_system,['Electron R1: ' num2str(spin_system.rlx.nott_r1e) ' Hz']);
        report(spin_system,['Electron R2: ' num2str(spin_system.rlx.nott_r2e) ' Hz']);
        report(spin_system,['Nuclear R1:  ' num2str(spin_system.rlx.nott_r1n) ' Hz']);
        report(spin_system,['Nuclear R2:  ' num2str(spin_system.rlx.nott_r2n) ' Hz']);
        report(spin_system,'========================================');
        
    case 'rlx_rates_weiz'
        report(spin_system,' ');
        report(spin_system,'==== Weizmann DNP relaxation theory ====');
        report(spin_system,['Electron R1: ' num2str(spin_system.rlx.weiz_r1e) ' Hz']);
        report(spin_system,['Electron R2: ' num2str(spin_system.rlx.weiz_r2e) ' Hz']);
        report(spin_system,['Nuclear R1:  ' num2str(spin_system.rlx.weiz_r1n) ' Hz']);
        report(spin_system,['Nuclear R2:  ' num2str(spin_system.rlx.weiz_r2n) ' Hz']);
        [rows,cols,vals]=find(spin_system.rlx.weiz_r1d);
        for n=1:numel(vals)
            report(spin_system,['Inter-nuclear dipolar R1(' num2str(rows(n)) ',' num2str(cols(n)) '): ' num2str(vals(n)) ' Hz']);
        end
        [rows,cols,vals]=find(spin_system.rlx.weiz_r2d);
        for n=1:numel(vals)
            report(spin_system,['Inter-nuclear dipolar R2(' num2str(rows(n)) ',' num2str(cols(n)) '): ' num2str(vals(n)) ' Hz']);
        end
        report(spin_system,'========================================');
        
    case 'symmetry'
        report(spin_system,header);
        report(spin_system,'=====================');
        report(spin_system,' Group    Spins      ');
        report(spin_system,'---------------------');
        for n=1:length(spin_system.comp.sym_spins)
            report(spin_system,['  ' spin_system.comp.sym_group{n} '     ' num2str(spin_system.comp.sym_spins{n})]);
        end
        report(spin_system,'=====================');
        
    case 'basis'
        nstates=size(spin_system.bas.basis,1);
        if nstates > spin_system.tols.basis_hush
            report(spin_system,['over ' num2str(spin_system.tols.basis_hush) ' states in the basis - printing suppressed.']);
        else
            report(spin_system,'final basis set summary (L,M quantum numbers in irreducible spherical tensor products).')
            report(spin_system,['N       ' blanks(length(num2str(spin_system.comp.nspins))) num2str(1:spin_system.comp.nspins,['%d ' blanks(7-length(num2str(spin_system.comp.nspins)))])]);
            for n=1:nstates
                current_line=blanks(7+8*spin_system.comp.nspins); spin_number=num2str(n);
                current_line(1:length(spin_number))=spin_number;
                for k=1:spin_system.comp.nspins
                    [L,M]=lin2lm(spin_system.bas.basis(n,k));
                    current_line(7+8*(k-1)+1)='(';
                    current_line(7+8*(k-1)+2)=num2str(L);
                    current_line(7+8*(k-1)+3)=',';
                    proj=num2str(M);
                    switch length(proj)
                        case 1
                            current_line(7+8*(k-1)+4)=proj;
                            current_line(7+8*(k-1)+5)=')';
                        case 2
                            current_line(7+8*(k-1)+4)=proj(1);
                            current_line(7+8*(k-1)+5)=proj(2);
                            current_line(7+8*(k-1)+6)=')';
                    end
                end
                report(spin_system,current_line);
            end
            report(spin_system,' ');
        end
        report(spin_system,['state space dimension ' num2str(nstates) ' (' num2str(100*nstates/(prod(spin_system.comp.mults)^2)) '% of the full state space).']);
        
    otherwise
        error('unknown topic.');
        
end

end

% 1: Blessed are the strong, for they shall possess the earth -- cursed are
% the weak, for they shall inherit the yoke.
%
% 2: Blessed are the powerful, for they shall be reverenced among men --
% cursed are the feeble, for they shall be blotted out.
%
% 3: Blessed are the bold, for they shall be masters of the world -- cursed
% are the humble, for they shall be trodden under hoofs.
%
% 4: Blessed are the victorious, for victory is the basis of right --
% cursed are the vanquished, for they shall be vassals forever.
%
% 5: Blessed are the iron-handed, for the unfit shall flee before them --
% cursed are the poor in spirit, for they shall be spat upon.
%
% 6: Blessed are the death-defiant, for their days shall be long in the
% lands -- cursed are the gazers toward a richer life beyond the grave, for
% they shall perish amidst plenty.
%
% 7: Blessed are the destroyers of false hope, for they are true Messiahs --
% cursed are the God-adorers, for they shall be shorn sheep.
%
% 8: Blessed are the valiant, for they shall obtain great treasure --
% cursed are the believers in good and evil, for they are frightened by
% shadows.
%
% 9: Blessed are those who believe in what is best for them, for never
% shall their minds be terrorized -- cursed are the "lambs of God", for
% they shall be bled whiter than snow.
%
% 10: Blessed is the man who has a sprinkling of enemies, for they shall
% make him a hero -- cursed is he who doeth good unto others who sneer upon
% him in return, for he shall be despised.
%
% 11: Blessed are the mighty-minded, for they shall ride the whirlwinds --
% cursed are they who teach lies for truth and truth for lies, for they are
% an abomination.
%
% 12: Thrice cursed are the weak whose insecurity makes them vile, for they
% shall serve and suffer.
%
% 13: The angel of self-deceit is camped in the souls of the "righteous" --
% the eternal flame of power through joy dwelleth within the flesh of a
% Satanist.
%
% Anton Szandor LaVey, "Satanic Bible"


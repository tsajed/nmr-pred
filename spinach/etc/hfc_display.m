% Ellipsoid plots of hyperfine coupling tensors.
%
% <http://spindynamics.org/wiki/index.php?title=Hfc_display.m>

function hfc_display(props,atoms,scaling_factor,conmatrix)

% Set up graphics
figure(); clf reset; hold on; colormap hot; opengl software;
light('Position',[-2,2,20]); light('Position',[10,10,10]);

% Draw the molecule
molplot(props.std_geom,conmatrix);

% Get a unit sphere
k=4; n=2^k-1; theta=pi*(-n:2:n)/n; phi=(pi/2)*(-n:2:n)'/n;
X=cos(phi)*cos(theta); Y=cos(phi)*sin(theta);
Z=sin(phi)*ones(size(theta));

% Loop over atoms
for n=1:props.natoms
    
    % Ignore atoms that are not on the user-supplied list
    if ismember(props.symbols{n},atoms)
        
        % Do not plot tensors smaller than 0.1 Gauss
        if norm(props.hfc.full.eigvals{n})>0.1
            
            % Stretch and scale sphere vertex coordinates
            coords=[reshape(X*props.hfc.full.eigvals{n}(1)*scaling_factor,1,length(phi)*length(theta));...
                    reshape(Y*props.hfc.full.eigvals{n}(2)*scaling_factor,1,length(phi)*length(theta));...
                    reshape(Z*props.hfc.full.eigvals{n}(3)*scaling_factor,1,length(phi)*length(theta))];
            
            % Rotate sphere vertex coordinates
            coords=props.hfc.full.eigvecs{n}*coords;
            
            % Fold back sphere vertex coordinates
            Xd=reshape(coords(1,:),length(phi),length(theta));
            Yd=reshape(coords(2,:),length(phi),length(theta));
            Zd=reshape(coords(3,:),length(phi),length(theta));
            
            % Color vertices by amplitude
            Cd=sqrt(Xd.^2+Yd.^2+Zd.^2);
            
            % Translate vertices
            Xd=Xd+props.std_geom(n,1);
            Yd=Yd+props.std_geom(n,2);
            Zd=Zd+props.std_geom(n,3);
            
            % Draw the surface
            surf(Xd,Yd,Zd,Cd,'FaceAlpha',0.5,'EdgeAlpha',0.1);
            
            % Scale eigenvectors
            vector_a=props.hfc.full.eigvecs{n}(:,1)*props.hfc.full.eigvals{n}(1)*scaling_factor;
            vector_b=props.hfc.full.eigvecs{n}(:,2)*props.hfc.full.eigvals{n}(2)*scaling_factor;
            vector_c=props.hfc.full.eigvecs{n}(:,3)*props.hfc.full.eigvals{n}(3)*scaling_factor;
            
            % Color eigenvectors
            col_a='r-'; col_b='r-'; col_c='r-';
            if props.hfc.full.eigvals{n}(1)<0, col_a='b-'; end
            if props.hfc.full.eigvals{n}(2)<0, col_b='b-'; end
            if props.hfc.full.eigvals{n}(3)<0, col_c='b-'; end
            
            % Draw eigenvectors
            plot3([props.std_geom(n,1)-vector_a(1) props.std_geom(n,1)+vector_a(1)],...
                  [props.std_geom(n,2)-vector_a(2) props.std_geom(n,2)+vector_a(2)],...
                  [props.std_geom(n,3)-vector_a(3) props.std_geom(n,3)+vector_a(3)],col_a,'LineWidth',1);
            plot3([props.std_geom(n,1)-vector_b(1) props.std_geom(n,1)+vector_b(1)],...
                  [props.std_geom(n,2)-vector_b(2) props.std_geom(n,2)+vector_b(2)],...
                  [props.std_geom(n,3)-vector_b(3) props.std_geom(n,3)+vector_b(3)],col_b,'LineWidth',1);
            plot3([props.std_geom(n,1)-vector_c(1) props.std_geom(n,1)+vector_c(1)],...
                  [props.std_geom(n,2)-vector_c(2) props.std_geom(n,2)+vector_c(2)],...
                  [props.std_geom(n,3)-vector_c(3) props.std_geom(n,3)+vector_c(3)],col_c,'LineWidth',1);
            
        end
    end
end

% Tidy up the picture
axis square; axis tight; axis equal;
set(gca,'Projection','perspective');
box on; grid on;

end

% The record for the shortest ever postdoc term in IK's group belongs
% to Gaby Slavcheva - 6 hours. She had turned up on a Monday morning, 
% went through a tour of the lab, a tour of the computing systems and
% an introduction to the SVN repositories that the group operates. IK
% explained that the repositories keep track of everyone's contributi-
% on to every project in a rigorous and objective way, making it easy
% to attribute credit, blame and to decide the author list of each pa-
% per. Gaby resigned her post that same evening, citing "draconian re-
% strictions" on her "academic freedom", and claiming two months' pay
% for "redundancy" from Southampton's hapless HR service.
%
% IK did not object. Last thing he heard many years down the line was
% that, with over 100 peer-reviewed papers and at the age of 56, Gaby
% was a postdoc at the University of Bath.

 
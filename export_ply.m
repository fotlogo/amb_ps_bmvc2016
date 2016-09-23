%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Description..... : 	Save the reconstruction PLY format, in 
%						mesh.ply, mesh.mtl and mesh.png
%						INPUT : XYZ, mask
%						
%	Author ......... : 	Fotios Logothetis (adapted from Yvain Queau) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function export_ply(XYZ,mask,filename)
	
	if (~exist('mask','var')|isempty(mask)) mask=ones(size(XYZ(:,:,1))); end;
	
	[nrows,ncols] = size(mask);
    %% SWITCH TO OPENGL AXIS
    XYZ(:,:,2) = - XYZ(:,:,2);
	XYZ(:,:,3) = - XYZ(:,:,3);
	%%	
	% Nuage de Points : 
	indices_mask = find(mask>0);
	[Imask,Jmask]=ind2sub(size(mask),indices_mask);
	indices = zeros(size(mask));
	indices(indices_mask) = 1:length(indices_mask);		
	mask=[mask;zeros(1,size(mask,2))];
	mask=[mask,zeros(size(mask,1),1)];
	X = XYZ(:,:,1);
	Y = XYZ(:,:,2);
	Z = XYZ(:,:,3);
	vertices = [X(indices_mask),Y(indices_mask),Z(indices_mask)];	
	
	disp('Meshing ...')
	indices_lower_triangle = find(mask(1:end-1,1:end-1)>0 & mask(2:end,2:end) & mask(2:end,1:end-1)); 
	[I_lt,J_lt] = ind2sub([nrows ncols],indices_lower_triangle);	
	indices_bas = sub2ind([nrows ncols],I_lt+1,J_lt);
	indices_bas_droite = sub2ind([nrows ncols],I_lt+1,J_lt+1);
	face_vertices = [indices(indices_lower_triangle),indices(indices_bas),indices(indices_bas_droite)];
	
	indices_upper_triangle = find(mask(1:end-1,1:end-1)>0 & mask(2:end,2:end) & mask(1:end-1,2:end));
	[I_ut,J_ut] = ind2sub([nrows ncols],indices_upper_triangle); 
	indices_droite = sub2ind([nrows ncols],I_ut,J_ut+1);
	indices_bas_droite = sub2ind([nrows ncols],I_ut+1,J_ut+1);
	face_vertices = [face_vertices;...
					 indices(indices_upper_triangle),indices(indices_bas_droite),indices(indices_droite)];	
    
    face_vertices=face_vertices-1;
                 
    [nverts,~]=size(vertices);
    [nfaces,~]=size(face_vertices);
    
    fprintf(1, 'writing %s ... \n',filename);
                 
	fileID = fopen(filename,'w');
    
    fprintf(fileID,'ply\nformat ascii 1.0\nelement vertex %d\nproperty float x\nproperty float y\nproperty float z\n',nverts);
    fprintf(fileID,'element face %d\nproperty list uint8 int32 vertex_indices\nend_header\n', nfaces);       
    
    for ii=1:nverts
        fprintf(fileID,'%f %f %f\n',vertices(ii,1),vertices(ii,2),vertices(ii,3));
    end
    
     for ii=1:nfaces
        fprintf(fileID,'3 %d %d %d\n',face_vertices(ii,1),face_vertices(ii,2),face_vertices(ii,3));
    end
   
    fclose(fileID); 
	

end



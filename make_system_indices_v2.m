function [I, J, indices_centre, I_xy, I_x1_y, I_x_y1,Gx,Gy ] = make_system_indices_v2(mask)
%% concides both fw and bw differences so it has 4 different omegas
  [nrows, ncols]=size(mask);
  
  Omega = zeros(size(mask,1),size(mask,2),4);
  Omega_padded = padarray(mask,[1 1],0); %add one row and one col of zeros in both side (size +2)
  
  order=1;
	
	if(order>1)
		Omega(:,:,1) = Omega_padded(3:end,2:end-1).*Omega_padded(2:end-1,3:end).*mask;
		Omega(:,:,2) = Omega_padded(3:end,2:end-1).*Omega_padded(2:end-1,1:end-2).*mask;
		Omega(:,:,3) = Omega_padded(1:end-2,2:end-1).*Omega_padded(2:end-1,1:end-2).*mask;
		Omega(:,:,4) = Omega_padded(1:end-2,2:end-1).*Omega_padded(2:end-1,3:end).*mask;
	elseif(order==1)
        NOmega=zeros(nrows,ncols,3,'logical');
        %fw differences in both dimensions
		Omega(:,:,1) = Omega_padded(3:end,2:end-1).*Omega_padded(2:end-1,3:end).*mask;
        NOmega(:,:,1) = ~Omega(:,:,1);
        %fw in y, bw in x 
        Omega(:,:,2) = NOmega(:,:,1).*Omega_padded(3:end,2:end-1).*Omega_padded(2:end-1,1:end-2).*mask;
        NOmega(:,:,2) = ~Omega(:,:,2);
        %both bw
        Omega(:,:,3) = NOmega(:,:,1).*NOmega(:,:,2).*Omega_padded(1:end-2,2:end-1).*Omega_padded(2:end-1,1:end-2).*mask;
        NOmega(:,:,3) = ~Omega(:,:,3);
        %bw in y, fw in x 
        Omega(:,:,4) = NOmega(:,:,1).*NOmega(:,:,2).*NOmega(:,:,3).*Omega_padded(1:end-2,2:end-1).*Omega_padded(2:end-1,3:end).*mask;
    end
	clear Omega_padded  NOmega
% 	Omega_rep = repmat(Omega,[1 1 1 m]);		
	indices_mask = (mask>0);

	nb_pixels_mask = sum(mask(:));%numel(indices_mask)
	mapping_matrix = zeros([nrows ncols]);
	mapping_matrix(indices_mask)=1:nb_pixels_mask; 
    
    % Regularization 	
	[X,Y]=find(mask>0);
	indices_centre=sub2ind(size(mask),X,Y);
	I_centre = mapping_matrix(indices_centre);	
    
    I = [];
	J = [];
    
	I=[I;I_centre(:)]; % Which lines are we talking about ? 
	J=[J;I_centre(:)]; % Which columns ?
	
	clear I_centre indices_centre
	
    indices_centre=cell(4,1);
    I_xy=cell(4,1);
    I_x1_y=cell(4,1);
    I_x_y1=cell(4,1);    
    
	% Terms in Omega(:,:,1)
	[Ii,Ji,indices_centre1,I_xy_1,I_voisins_x1,I_voisins_y1]=indices_omega(Omega(:,:,1),mapping_matrix,1,1);    
    I=[I;Ii];
    J=[J;Ji];
    indices_centre{1}=indices_centre1;
    I_xy{1}=I_xy_1;
    I_x1_y{1}=I_voisins_x1;
    I_x_y1{1}=I_voisins_y1;    
    % Terms in Omega(:,:,2)
	[Ii,Ji,indices_centre2,I_centre2,I_voisins_x2,I_voisins_y2]=indices_omega(Omega(:,:,2),mapping_matrix,1,-1);
    I=[I;Ii];
    J=[J;Ji];
    indices_centre{2}=indices_centre2;
    I_xy{2}=I_centre2;
    I_x1_y{2}=I_voisins_x2;
    I_x_y1{2}=I_voisins_y2;    
    % Terms in Omega(:,:,3)
	[Ii,Ji,indices_centre3,I_centre3,I_voisins_x3,I_voisins_y3]=indices_omega(Omega(:,:,3),mapping_matrix,-1,-1);
    I=[I;Ii];
    J=[J;Ji];
    indices_centre{3}=indices_centre3;
    I_xy{3}=I_centre3;
    I_x1_y{3}=I_voisins_x3;
    I_x_y1{3}=I_voisins_y3;
    % Terms in Omega(:,:,4)
    [Ii,Ji,indices_centre4,I_centre4,I_voisins_x4,I_voisins_y4]=indices_omega(Omega(:,:,4),mapping_matrix,-1,1);    
    I=[I;Ii];
    J=[J;Ji];
    indices_centre{4}=indices_centre4;  
    I_xy{4}=I_centre4;
    I_x1_y{4}=I_voisins_x4;
    I_x_y1{4}=I_voisins_y4;
    
    %% MAKE grand matrices, i.e. zx=Gx*z (z,zx mask vectors)
    %to know when are we doing fw or backwards differences
    b1f=[-1,-1,+1,+1];
    b2f=[-1,+1,+1,-1];  
    NN=size(I_xy{1},1)+size(I_xy{2},1)+size(I_xy{3},1)+size(I_xy{4},1);
    II=zeros(2*NN,1); %same for both x and y
    JJx=zeros(2*NN,1); KKx=zeros(2*NN,1); JJy=zeros(2*NN,1); KKy=zeros(2*NN,1);
    %
    curr_indx_2=1;
    for indx=1:4
        N_omega_current=size(I_xy{indx},1);
        II(curr_indx_2:curr_indx_2+N_omega_current-1)=I_xy{indx};
        %
        JJx(curr_indx_2:curr_indx_2+N_omega_current-1)=I_x1_y{indx}; 
        KKx(curr_indx_2:curr_indx_2+N_omega_current-1)=-b1f(indx)*ones(N_omega_current,1);  
        JJy(curr_indx_2:curr_indx_2+N_omega_current-1)=I_x_y1{indx}; 
        KKy(curr_indx_2:curr_indx_2+N_omega_current-1)=-b2f(indx)*ones(N_omega_current,1);   
        curr_indx_2=curr_indx_2+N_omega_current;
        
        II(curr_indx_2:curr_indx_2+N_omega_current-1)=I_xy{indx};
        %
        JJx(curr_indx_2:curr_indx_2+N_omega_current-1)=I_xy{indx}; 
        KKx(curr_indx_2:curr_indx_2+N_omega_current-1)=b1f(indx)*ones(N_omega_current,1);   
        JJy(curr_indx_2:curr_indx_2+N_omega_current-1)=I_xy{indx}; 
        KKy(curr_indx_2:curr_indx_2+N_omega_current-1)=b2f(indx)*ones(N_omega_current,1);    
        curr_indx_2=curr_indx_2+N_omega_current;
    end
    %
    Gx=sparse(II,JJx,KKx,nb_pixels_mask,nb_pixels_mask);
    Gy=sparse(II,JJy,KKy,nb_pixels_mask,nb_pixels_mask);
   
  function [I,J,indices_centre1,I_xy,I_x1_y,I_x_y1]=indices_omega(Omega,mapping_matrix,sign_x,sign_y)
        [xx,yy]=find(Omega>0);
        indices_centre1=sub2ind(size(Omega),xx,yy);
        I_xy = mapping_matrix(indices_centre1);
        indices_voisins_x1=sub2ind(size(Omega),xx+sign_x,yy);
        I_x1_y = mapping_matrix(indices_voisins_x1);
        indices_voisins_y1=sub2ind(size(Omega),xx,yy+sign_y);
        I_x_y1 = mapping_matrix(indices_voisins_y1);  
        %%
	nb_pixels_omega=length(I_xy);
    nb_pixels_I=9*nb_pixels_omega;   
    
    I=zeros(nb_pixels_I,1);
    J=zeros(nb_pixels_I,1);

    curr_indx=1;
    %1
    I(curr_indx:curr_indx+nb_pixels_omega-1)=I_xy(:);
    J(curr_indx:curr_indx+nb_pixels_omega-1)=I_xy(:);    
    curr_indx=curr_indx+nb_pixels_omega;      
    %2
    I(curr_indx:curr_indx+nb_pixels_omega-1)=I_xy(:);
    J(curr_indx:curr_indx+nb_pixels_omega-1)=I_x1_y(:);    
    curr_indx=curr_indx+nb_pixels_omega;       
    %3
    I(curr_indx:curr_indx+nb_pixels_omega-1)=I_xy(:);
    J(curr_indx:curr_indx+nb_pixels_omega-1)=I_x_y1(:);    
    curr_indx=curr_indx+nb_pixels_omega;        
    %4
    I(curr_indx:curr_indx+nb_pixels_omega-1)=I_x1_y(:);
    J(curr_indx:curr_indx+nb_pixels_omega-1)=I_xy(:);    
    curr_indx=curr_indx+nb_pixels_omega;        
    %5
    I(curr_indx:curr_indx+nb_pixels_omega-1)=I_x1_y(:);
    J(curr_indx:curr_indx+nb_pixels_omega-1)=I_x1_y(:);    
    curr_indx=curr_indx+nb_pixels_omega;    
    %6
    I(curr_indx:curr_indx+nb_pixels_omega-1)=I_x1_y(:);
    J(curr_indx:curr_indx+nb_pixels_omega-1)=I_x_y1(:);    
    curr_indx=curr_indx+nb_pixels_omega;         
    %7
    I(curr_indx:curr_indx+nb_pixels_omega-1)=I_x_y1(:);
    J(curr_indx:curr_indx+nb_pixels_omega-1)=I_xy(:);    
    curr_indx=curr_indx+nb_pixels_omega;           
    %8
    I(curr_indx:curr_indx+nb_pixels_omega-1)=I_x_y1(:);
    J(curr_indx:curr_indx+nb_pixels_omega-1)=I_x1_y(:);    
    curr_indx=curr_indx+nb_pixels_omega;   
    %9
    I(curr_indx:curr_indx+nb_pixels_omega-1)=I_x_y1(:);
    J(curr_indx:curr_indx+nb_pixels_omega-1)=I_x_y1(:); 
    end

end


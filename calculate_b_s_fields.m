function [b,s]=calculate_b_s_fields(I,mask_indx,x,y,Z,f,C,A,H,thresholds, ambient)
%calculate_b_s_fields Under perspective projection with or without ambient
%chose one of the 2 versions 
if ambient
    [b,s]=calculate_b_s_fields_amb(I,mask_indx,x,y,Z,f,C,A,H,thresholds);
else
    [b,s]=calculate_b_s_fields_dark(I,mask_indx,x,y,Z,f,C,A,H,thresholds);
end
    %% PUT a flat patch in places where there are no data    
    ss=sum(s,3);
    holes=(ss==0);
    b1=b(:,:,1,1);
    b1(holes)=1;
    b(:,:,1,1)=b1;
end
function [b,s]=calculate_b_s_fields_dark(I,indices_mask,x,y,Z,f,C,A,H,thresholds)
shadow_threshold=thresholds(1);
saturation_thresh=thresholds(2);
% Fundamental equation fields
 [nrows,ncols,nb_images] = size(I);
% [h1,h2,h3,h4]=size(H)
nb_combinations = nchoosek(nb_images,2);
b = zeros(nrows,ncols,nb_combinations,2);
s = zeros(nrows,ncols,nb_combinations);

b1=zeros(nrows,ncols);
b2=zeros(nrows,ncols);
s12=zeros(nrows,ncols);

current_index = 1;
    for k = 1:nb_images-1
        for i = k+1:nb_images 
            %unfortunatelly, matlab does not allow I(mask,:) so we need
            %temp varibles
            A_i=A(:,:,i);
            A_k=A(:,:,k);
            I_i = I(:,:,i);
            I_k = I(:,:,k);
            H_1_i=H(:,:,1,i);
            H_2_i=H(:,:,2,i);
            H_3_i=H(:,:,3,i);
            
            H_1_k=H(:,:,1,k);
            H_2_k=H(:,:,2,k);
            H_3_k=H(:,:,3,k);                   
            
           A_i_I_k_c= (A_i(indices_mask).*I_k(indices_mask)).^C(indices_mask);
           A_k_I_i_c= (A_k(indices_mask).*I_i(indices_mask)).^C(indices_mask);  
           %% FIXME make sure the order of b1, b2 is correct
           %it looks like it should be reversed cause of matlab I(y,x)
           %convention
           b2(indices_mask) = A_i_I_k_c.*(H_1_i(indices_mask)-x(indices_mask).*H_3_i(indices_mask)/f) - ...
                               A_k_I_i_c.*(H_1_k(indices_mask)-x(indices_mask).*H_3_k(indices_mask)/f);
           b1(indices_mask) = A_i_I_k_c.*(H_2_i(indices_mask)-y(indices_mask).*H_3_i(indices_mask)/f)-...
                               A_k_I_i_c.*(H_2_k(indices_mask)-y(indices_mask).*H_3_i(indices_mask)/f);
           s12(indices_mask)= (1+Z(indices_mask)/f).*(A_i_I_k_c.*H_3_i(indices_mask)-A_k_I_i_c.*H_3_k(indices_mask));          
          %TODO instead of the 0=0 entries, reduce the size of the matrices
            b1(I_k<shadow_threshold) = 0;
            b2(I_k<shadow_threshold) = 0;
            s12(I_k<shadow_threshold) = 0;
            
            b1(I_i<shadow_threshold) = 0;
            b2(I_i<shadow_threshold) = 0;
            s12(I_i<shadow_threshold) = 0;
          %% SATURATION Threshold       
%           sum(sum((I_k>0.45)))
            b1(I_k>saturation_thresh) = 0;
            b2(I_k>saturation_thresh) = 0;
            s12(I_k>saturation_thresh) = 0;
            
            b1(I_i>saturation_thresh) = 0;
            b2(I_i>saturation_thresh) = 0;
            s12(I_i>saturation_thresh) = 0;            
            
            b(:,:,current_index,1) = b1;
            b(:,:,current_index,2) = b2;
            s(:,:,current_index) = s12;
            current_index = current_index+1;
        end
    end   
end
function [b,s]=calculate_b_s_fields_amb(I,mask_indx,x,y,Z,f,C,A,H,thresholds)

shadow_threshold=thresholds(1);
saturation_thresh=thresholds(2);

[nrows,ncols,nb_images] = size(I);

nb_combinations = nchoosek(nb_images,4); %%IF A LOT OF IMAGES THIS IS HUGE

b = zeros(nrows,ncols,nb_combinations,2);
s = zeros(nrows,ncols,nb_combinations);

b1=zeros(nrows,ncols);
b2=zeros(nrows,ncols);
s12=zeros(nrows,ncols);

current_index = 1;

for ii = 1 : nb_images-3 
        %unfortunatelly, matlab does not allow I(mask,:) so we need
        %temp varibles 
        I_i=I(:,:,ii);
        
        A_i=A(:,:,ii);        
        H_1_i=H(:,:,1,ii);
        H_1_i=H_1_i(mask_indx);
        H_2_i=H(:,:,2,ii);
        H_2_i=H_2_i(mask_indx);
        H_3_i=H(:,:,3,ii);
        H_3_i=H_3_i(mask_indx);
        
        gamma_i= ( A_i(mask_indx).^C(mask_indx) ) ./  ( I_i(mask_indx).^( C(mask_indx)-1)) ;  
        
        for jj = ii+1 : nb_images-2           
            I_j=I(:,:,jj);
            
            A_j=A(:,:,jj);            
            H_1_j=H(:,:,1,jj);
            H_1_j=H_1_j(mask_indx);
            H_2_j=H(:,:,2,jj);
            H_2_j=H_2_j(mask_indx);
            H_3_j=H(:,:,3,jj);
            H_3_j=H_3_j(mask_indx);
            
            gamma_j= ( A_j(mask_indx).^C(mask_indx) ) ./  ( I_j(mask_indx).^( C(mask_indx)-1)) ;
                     
            dI_i_j= I_i(mask_indx) - I_j(mask_indx);  %naming convention: d:difference, I image, j and i  
            
            for q = jj+1 : nb_images-1   ;          
                I_q=I(:,:,q);
                
                A_q=A(:,:,q);
                H_1_q=H(:,:,1,q);
                H_1_q=H_1_q(mask_indx);
                H_2_q=H(:,:,2,q);
                H_2_q=H_2_q(mask_indx);
                H_3_q=H(:,:,3,q);
                H_3_q=H_3_q(mask_indx);
                
                gamma_q= ( A_q(mask_indx).^C(mask_indx) ) ./  ( I_q(mask_indx).^( C(mask_indx)-1)) ;                

                for r = q+1 : nb_images
                    I_r=I(:,:,r);                    
                    
                    A_r=A(:,:,r);
                    H_1_r=H(:,:,1,r);
                    H_1_r=H_1_r(mask_indx);
                    H_2_r=H(:,:,2,r);
                    H_2_r=H_2_r(mask_indx);
                    H_3_r=H(:,:,3,r);
                    H_3_r=H_3_r(mask_indx);
                    
                    gamma_r= ( A_r(mask_indx).^C(mask_indx) ) ./  ( I_r(mask_indx).^( C(mask_indx)-1)) ;
                    
                    dI_q_r= I_q(mask_indx) - I_r(mask_indx);   
               %%set b1 (actually b2)      
               b2(mask_indx) =      dI_q_r.* ( gamma_i.* (H_1_i- x(mask_indx).*H_3_i/f) - gamma_j.* (H_1_j- x(mask_indx).*H_3_j/f)) - ...
                           dI_i_j.* ( gamma_q.* (H_1_q- x(mask_indx).*H_3_q/f) - gamma_r.* (H_1_r- x(mask_indx).*H_3_r/f));
               %%set b2 (actually b1)      
               b1(mask_indx) =      dI_q_r.* ( gamma_i.* (H_2_i- y(mask_indx).*H_3_i/f) - gamma_j.* (H_2_j- y(mask_indx).*H_3_j/f)) - ...
                           dI_i_j.* ( gamma_q.* (H_2_q- y(mask_indx).*H_3_q/f) - gamma_r.* (H_2_r- y(mask_indx).*H_3_r/f));          
              %set s         
               s12(mask_indx)= (1+Z(mask_indx)./f).* ( dI_q_r.*(gamma_i.*H_3_i - gamma_j.*H_3_j) - dI_i_j.*(gamma_q.*H_3_q - gamma_r.*H_3_r));                    
            
            %% SHADOW/LOW DIFF Threshold    
            %TODO instead of the 0=0 entries, reduce the size of the matrices
       
          
            b1( abs(I_i-I_j)<shadow_threshold) = 0;
            b2( abs(I_i-I_j)<shadow_threshold) = 0;
            s12( abs(I_i-I_j)<shadow_threshold) = 0;      
                    
            b1(abs(I_q-I_r)<shadow_threshold) = 0;
            b2(abs(I_q-I_r)<shadow_threshold) = 0;
            s12(abs(I_q-I_r)<shadow_threshold) = 0;%  
          
          %% SATURATION Threshold  
                   b1(I_i>saturation_thresh) = 0;
                   b2(I_i>saturation_thresh) = 0;
                   s12(I_i>saturation_thresh) = 0;
% %             
                   b1(I_j>saturation_thresh) = 0;
                   b2(I_j>saturation_thresh) = 0;
                   s12(I_j>saturation_thresh) = 0;%    
                   
                   b1(I_q>saturation_thresh) = 0;
                   b2(I_q>saturation_thresh) = 0;
                   s12(I_q>saturation_thresh) = 0;
%                    
                   b1(I_r>saturation_thresh) = 0;
                   b2(I_r>saturation_thresh) = 0;
                   s12(I_r>saturation_thresh) = 0;        
          
            b(:,:,current_index,1) = b1;
            b(:,:,current_index,2) = b2;
            s(:,:,current_index) = s12;
            current_index = current_index+1;
                end
            end
        end
end

end 

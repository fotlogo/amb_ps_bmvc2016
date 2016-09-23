function [b,s]=calculate_b_s_fields(I,indices_mask,x,y,Z,f,C,A,H,thresholds,ambient)
%calculate_b_s_fields Under perspective projection with or without ambient
%chose one of the 2 versions 
if ambient==0
    [b,s]=calculate_b_s_fields_simple(I,indices_mask,x,y,Z,f,C,A,H,thresholds);
else
%       [bv1,sv1]=calculate_b_s_fields_ambient(I,indices_mask,x,y,Z,f,C,A,H,thresholds);
%       [bv2,sv2]=calculate_b_s_fields_ambient_v2(I,indices_mask,x,y,Z,f,C,A,H,thresholds);
%       
%       b=cat(3,bv1,0.1*bv2);
%       s=cat(3,sv1,0.1*sv2);
 [b,s]=calculate_b_s_fields_ambient(I,indices_mask,x,y,Z,f,C,A,H,thresholds);
end
    %% PUT a flat patch in places where there are no data    
ss=sum(s,3);
holes=(ss==0);
b1=b(:,:,1,1);
b1(holes)=1;
b(:,:,1,1)=b1;
end
function [b,s]=calculate_b_s_fields_simple(I,indices_mask,x,y,Z,f,C,A,H,thresholds)

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

cm=C(indices_mask);

current_index = 1;
    for ii = 1:nb_images-1 
         %unfortunatelly, matlab does not allow I(mask,:) so we need
         %temp varibles
         I_i = I(:,:,ii);
         A_i=A(:,:,ii);
         
         H_1_i=H(:,:,1,ii);
         H_2_i=H(:,:,2,ii);
         H_3_i=H(:,:,3,ii);  
            
        for jj = ii+1:nb_images 
           I_j = I(:,:,jj);
           A_j=A(:,:,jj);         
           
           H_1_j=H(:,:,1,jj);
           H_2_j=H(:,:,2,jj);
           H_3_j=H(:,:,3,jj);                           
            
           A_j_I_i_c= (A_j(indices_mask).*I_i(indices_mask)).^cm;
           A_i_I_j_c= (A_i(indices_mask).*I_j(indices_mask)).^cm;  
           %% FIXME make sure the order of b1, b2 is correct
           %it looks like it should be reversed cause of matlab I(y,x)
           %convention
           b2(indices_mask) = A_j_I_i_c.*(H_1_j(indices_mask)-x(indices_mask).*H_3_j(indices_mask)/f) - ...
                               A_i_I_j_c.*(H_1_i(indices_mask)-x(indices_mask).*H_3_i(indices_mask)/f);
           b1(indices_mask) = A_j_I_i_c.*(H_2_j(indices_mask)-y(indices_mask).*H_3_j(indices_mask)/f)-...
                               A_i_I_j_c.*(H_2_i(indices_mask)-y(indices_mask).*H_3_j(indices_mask)/f);
           s12(indices_mask)= (1+Z(indices_mask)/f).*(A_j_I_i_c.*H_3_j(indices_mask)-A_i_I_j_c.*H_3_i(indices_mask));          
          %TODO instead of the 0=0 entries, reduce the size of the matrices
            b1(I_i<shadow_threshold) = 0;
            b2(I_i<shadow_threshold) = 0;
            s12(I_i<shadow_threshold) = 0;
            
            b1(I_j<shadow_threshold) = 0;
            b2(I_j<shadow_threshold) = 0;
            s12(I_j<shadow_threshold) = 0;             
%              [XX YY ZZ] = meshgrid(1:size(x,1), 1:size(x,2), 1:2);
             
%             figure;
% %             streamslice(XX,YY,ZZ,b2,b1,s12,[],[],[]);
%             quiver(x,y,b2,b1);
%             pause
            
            
          %% SATURATION Threshold       
%           sum(sum((I_k>0.45)))
            b1(I_i>saturation_thresh) = 0;
            b2(I_i>saturation_thresh) = 0;
            s12(I_i>saturation_thresh) = 0;
            
            b1(I_j>saturation_thresh) = 0;
            b2(I_j>saturation_thresh) = 0;
            s12(I_j>saturation_thresh) = 0;      
            
%             ss=(I_j>saturation_thresh | I_i>saturation_thresh) ;
%             ns=sum(sum(ss));
%             nsp=100*ns/size(A_j_I_i_c,1);
%             fprintf(1,'removed %d pixels (%.2f %%)due to saturation (>%.2f)\n', ns, nsp,saturation_thresh);     
            
            b(:,:,current_index,1) = b1;
            b(:,:,current_index,2) = b2;
            s(:,:,current_index) = s12;
            current_index = current_index+1;
        end
    end      
end 
%get the index of the next LED in the outer circle
function jj=out_circle_next(ii,n)
    jj=ii+n;
    if(jj==16)
        jj=24;
    end
    if(jj==25)
        jj=1;
    end
end
%get the index of the next LED in the inner circle
function r=in_circle_next(q,n)
    r=q+n;
    if(r==24)
        r=16;
    end
end

function [b,s]=calculate_b_s_fields_ambient(I,mask_indx,x,y,Z,f,C,A,H,thresholds)

shadow_threshold=thresholds(1);
saturation_thresh=thresholds(2);
%% ENUMERATING CONSTS
COMBINATIONS_24=0;
COMBINATIONS_3=1;
COMBINATIONS_FULL=2;
% Fundamental equation fields
 [nrows,ncols,nb_images] = size(I);
%If 24 images it is our 
if nb_images==24
    combinations_case= COMBINATIONS_24; 
    nb_combinations = 8*4;
elseif nb_images==3
    combinations_case= COMBINATIONS_3;
    nb_combinations =2; 
else
    combinations_case= COMBINATIONS_FULL; 
    nb_combinations = nchoosek(nb_images,4);
end

% combinations_case= COMBINATIONS_FULL; 
%     nb_combinations = nchoosek(nb_images,4);


% if combinations_case
%     nb_combinations = nchoosek(nb_images,4);
% else
%     nb_combinations = 8*4;
% end
b = zeros(nrows,ncols,nb_combinations,2);
s = zeros(nrows,ncols,nb_combinations);

b1=zeros(nrows,ncols);
b2=zeros(nrows,ncols);
s12=zeros(nrows,ncols);

current_index = 1;

if combinations_case== COMBINATIONS_24
	combinations_i= 1:8;
elseif combinations_case== COMBINATIONS_3
	combinations_i= 1;
else % combinations_case= COMBINATIONS_FULL; 
	combinations_i= 1 : nb_images-3 ;
end

for ii = combinations_i
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

        if combinations_case== COMBINATIONS_24           
            combinations_j=out_circle_next(ii,8);           
        elseif combinations_case== COMBINATIONS_3
            combinations_j= 2;
        else % combinations_case= COMBINATIONS_FULL; 
           combinations_j= ii+1 : nb_images-2;
        end

        for jj = combinations_j           
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

            if combinations_case== COMBINATIONS_24
                combinations_q= 16:19 ;
            elseif combinations_case== COMBINATIONS_3
                combinations_q= [1,2];
            else % combinations_case= COMBINATIONS_FULL; 
                combinations_q= jj+1 : nb_images-1 ;
            end
            
            for q = combinations_q   ;          
                I_q=I(:,:,q);
                
                A_q=A(:,:,q);
                H_1_q=H(:,:,1,q);
                H_1_q=H_1_q(mask_indx);
                H_2_q=H(:,:,2,q);
                H_2_q=H_2_q(mask_indx);
                H_3_q=H(:,:,3,q);
                H_3_q=H_3_q(mask_indx);
                
                gamma_q= ( A_q(mask_indx).^C(mask_indx) ) ./  ( I_q(mask_indx).^( C(mask_indx)-1)) ;
                
                if combinations_case== COMBINATIONS_24                   
                   combinations_r=in_circle_next(q,4);                  
                elseif combinations_case== COMBINATIONS_3
                    combinations_r= 3;
                else % combinations_case= COMBINATIONS_FULL; 
                     combinations_r= q+1 : nb_images;
                end
                for r = combinations_r
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
                  
             if(nb_images==4)  
                    m1= mean(abs(b1(:)));
                    m2= mean(abs(b2(:)));
                    m3= mean(abs(s12(:)));
                    
                    fprintf(1,'drawing isocontours: mean b1 %f mean b2 %f mean s (i.e. zero) %f \n',m1,m2,m3);%               
% %               figure;
% %               imshow(b1./m1);
% %                figure;
% %               imshow(b2./m2);   

                step=20;
                mx=min(x(:));
                my=min(y(:));
                Mx=max(x(:));
                My=max(y(:));      
        
              
                figure            
                hold on;
                starty =mx:step:Mx; %cover all y range
                 %rotate in gimp with shit+F
                for kk=my:step:My                        
            %                startx = zeros(1,ncols);%1:min(ncols,nrows);
                    startx = kk*ones(1,size(starty,2));%1:min(ncols,nrows);
                    my_streamline(x,y, b2,b1,startx,starty);  
%                   my_streamline(x,y,b1,flipud(b2),startx,starty);  
           
                end
             end   
%             quiver(x,y,b2,b1);   
%                 my_streamline(x,y,b2,b1,startx,starty); 
            %% SHADOW Threshold    
            %TODO instead of the 0=0 entries, reduce the size of the matrices
            %TODO shadow thresshold meaningless if ambient exists
            if combinations_case~= COMBINATIONS_3
                    b1( abs(I_i-I_j)<shadow_threshold) = 0;
                    b2( abs(I_i-I_j)<shadow_threshold) = 0;
                    s12( abs(I_i-I_j)<shadow_threshold) = 0;   
%                 nr= sum(abs(dI_q_r)<shadow_threshold | abs(dI_i_j)<shadow_threshold);
%                 nrp=100*nr/size(dI_q_r,1);
%                 fprintf(1,'removed %d pixels (%.2f %%)due to low dif (<%.2f)\n', nr, nrp,shadow_threshold);            
%             
%                     b1(I_j<shadow_threshold) = 0;
%                     b2(I_j<shadow_threshold) = 0;
%                     s12(I_j<shadow_threshold) = 0;%     
%                     
%                     b1(I_q<shadow_threshold) = 0;
%                     b2(I_q<shadow_threshold) = 0;
%                     s12(I_q<shadow_threshold) = 0;%     
                    
                    b1(abs(I_q-I_r)<shadow_threshold) = 0;
                    b2(abs(I_q-I_r)<shadow_threshold) = 0;
                    s12(abs(I_q-I_r)<shadow_threshold) = 0;%  
            end
          %% SATURATION Threshold       
% %           sum(sum((I_k>0.45)))
                if combinations_case~= COMBINATIONS_3
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
                   
%                    ss=(I_i>saturation_thresh | I_j>saturation_thresh | I_q>saturation_thresh | I_r>saturation_thresh ) ;
%                    ns=sum(sum(ss));
%                    nsp=100*ns/size(dI_q_r,1);
%                    fprintf(1,'removed %d pixels (%.2f %%)due to saturation (>%.2f)\n', ns, nsp,saturation_thresh);            
               end
            b(:,:,current_index,1) = b1;
            b(:,:,current_index,2) = b2;
            s(:,:,current_index) = s12;
            current_index = current_index+1;
                end
            end
        end
end 
% figure;
% % surf(ss);
% imshow(holes);
% pause

end 
function [b,s]=calculate_b_s_fields_ambient_v2(I,mask_indx,x,y,Z,f,C,A,H,thresholds)

shadow_threshold=thresholds(1);
saturation_thresh=thresholds(2);

 [nrows,ncols,nb_images] = size(I);
 assert(nb_images==24);
%If 24 images it is our 
nb_combinations = 8*4;

% if combinations_case
%     nb_combinations = nchoosek(nb_images,4);
% else
%     nb_combinations = 8*4;
% end
b = zeros(nrows,ncols,nb_combinations,2);
s = zeros(nrows,ncols,nb_combinations);

b1=zeros(nrows,ncols);
b2=zeros(nrows,ncols);
s12=zeros(nrows,ncols);

current_index = 1;

for ii = 1:8
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
%% J
        jj = out_circle_next(ii,5);
        
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
       
       D_i_j=2./( ( C(mask_indx)-1).*( 1./I_i(mask_indx)-1./I_j(mask_indx)));
 %% K      
        k = out_circle_next(jj,5);
%         for k = 16:19   ;         
       I_k=I(:,:,k);
            
       A_k=A(:,:,k);            
       H_1_k=H(:,:,1,k);
       H_1_k=H_1_k(mask_indx);
       H_2_k=H(:,:,2,k);
       H_2_k=H_2_k(mask_indx);
       H_3_k=H(:,:,3,k);
       H_3_k=H_3_k(mask_indx);
            
        gamma_k= ( A_k(mask_indx).^C(mask_indx) ) ./  ( I_k(mask_indx).^( C(mask_indx)-1)) ;
                     
        dI_i_k= I_i(mask_indx) - I_k(mask_indx); 
       
        D_i_k=2./( ( C(mask_indx)-1).*( 1./I_i(mask_indx)-1./I_k(mask_indx)));
   %% i-j-k        
        B_i_j_k=D_i_j.*dI_i_j-D_i_k.*dI_i_k; 

        for p = 16:19   ; 
            I_p=I(:,:,p);
                
            A_p=A(:,:,p);
            H_1_p=H(:,:,1,p);
            H_1_p=H_1_p(mask_indx);
            H_2_p=H(:,:,2,p);
            H_2_p=H_2_p(mask_indx);
            H_3_p=H(:,:,3,p);
            H_3_p=H_3_p(mask_indx);
                
            gamma_p= ( A_p(mask_indx).^C(mask_indx) ) ./  ( I_p(mask_indx).^( C(mask_indx)-1)) ;
%% Q
            q = in_circle_next(p,3);
            I_q=I(:,:,q);
                
            A_q=A(:,:,q);
            H_1_q=H(:,:,1,q);
            H_1_q=H_1_q(mask_indx);
            H_2_q=H(:,:,2,q);
            H_2_q=H_2_q(mask_indx);
            H_3_q=H(:,:,3,q);
            H_3_q=H_3_q(mask_indx);
                
            gamma_q= ( A_q(mask_indx).^C(mask_indx) ) ./  ( I_q(mask_indx).^( C(mask_indx)-1)) ;
            
            dI_p_q= I_p(mask_indx) - I_q(mask_indx);   
            
            D_p_q=2./( ( C(mask_indx)-1).*( 1./I_p(mask_indx)-1./I_q(mask_indx)));
 %% R           
            r = in_circle_next(q,2);
            I_r=I(:,:,r);                    
                    
            A_r=A(:,:,r);
            H_1_r=H(:,:,1,r);
            H_1_r=H_1_r(mask_indx);
            H_2_r=H(:,:,2,r);
            H_2_r=H_2_r(mask_indx);
            H_3_r=H(:,:,3,r);
            H_3_r=H_3_r(mask_indx);
                    
            gamma_r= ( A_r(mask_indx).^C(mask_indx) ) ./  ( I_r(mask_indx).^( C(mask_indx)-1)) ;
                    
            dI_p_r= I_p(mask_indx) - I_r(mask_indx);   
                    
            D_p_r=2./( ( C(mask_indx)-1).*( 1./I_p(mask_indx)-1./I_r(mask_indx)));        
       %% p-q-r        
        B_p_q_r=D_p_q.*dI_p_q-D_p_r.*dI_p_r;     
        
        x_f=x(mask_indx)./f;
        y_f=y(mask_indx)./f;
        z_f=(Z(mask_indx)+f)./f;
                    
        %%set b1 (actually b2)      
        b2(mask_indx) = B_i_j_k.*(D_p_q.*gamma_p.*(H_1_p - x_f.*H_3_p)-D_p_q.*gamma_q.*(H_1_q - x_f.*H_3_q) ...
                                -D_p_r.*gamma_p.*(H_1_p - x_f.*H_3_p)+D_p_r.*gamma_r.*(H_1_r - x_f.*H_3_r) ) ...
                        -B_p_q_r.*(D_i_j.*gamma_i.*(H_1_i - x_f.*H_3_i)-D_i_j.*gamma_j.*(H_1_j - x_f.*H_3_j) ...
                                -D_i_k.*gamma_i.*(H_1_i - x_f.*H_3_i)+D_i_k.*gamma_k.*(H_1_k - x_f.*H_3_k) );   
                            
        %%set b2 (actually b1)                     
        b1(mask_indx) = B_i_j_k.*(D_p_q.*gamma_p.*(H_2_p - y_f.*H_3_p)-D_p_q.*gamma_q.*(H_2_q - y_f.*H_3_q) ...
                                -D_p_r.*gamma_p.*(H_2_p - y_f.*H_3_p)+D_p_r.*gamma_r.*(H_2_r - y_f.*H_3_r) ) ...
                        -B_p_q_r.*(D_i_j.*gamma_i.*(H_2_i - y_f.*H_3_i)-D_i_j.*gamma_j.*(H_2_j - y_f.*H_3_j) ...
                                -D_i_k.*gamma_i.*(H_2_i - y_f.*H_3_i)+D_i_k.*gamma_k.*(H_2_k - y_f.*H_3_k) );           
        %set s         
	    s12(mask_indx)=z_f.*B_i_j_k.*(D_p_q.*gamma_p.*H_3_p - D_p_q.*gamma_q.*H_3_q - D_p_r.*gamma_p.*H_3_p + D_p_r.*gamma_r.*H_3_r) ...
                          -z_f.*B_p_q_r.*(D_i_j.*gamma_i.*H_3_i - D_i_j.*gamma_j.*H_3_j - D_i_k.*gamma_i.*H_3_i + D_i_k.*gamma_r.*H_3_k);         
          

%                       
%                       nanmean(b2(:))
%                       nanmean(b1(:))
%                       nanmean(s12(:))
                      
                      
                    b1( abs(I_i-I_j)<shadow_threshold) = 0;
                    b2( abs(I_i-I_j)<shadow_threshold) = 0;
                    s12( abs(I_i-I_j)<shadow_threshold) = 0;   
                    
                    b1( abs(I_i-I_k)<shadow_threshold) = 0;
                    b2( abs(I_i-I_k)<shadow_threshold) = 0;
                    s12( abs(I_i-I_k)<shadow_threshold) = 0;   
                nr= sum(abs(dI_p_q)<shadow_threshold | abs(dI_i_j)<shadow_threshold);
                nrp=100*nr/size(dI_p_q,1);
                fprintf(1,'removed %d pixels (%.2f %%)due to low dif (<%.2f)\n', nr, nrp,shadow_threshold);            
% %             
% %                     b1(I_j<shadow_threshold) = 0;
% %                     b2(I_j<shadow_threshold) = 0;
% %                     s12(I_j<shadow_threshold) = 0;%     
% %                     
% %                     b1(I_q<shadow_threshold) = 0;
% %                     b2(I_q<shadow_threshold) = 0;
% %                     s12(I_q<shadow_threshold) = 0;% 
                    b1(abs(I_p-I_q)<shadow_threshold) = 0;
                    b2(abs(I_p-I_q)<shadow_threshold) = 0;
                    s12(abs(I_p-I_q)<shadow_threshold) = 0;%  
%                     
                    b1(abs(I_p-I_r)<shadow_threshold) = 0;
                    b2(abs(I_p-I_r)<shadow_threshold) = 0;
                    s12(abs(I_p-I_r)<shadow_threshold) = 0;%  
%           
          %% SATURATION Threshold       
% %           sum(sum((I_k>0.45)))
%              
%                    b1(I_i>saturation_thresh) = 0;
%                    b2(I_i>saturation_thresh) = 0;
%                    s12(I_i>saturation_thresh) = 0;
% % %             
%                    b1(I_j>saturation_thresh) = 0;
%                    b2(I_j>saturation_thresh) = 0;
%                    s12(I_j>saturation_thresh) = 0;%    
%                    
%                    b1(I_q>saturation_thresh) = 0;
%                    b2(I_q>saturation_thresh) = 0;
%                    s12(I_q>saturation_thresh) = 0;
% %                    
%                    b1(I_r>saturation_thresh) = 0;
%                    b2(I_r>saturation_thresh) = 0;
%                    s12(I_r>saturation_thresh) = 0;
%                    
%                    ss=(I_i>saturation_thresh | I_j>saturation_thresh | I_q>saturation_thresh | I_r>saturation_thresh ) ;
%                    ns=sum(sum(ss));
%                    nsp=100*ns/size(dI_q_r,1);
%                    fprintf(1,'removed %d pixels (%.2f %%)due to saturation (>%.2f)\n', ns, nsp,saturation_thresh);

            fprintf(1, '%d %d %d === %d %d %d \n', ii,jj,k,p,q,r);
             
            b(:,:,current_index,1) = b1;
            b(:,:,current_index,2) = b2;
            s(:,:,current_index) = s12;
            current_index = current_index+1;
        end                 
	end   
end 

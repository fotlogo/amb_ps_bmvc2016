function [C_refined]=re_estimate_c(I,C,N,H,A,mask,shadow_threshold,ambient)  
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if ambient==0
        C_refined=re_estimate_c_dark(I,C,N,H,A,mask,shadow_threshold);
    else
        C_refined=re_estimate_c_ambient(I,C,N,H,A,mask,shadow_threshold);
    end
end
function [C_refined]=re_estimate_c_dark(I,C,N,H,A,mask,shadow_threshold)    
    disp('C estimation...');
    [nrows,ncols,nb_images] = size(I);
    indices_mask = (mask>0);
% 	t0 = tic; 
%% calculate n.H for each image per pixel and rank images accordingly. This prety much determines photometric parallax
    NDotH=zeros(nrows,ncols,nb_images);   
    
    Cmap=zeros(nrows,ncols);
    
	for ii = 1:nb_images          
        Hi_dot_N=zeros(nrows,ncols);
        for k=1:size(N,3) %i.e. 3           
            Nk=N(:,:,k);
            Hki=H(:,:,k,ii);
            Hi_dot_N(indices_mask)= Hi_dot_N(indices_mask)+ Hki(indices_mask).*Nk(indices_mask);
        end  
        
        Hi_dot_N(I(:,:,ii)<shadow_threshold)=-1;
        
        NDotH(:,:,ii)=Hi_dot_N;
    end
%      m=  nanmean(NDotH,3);
%      ss=  nanstd(NDotH,[],3);
%      
%      figure;
%      imshow(m);
%      figure;
%      imshow(ss);
%      pause    
    %% sort according to photo parallax
     [~,ind] = sort(NDotH,3,'descend');    
    
    %% loop over pixles : BAD PERFORMANCE
    for jj=1: ncols   
        for ii=1: nrows   
            if( mask(ii,jj) ==0)
                continue;
            end  
            if( isnan(N(ii,jj,1)))
                continue;
            end  

            N_H = squeeze(NDotH(ii,jj,ind(ii,jj,:)));  %    squeeze(NDotH(ind(ii,jj,:)))    
            
            last_indx=nb_images; 
            %% Discard light sources where they are in shadow
            while( N_H(last_indx)< 3*shadow_threshold)
                last_indx=last_indx-1;      
                
                if(last_indx<2)
                    break;
                end
            end   
            
            if(last_indx<2)
                continue;
            end            
            num=0;
            denom=0;
            
            for q=1 : floor( (last_indx+1)/4 )
                r=last_indx+1-q;
                
                indx_q=ind(ii,jj,q);
                indx_r=ind(ii,jj,r);  
                
                current_num=log(NDotH(ii,jj,indx_q)) - log(NDotH(ii,jj,indx_r));
                current_denom=log(I(ii,jj,indx_q)) - log(I(ii,jj,indx_r)) - log(A(ii,jj,indx_q)) + log(A(ii,jj,indx_r)) ;  
                %if no point passes constraints then it will not get
                %updated
                if(current_denom>0) %numerator is always positive cause of sorting
                    if (current_num>current_denom) %it has to be true
                        num=num + current_num;
                        denom=denom +current_denom;    
                    end
                end
%                 fprintf(1, 'num: (%f) - (%f) \n', log(NDotH(ii,jj,indx_q)),log(NDotH(ii,jj,indx_r)));
%                 fprintf(1, 'denomnum: (%f) - (%f) -(%f) +(%f) \n', log(I(ii,jj,indx_q)),log(I(ii,jj,indx_r)),log(A(ii,jj,indx_q)),log(A(ii,jj,indx_r)));
%                 fprintf(1,'c est %f \n', num/denom);
%                 pause
            end
            d=num/(num-denom);
            Cmap(ii,jj)=d/(d+1);
        
        end
    end   
    
    %% JUST DONT UPDATE AT UNKNOWNS
    indx=Cmap > 1.0;
    Cmap (indx)= C(indx);    
    indx=Cmap < 0.0;
    Cmap (indx)= C(indx);   
    indx=isnan(Cmap);
    Cmap (indx)= C(indx);      
    
	Cmap=imresize(Cmap,[nrows,ncols]/4); 
    Cmap=imresize(Cmap,[nrows,ncols]); 
    
 
%     figure;
%     hist(Cmap(:),100);
%     pause
   
% 	pause		
%  	C_median = nanmedian(Cmap(indices_mask))
%     C_mean = nanmean(Cmap(indices_mask))

%     c2=zeros(nrows,ncols,2);
%     c2(:,:,1)=C;
%     c2(:,:,2)=Cmap;
% 
%     C_refined=nanmean(c2,3); %C_median*ones(nrows,ncols);
      C_refined=0.15*C+0.85*Cmap;
    
%     figure;
%     imshow(C_refined);
%     pause
    
      C_mean_new = nanmean(C_refined(indices_mask))
%     C_refined= Cmap;
    
	%~ C(find(isnan(C))) = 1;
	%~ C(mask==0) = 1;
	%~ C = medfilt2(real(C),[5 5],'symmetric');			
	%~ C(find(C<2*eps)) = 1;
	%~ clear Cmap num denom
		
% 	t1 = toc(t0);
% 	fprintf('Elapsed time : %.2f\n',t1);		
end
function [C_refined]=re_estimate_c_ambient(I,C,N,H,A,mask,shadow_threshold)    
    disp('C estimation (ambient version)...');
    [nrows,ncols,nb_images] = size(I);
    indices_mask = (mask>0);
% 	t0 = tic; 
%% calculate n.H for each image per pixel and rank images accordingly. This prety much determines photometric parallax
    NDotH=zeros(nrows,ncols,nb_images);   
    
    Cmap=zeros(nrows,ncols);
    
	for yy = 1:nb_images          
        Hi_dot_N=zeros(nrows,ncols);
        for k=1:size(N,3) %i.e. 3           
            Nk=N(:,:,k);
            Hki=H(:,:,k,yy);
            Hi_dot_N(indices_mask)= Hi_dot_N(indices_mask)+ Hki(indices_mask).*Nk(indices_mask);
        end  
        
        Hi_dot_N(I(:,:,yy)<shadow_threshold)=-1;
        
        NDotH(:,:,yy)=Hi_dot_N;
    end
%      m=  nanmean(NDotH,3);
%      ss=  nanstd(NDotH,[],3);
%      
%      figure;
%      imshow(m);
%      figure;
%      imshow(ss);
%      pause    
    %% sort according to photo parallax
     [~,ind] = sort(NDotH,3,'descend');    
    
    %% loop over pixles : BAD PERFORMANCE
    for xx=1: ncols   
        for yy=1: nrows   
            if( mask(yy,xx) ==0)
                continue;
            end  
            if( isnan(N(yy,xx,1)))
                continue;
            end  
            N_H = squeeze(NDotH(yy,xx,ind(yy,xx,:)));  %    squeeze(NDotH(ind(ii,jj,:)))                
            last_indx=nb_images; 
            %% Discard light sources where they are in shadow
            while( N_H(last_indx)< 3*shadow_threshold)
                last_indx=last_indx-1;      
                
                if(last_indx<2)
                    break;
                end
            end   
            
            if(last_indx<2)
                continue;
            end            
           
            nunmber_of_quadruples=floor( (last_indx+1)/4 );
            
            crev=zeros(nunmber_of_quadruples,1);
            
            for ii=1 : nunmber_of_quadruples
%                 ii
                jj= 2*ii;
                q = last_indx+1-2*ii;
                r = last_indx+1-ii;                
                
                indx_i=ind(yy,xx,ii);
                indx_j=ind(yy,xx,jj);
                indx_q=ind(yy,xx,q);  
                indx_r=ind(yy,xx,r);  
                
                dI_ij= I(yy,xx,indx_i)-I(yy,xx,indx_j);
                dI_qr= I(yy,xx,indx_q)-I(yy,xx,indx_r);
                
                %write down the equation as sum( Ki*Mi^x)=0 where x=1/c
                %K1=dI_qr*ai Mi=NDotH_i etc                
                K=zeros(4,1);
                M=zeros(4,1);
                
                K(1)= dI_qr.*A(yy,xx,indx_i);
                K(2)= -dI_qr.*A(yy,xx,indx_j);
                K(3)= -dI_ij.*A(yy,xx,indx_q);
                K(4)= dI_ij.*A(yy,xx,indx_r);
                
                M(1)=NDotH(yy,xx,indx_i);
                M(2)=NDotH(yy,xx,indx_j);
                M(3)=NDotH(yy,xx,indx_q);
                M(4)=NDotH(yy,xx,indx_r);
                
                crev(ii)=newton_solve(K,M,1./C(yy,xx));    
            end
           
            Cmap(yy,xx)=1/nanmean(crev);
        
        end
    end   
    
%     figure;
%     imshow(Cmap);
% %     pause
    
    %% JUST DONT UPDATE AT UNKNOWNS
    indx=Cmap > 1.0;
    Cmap (indx)= C(indx);    
    indx=Cmap < 0.0;
    Cmap (indx)= C(indx);   
    indx=isnan(Cmap);
    Cmap (indx)= C(indx);      
    
	Cmap=imresize(Cmap,[nrows,ncols]/4); 
    Cmap=imresize(Cmap,[nrows,ncols]); 
    
 
%     figure;
%     hist(Cmap(:),100);
%     pause
   
% 	pause		
%  	C_median = nanmedian(Cmap(indices_mask))
%     C_mean = nanmean(Cmap(indices_mask))

%     c2=zeros(nrows,ncols,2);
%     c2(:,:,1)=C;
%     c2(:,:,2)=Cmap;
% 
%     C_refined=nanmean(c2,3); %C_median*ones(nrows,ncols);
       C_refined=0.15*C+0.85*Cmap;
%         C_refined=0.5*C+0.5*Cmap;
%        C_refined=Cmap;
      
       C_mean_new = nanmean(C_refined(indices_mask))
    
%     figure;
%     imshow(C_refined);
%     pause
    
      
%     C_refined= Cmap;
    
	%~ C(find(isnan(C))) = 1;
	%~ C(mask==0) = 1;
	%~ C = medfilt2(real(C),[5 5],'symmetric');			
	%~ C(find(C<2*eps)) = 1;
	%~ clear Cmap num denom
		
% 	t1 = toc(t0);
% 	fprintf('Elapsed time : %.2f\n',t1);		
end
%solve for sum( Ki*Mi^x)=0 
%enforce that x>1
%use a regulariser
function xk=newton_solve(K,M,x0)

% lamda=1e-2;
lamda=0;

xk=x0;
for k=1:10   
    fx=K(1).*(M(1).^xk)+K(2).*(M(2).^xk)+K(3).*(M(3).^xk)+K(4).*(M(4).^xk) ;
    fx_reg=fx.^2+lamda*((xk-1).^2);
    fdotx=K(1).*log(M(1)).*(M(1).^xk)+K(2).*log(M(2)).*(M(2).^xk)+K(3).*log(M(3)).*(M(3).^xk)+K(4).*log(M(4)).*(M(4).^xk);
    fdotx_reg=2*fx.*fdotx+2*(xk-1);
%     xk=xk-fx./fdotx;
    xk=xk-fx_reg./fdotx_reg;
    if(abs(fx)<1e-16)
        break;
    end
    if(xk<1)
        break;
    end
end
    if(xk<1)
        xk=NaN;
    end
%      if(isnan(xk))
%          
%      else
%      fprintf(1, 'x=%f f=%e \n',xk,fx);
%      end

end
%%version with albedo
function [C_refined]=re_estimate_c2(I,C,N,H,A,mask,shadow_threshold)    
    disp('C estimation...');
    [nrows,ncols,nb_images] = size(I);
    indices_mask = (mask>0);
    C_rev=1./C;
% 	t0 = tic; 
%% calculate n.H for each image per pixel and rank images accordingly. This prety much determines photometric parallax
    NDotH=zeros(nrows,ncols,nb_images);     
    
	for ii = 1:nb_images          
        Hi_dot_N=zeros(nrows,ncols);
        for k=1:size(N,3) %i.e. 3           
            Nk=N(:,:,k);
            Hki=H(:,:,k,ii);
            Hi_dot_N(indices_mask)= Hi_dot_N(indices_mask)+ Hki(indices_mask).*Nk(indices_mask);
        end  
        
        Hi_dot_N(I(:,:,ii)<shadow_threshold)=NaN;
        
        NDotH(:,:,ii)=Hi_dot_N;
    end

    rhomap=zeros(nrows,ncols,nb_images);
    rho_i=zeros(nrows,ncols);
    
    for ii = 1:nb_images        
        Ii=I(:,:,ii);
        Ii=Ii(indices_mask);
        
        ai=A(:,:,ii);
        ai=ai(indices_mask);   
        
        NDotHi=NDotH(:,:,ii);
        NDotHi=NDotHi(indices_mask);   
        NDotHi=max(NDotHi,0.001);        
        
        rho_i(indices_mask)= Ii./( ai.* (NDotHi.^C_rev(indices_mask)));
        rhomap(:,:,ii)=rho_i;
    end
    
    rho=nanmean(rhomap,3);
    rho=max(rho,0.01);    
    rho(rho>1)=NaN;
    
%     am=nanmean(A,3);
    
    figure;
%     imshow(am);
     imshow(rho);  
%     pause
%         
    mrho=nanmean(rho(:))
    srho=nanstd(rho(:))
   
    
    Cmap_full=zeros(nrows,ncols,nb_images);
    Cmap_i=nan(nrows,ncols);
    
    for ii = 1:nb_images        
        Ii=I(:,:,ii);
        Ii=Ii(indices_mask);
        Ii=max(Ii,0.001);
        
        ai=A(:,:,ii);
        ai=ai(indices_mask);
        ai=max(ai,0.001);
        
        NDotHi=NDotH(:,:,ii);
        NDotHi=NDotHi(indices_mask);        
        NDotHi=max(NDotHi,0.001);
        
        %rho same everywhere
        rho_i=rho(indices_mask); 
        rho_i=max(rho_i,0.001);       
     
        
        Cmap_i(indices_mask)= log(NDotHi)./( log(Ii)-log(rho_i)-log(ai));
%         find(Ii<0)
%         find(rho_i<0)
%         find(ai<0)
        Cmap_full(:,:,ii)=Cmap_i;        
    end
    
    Cmap=nanmedian(Cmap_full,3);
    

%% JUST DONT UPDATE AT UNKNOWNS
    indx=Cmap > 1.0;
    Cmap (indx)= C(indx);    
    indx=Cmap < 0.0;
    Cmap (indx)= C(indx);   
    indx=isnan(Cmap);
    Cmap (indx)= C(indx);      
    
	Cmap=imresize(Cmap,[nrows,ncols]/8); 
    Cmap=imresize(Cmap,[nrows,ncols]); 
    
 
%     figure;
%     hist(Cmap(:),100);
%     pause
   
% 	pause		
%  	C_median = nanmedian(Cmap(indices_mask))
%     C_mean = nanmean(Cmap(indices_mask))

%     c2=zeros(nrows,ncols,2);
%     c2(:,:,1)=C;
%     c2(:,:,2)=Cmap;
% 
%     C_refined=nanmean(c2,3); %C_median*ones(nrows,ncols);
      C_refined=0.5*C+0.5*Cmap;
    
    figure;
    imshow(C_refined);
%     pause
    
      C_mean_new = nanmean(C_refined(indices_mask))
%     C_refined= Cmap;
    
	%~ C(find(isnan(C))) = 1;
	%~ C(mask==0) = 1;
	%~ C = medfilt2(real(C),[5 5],'symmetric');			
	%~ C(find(C<2*eps)) = 1;
	%~ clear Cmap num denom
		
% 	t1 = toc(t0);
% 	fprintf('Elapsed time : %.2f\n',t1);		
end

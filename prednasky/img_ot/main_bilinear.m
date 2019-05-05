function main_priklad

% nageneruju data
p = 225;
n = 20;
er = 2;
data = 2;
switch data
    case 1
        A = zeros(p,er);
        pol = 4;
        x = 7;
        y = 6;
        kolo1 = zeros(sqrt(p),sqrt(p));
        for i = 1:sqrt(p)
            for j = 1:sqrt(p)
                if sqrt( (i-x)^2 + (j-y)^2 ) <= pol
                    kolo1(i,j) = 1;
                end
            end
        end
        x = 10;
        y = 10;
        kolo2 = zeros(sqrt(p),sqrt(p));
        for i = 1:sqrt(p)
            for j = 1:sqrt(p)
                if sqrt( (i-x)^2 + (j-y)^2 ) <= pol
                    kolo2(i,j) = 1;
                end
            end
        end
    %     figure(1); imagesc(kolo1+kolo2);

        A(:,1) = kolo1(:);
        A(:,2) = kolo2(:);

        X = zeros(n,er);
        X(:,1) = [(0:(0.1):0.9)'; ones(5,1); zeros(5,1)];
        X(:,2) = [ones(5,1); ones(5,1); (1:(-0.1):0.1)'];

    %     tisk(A,X,1)
        D = A*X' + 0.1*randn(size(A*X'));
        D(D<0) = 0;

    %     tisk_D(D,2);

    %     save ./cvicna_data.mat D A X
    case 2
        load ./cvicna_data.mat
    case 3
        D=imgread('/home/otichy/ot/scintigraphy/data/rnu11.img',64*64);
    case 4
        D=imgread('/home/otichy/ot/scintigraphy/data/studie_rakousko_crop_correct/IM3crop.img',64*64);
    case 5
        D=imgread('/home/otichy/ot/scintigraphy/data/studie_rakousko_crop_correct/IM2.img',64*64); 
    case 6
        D=imgread('/home/otichy/ot/scintigraphy/data/studie_rakousko_crop_correct/IM1.img',64*64);
    case 7
        Dorig=multibandread('/home/otichy/ot/scintigraphy/SNPA/Terrain/TERRAIN',[500 307 210],'uint16',0,'bil','ieee-be');
        
        sqp = 150;
        en = 1:210;
        
        D = zeros(sqp^2,size(en,2));
        
        for j = en
            Do = Dorig(1:sqp,1:sqp,j);
            D(:,(j-en(1)+1)) = Do(:);
        end
        save ../8_cviko/cviko_hyperspektral.mat D
        error('jo')
end


%% tisk fantomu - dat i zdroju
if 0
    tisk(A,X,1)
    tisk_D(D,2);
    
    figure; 
    imagesc(D); 
%     set(gca,'YTick',[]); set(gca,'XTick',[]); 
    colormap(1-gray)

    error('bagr')
end

[p,n] = size(D);

%% realna data, otoceni, tisk D
if 0
    Do = zeros(sqrt(p));
    for j=1:n
        Do(:) = D(:,j);
        Do = Do';
        D(:,j) = Do(:);
    end
    size(D)
%     save ./data_rnu11.mat D
    
end

if 0
    figure(24);
    radku = 4;
    sloupcu = 8;
    colormap(1-gray);
    Do = zeros(sqrt(p));
    for i=1:(radku*sloupcu)
        j = (i-1)*3+1;%i*2 - 0;
%             if j > 40 % PET, d255
%                 j = 40;
%             end
%         j=i;
        display(j);
        Do(:) = D(:,j);
        subplot(radku,sloupcu,i); imagesc(Do); set(gca,'YTick',[]); set(gca,'XTick',[]); axis('square');
    end
    
    error('baf')
    
    if 0 % scaling of the data
        sOn = 1./sum(D,1)';
        sOp = 1./sum(D,2);
        scD = D.*(sqrt(sOp)*sqrt(sOn)');
    end
end


%% PCA
if 0 % nejak nefunguje, asi neco blbe...
    r = 3;
    C = D'*D/(n-1);
    [V,L] = ldl(C);
    pc = D*V;
%     keyboard
    figure; 
    rows = 1;
    cols = r;
    for i = 1:r
        subplot(rows,cols,i)
            imagesc(reshape(pc(:,i),sqrt(p),sqrt(p)));
            set(gca,'YTick',[]); set(gca,'XTick',[]);
            colormap(1-gray)
            axis('square');
    end
end


%% SVD
if 0
    r = 2;
    [U,S,V] = svd(D);
    pc = U*S(:,1:r);
%     keyboard

    tisk(pc(:,1:r),V(:,1:r),3);
    
    % zpetna rekonstrukce
    tisk_D(pc(:,1:r)*V(:,1:r)',4)


%     figure; 
%     rows = 2;
%     cols = r;
%     for i = 1:r
%         subplot(rows,cols,i)
%             imagesc(reshape(pc(:,i),sqrt(p),sqrt(p)));
%             set(gca,'YTick',[]); set(gca,'XTick',[]);
%             colormap(1-gray)
%             axis('square');
%         
%     end
end


%% NMF
if 0
    r = 3;
    hA = rand(p,r);
    hX = rand(n,r);
    for i = 1:100
        hX = (hX'.*( (hA'*D)./(hA'*hA*hX') ) )';
        hA = hA.*( (D*hX)./(hA*hX'*hX) );
    end
    tisk(hA,hX,4)
end

%% VB - jen I_r vsude
if 1
    r = 3;
    
    theta0 = 1e-10;
    rho0 = 1e-10;
    hat_omega = 1;
    hat_XtX = eye(r);
    hat_X = ones(n,r);
%     hat_X = rand(n,r);
    
    hat_XX = zeros(n,r);
    hat_AA = zeros(p,r);
    hat_A = zeros(p,r);
    
    for i = 1:100
%         i
        % A
        sA = ( hat_omega * hat_XtX + eye(r))^(-1);
        mA = hat_omega * D * hat_X * sA; 
        
        hat_A = mA;
        hat_AtA = mA'*mA + p*sA;
        
        if 1
            hat_A = mA;
%             hat_A(hat_A<0) = 0.00001;
            hat_AtA = mA'*mA + p*sA;
        else
            blokPhiA = ones(p,1)*diag(sA)';
            diagkronPhiA = reshape(blokPhiA,p*r,1);
            [hat_A(:) hat_AA(:)] = momtrun_low(mA(:),sqrt(diagkronPhiA),0);
            varAA = hat_AA - hat_A.^2;
            hat_AtA = hat_A'*hat_A + diag(sum(varAA));    
        end
        
        % X
        sX = ( hat_omega * hat_AtA + eye(r))^(-1);
        mX = hat_omega * D' * hat_A * sX;
        
        if 1
            hat_X = mX;
%             hat_X(hat_X<0) = 0.00001;
            hat_XtX = mX'*mX + n*sX;
        else
            [hat_X(:) hat_XX(:)]= momtrun_low(mX(:),sqrt(diag( kron(sX,eye(n)) )),0);
            varXX = hat_XX - hat_X.^2;
            hat_XtX = hat_X'*hat_X + diag(sum(varXX)); 
        end
        
        % omega
        theta = theta0 + n*p/2;
        rho = rho0 + (1/2)*trace(D'*D - D'*hat_A*hat_X' - hat_X*hat_A'*D) + (1/2)*trace(hat_AtA*hat_XtX);
%         rho = rho0 + (1/2)*(trace(D*D' - hat_A*hat_X'*D' - D*hat_X*hat_A') + (1/2)*trace(hat_AtA*hat_XtX));
        hat_omega = theta/rho;
    end
    tisk(hat_A,hat_X,5);
%     keyboard
end


%% VB - ARD na A
if 0
    display('poustim VB s ARD')
    r = 3;
    
    alpha0 = 1e-10;
    beta0 = 1e-10;
    theta0 = 1e-10;
    rho0 = 1e-10;
    hat_omega = 1;
    hat_XtX = eye(r);
    hat_X = ones(n,r);
    hat_v = ones(r,1);
    
    hat_XX = zeros(n,r);
    hat_AA = zeros(p,r);
    hat_A = zeros(p,r);
    
    info_v = zeros(100,r);
    
    for i = 1:100
        
        % A
        sA = ( hat_omega * hat_XtX + diag(hat_v))^(-1);
        mA = hat_omega * D * hat_X * sA; 
        
        if 0
            hat_A = mA;
%             hat_A(hat_A<0) = 0.0000;
            hat_AtA = mA'*mA + p*sA;
        else
            blokPhiA = ones(p,1)*diag(sA)';
            diagkronPhiA = reshape(blokPhiA,p*r,1);
            [hat_A(:) hat_AA(:)] = momtrun_low(mA(:),sqrt(diagkronPhiA),0);
            varAA = hat_AA - hat_A.^2;
            hat_AtA = hat_A'*hat_A + diag(sum(varAA));    
        end
        
        % v na A
        alpha = alpha0 + (p/2)*ones(r,1);
        beta = beta0 + (1/2)*diag(hat_AtA);
        hat_v = (alpha./beta);
        info_v(i,:) = hat_v; 
    
        % X
        sX = ( hat_omega * hat_AtA + eye(r))^(-1);
        mX = hat_omega * D' * hat_A * sX;
        
        if 0
            hat_X = mX;
%             hat_X(hat_X<0) = 0.0000;
            hat_XtX = mX'*mX + n*sX;
        else
            [hat_X(:) hat_XX(:)]= momtrun_low(mX(:),sqrt(diag( kron(sX,eye(n)) )),0);
            varXX = hat_XX - hat_X.^2;
            hat_XtX = hat_X'*hat_X + diag(sum(varXX)); 
        end
        
        % omega
        theta = theta0 + n*p/2;
        rho = rho0 + (1/2)*trace(D'*D - D'*hat_A*hat_X' - hat_X*hat_A'*D) + (1/2)*trace(hat_AtA*hat_XtX);
        hat_omega = theta/rho;
    end
    tisk(hat_A,hat_X,6);

    if 0
        fig = figure(100);
        set(fig, 'Position', [1800, 200, 400, 200]);
        semilogy(info_v);
    end
    
    if 0
        A = hat_A;
        X = hat_X;
        [p,r] = size(A);
        [n,r] = size(X);
        sqp = sqrt(p);

        fig = figure(6);
        set(fig, 'Position', [1800, 200, 200*r, 400]);
        rows = 2;
        cols = r;

        for i = 1:r
            subplot(rows,cols,i)
                imagesc(reshape(A(:,i),sqp,sqp))
                set(gca,'YTick',[]); set(gca,'XTick',[]);
                colormap(1-gray)
                axis('square');
                title(['Zdroj ' num2str(i)])
    %             axis off;
            subplot(rows,cols,r+i)
                plot(X(:,i),'LineWidth',2)
                xlim([1 n])
                ylim([min(0,min(min(X(:,i)))) 1.1*max(max(X(:,i)))])
        end
    end
end



end


function tisk_ard(A,X,fig_no)
    [p,r] = size(A);
    [n,r] = size(X);
    sqp = sqrt(p);
    
    fig = figure(fig_no);
    set(fig, 'Position', [1800, 200, 200*r, 400]);
    rows = 2;
    cols = r;
    
    for i = 1:r
        subplot(rows,cols,i)
            imagesc(reshape(A(:,i),sqp,sqp))
            set(gca,'YTick',[]); set(gca,'XTick',[]);
            colormap(1-gray)
            axis('square');
            title(['Zdroj ' num2str(i)])
%             axis off;
        subplot(rows,cols,r+i)
            if 1
                plot(max(A(:,i))*X(:,i),'LineWidth',2)
                xlim([1 n])
                ylim([min(0,min(min(max(A(:,i))*X(:,i)))) 1.1*max(max(max(A(:,i))*X(:,i)))])
            else
                plot(X(:,i),'LineWidth',2)
                xlim([1 n])
                ylim([min(0,min(min(X(:,i)))) 1.1*max(max(X(:,i)))])
            end
    end
end


function tisk(A,X,fig_no)
    [p,r] = size(A);
    [n,r] = size(X);
    sqp = sqrt(p);
    
    fig = figure(fig_no);
    set(fig, 'Position', [1800, 200, 200*r, 400]);
    rows = 2;
    cols = r;
    
    for i = 1:r
        subplot(rows,cols,i)
            imagesc(reshape(A(:,i),sqp,sqp))
            set(gca,'YTick',[]); set(gca,'XTick',[]);
            colormap(1-gray)
            axis('square');
            title(['Zdroj ' num2str(i)])
%             axis off;
        subplot(rows,cols,r+i)
            plot(X(:,i),'LineWidth',2)
            xlim([1 n])
            ylim([min(0,min(min(X(:,i)))) 1.1*max(max(X(:,i)))])
    end
end

function tisk_D(D,fig_no)
    [p,n] = size(D);
    mind = min(min(D));
    maxd = max(max(D));
    sqp = sqrt(p);
    fig = figure(fig_no);
    set(fig, 'Position', [1800, 200, 800, 500]);
    rows = 4;
    cols = 5;
    for i = 1:20
        subplot(rows,cols,i)
            imagesc(reshape(D(:,i),sqp,sqp),[mind maxd])
            set(gca,'YTick',[]); set(gca,'XTick',[]);
            colormap(1-gray)
            axis('square');
            title(i)
%             axis off;
    end
end





function DH_simulator(input_filepath, processDataButtonHandle)
    originalButtonText = processDataButtonHandle.Text;
    cleanup = onCleanup(@()set(processDataButtonHandle,'Text',originalButtonText,'Icon',''));
    processDataButtonHandle.Text = 'Processing...';
    processDataButtonHandle.IconAlignment = 'bottom';
    wbar = permute(repmat(processDataButtonHandle.BackgroundColor,15,1,200),[1,3,2]);
    wbar([1,end],:,:) = 0;
    wbar(:,[1,end],:) = 0;
    processDataButtonHandle.Icon = wbar;


    % シミュレーション変数
    J              = 1i;                            % 虚数単位

    % 実験変数
    lambda     = 532.0e-9;          % レーザーの波長[m]
    pix_imager = 3.45e-6;           % カメラピクセルサイズ（矩形1辺）[m]
    N_imager   = 2048;              % 計算格子ピクセル数（一辺）[pix]

    % マスク切り出し変数
    N_mask      = 256;      % 切り出しマスクサイズ[pix]
    use_left    = false;     % 切り出しマスクの位置，左:true，右:false

    % 伝播計算区間
    d_start = -0.08; % 開始地点[m]，（負の値は逆伝播を表す）
    d_step  = -0.02; % 刻み幅[m]
    d_end   = -0.18; % 終了地点[m]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get hologram
    file_name  = strcat(input_filepath);
    DH_letter  = imread(file_name);
    DH_letter  = im2double(DH_letter);
    hologram   = rgb2gray(DH_letter);

    % zero padding(1440,1080) -> (N_imager,N_imager)
    hologram      = padarray(hologram,[484 304],0,'both');
    figure('Visible', 'off'), imagesc(hologram);
    saveas(gcf,'Figs/hologram.png');

    % ホログラムのフーリエ変換
    FT_hologram = fftshift(fft2(hologram));

    figure('Visible', 'off'), imagesc(abs(FT_hologram).^2), colorbar, caxis([0 1e8]);
    saveas(gcf,'Figs/FT_hologram.png');


    % 切り出しマスクの中心座標の探索
    if use_left
        % 左半分
        FT_hologram_masked_left = FT_hologram(:,1:(N_imager/2));
        FT_hologram_masked_left((N_imager/2-N_mask/2):(N_imager/2+N_mask/2),(N_imager/2-N_mask/2):N_imager/2)=0;

        [~,max_idx]=max(abs(FT_hologram_masked_left(:)));             % ベクトル化して最大値とインデックスを探索
        [Cy, Cx] = ind2sub(size(FT_hologram_masked_left),max_idx);
    else
        % 右半分
        FT_hologram_masked_right = FT_hologram(:,(N_imager/2)+1:N_imager);
        FT_hologram_masked_right((N_imager/2-N_mask/2):(N_imager/2+N_mask/2),1:N_mask/2)=0;
        [~,max_idx]=max(abs(FT_hologram_masked_right(:)));
        [Cy, Cx] = ind2sub(size(FT_hologram_masked_right),max_idx);
        Cx = Cx + N_imager/2;
    end


    % 成分を切り出し
    FT_hologram_masked_real = imcrop(real(FT_hologram), [Cx-N_mask/2+1 Cy-N_mask/2+1 N_mask-1 N_mask-1]);
    FT_hologram_masked_imag = imcrop(imag(FT_hologram), [Cx-N_mask/2+1 Cy-N_mask/2+1 N_mask-1 N_mask-1]);
    FT_hologram_masked_tmp  = FT_hologram_masked_real + J*FT_hologram_masked_imag;

    % ゼロパディング
    FT_hologram_masked  = padarray(FT_hologram_masked_tmp,[(N_imager-N_mask)/2 (N_imager-N_mask)/2],0,'both');

    figure('Visible', 'off'), imagesc(abs(FT_hologram_masked).*abs(FT_hologram_masked)), colorbar, caxis([0 1e8]);
    saveas(gcf,'Figs/FT_hologram_masked.png');

    clear FT_hologram_masked_real FT_hologram_masked_imag FT_hologram_masked_tmp

    % 逆フーリエ変換
    camp_read_imager = ifft2(FT_hologram_masked);

    figure('Visible', 'off'), imagesc(abs(camp_read_imager).^2);
    saveas(gcf,'Figs/Reconst_intensity_0cm.png');

    figure('Visible', 'off'), imagesc(angle(camp_read_imager)), colorbar, caxis([0 6.28]), colormap(hsv);
    saveas(gcf,'Figs/Reconst_phase_0cm.png');

    % 【Read-4】伝搬計算
    fresnel_read = zeros(N_imager,N_imager);
    step = 0;
    total_step = (d_end - d_start)/d_step;
    for d_read = d_start:d_step:d_end
        for j=1:N_imager
            for i=1:N_imager
                fresnel_read(j,i) = exp((J*pi/(lambda*d_read))*(((i-N_imager/2)*pix_imager)^2+((j-N_imager/2)*pix_imager)^2));
                % x = (i-N_imager/2)*pix_imager, y = (j-N_imager/2)*pix_imager
            end
        end
        camp_read	= ((pix_imager*pix_imager)/(J*lambda*d_read))*fftshift(ifft2(FT_hologram_masked.*fftshift(fft2(fresnel_read)))) ;
        figure('Visible', 'off'), imagesc(angle(fresnel_read));
        saveas(gcf,'Figs/fresnel_read.png');
        % ★★★ FT_hologram_maskedと二次位相分布のFTを掛け算したのち逆変換
        clear fresnel_read;

        % ---Figure---
        figure('Visible', 'off'), imagesc(abs(camp_read).^2);
        FN_fig = strcat('Figs/Reconst_intensity_',num2str(d_start*100),'cm.png');
        saveas(gcf,FN_fig);

        currentProg = min(round((size(wbar,2)-2)*(step/total_step)),size(wbar,2)-2);
        RGB = processDataButtonHandle.Icon;
        RGB(2:end-1, 2:currentProg+1, 1) = 0.25391; % (royalblue)
        RGB(2:end-1, 2:currentProg+1, 2) = 0.41016;
        RGB(2:end-1, 2:currentProg+1, 3) = 0.87891; 
        processDataButtonHandle.Icon = RGB;
    %     figure('Visible', 'off'), imagesc(angle(camp_read)), colorbar, caxis([0 6.28]), colormap(hsv);
    %     FN_fig = strcat('Figs/Reconst_phase-',num2str(file_d),'cm.png');
    %     saveas(gcf,FN_fig);
        d_start = d_start + d_step;
        step = step + 1;
    end
end
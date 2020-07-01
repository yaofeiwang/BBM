% -------------------------------------------------------------------------
% BBM | Feb 2020 | version 1.0 |
% -------------------------------------------------------------------------
% Contact: yaofei@mail.ustc.edu.cn, Yaofei Wang
% -------------------------------------------------------------------------
% Non-Additive Cost Functions for JPEG Steganography Based on Block Boundary Maintenance
% -------------------------------------------------------------------------

function BBCupBBMexample()
    addpath(genpath(pwd));
    poolnum = str2double(getenv('SLURM_CPUS_PER_TASK'));
    costfun = 'UERD';
    Scale_fun = 'exp';
    payload = 0.4;
    embed_order = 6;
    quality_factor = 75

    DP1 = cell(1, 3);

    for i = 1:numel(DP1)
        DP1{i} = zeros(8, 8);
    end

    DP1{1}(1:2, 1:2) = 1;
    DP1{2}(1:2, 3:4) = 1;
    DP1{2}(3:4, 1:2) = 1;
    DP1{3} = 1 - DP1{1} - DP1{2};

    cover_dir = 'cover';
    stego_dir = ['stego', costfun, 'BBCupBBM_', Scale_fun, num2str(embed_order), '_payload', num2str(payload)];
    if ~exist(stego_dir, 'dir'); mkdir(stego_dir); end

    imgs = dir(cover_dir);
    len = length(imgs)
    default_gray_jpeg_obj = jpeg_read([cover_dir, '/', imgs(3).name]);
    q_tab = default_gray_jpeg_obj.quant_tables{1}
    [W, L] = size(default_gray_jpeg_obj.coef_arrays{1});

    [D4, T_hor, T_ver, R, F_hor2, F_hor4, F_hor6, F_ver2, F_ver4, F_ver6] = TRFDparameter(q_tab, W / 2, L / 2);

    parpool(poolnum);

    parfor i = 3:len
        img_name = imgs(i).name;
        cover_img = [cover_dir, '/', img_name];
        stego_img = [stego_dir, '/', img_name];
        BBCup_BBM(payload, cover_img, stego_img, DP1, D4, T_hor, T_ver, R, F_hor2, F_hor4, F_hor6, F_ver2, F_ver4, F_ver6, embed_order, costfun);
    end

    poolobj = gcp('nocreate');
    delete(poolobj);

end

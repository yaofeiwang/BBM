% -------------------------------------------------------------------------
% UERD | Feb 2017 | version 1.0 | 
% -------------------------------------------------------------------------
% Contact: chenkj@mail.ustc.edu.cn, Kejiang Chen
% -------------------------------------------------------------------------
% Using Statistical Image Model for JPEG Steganography: Uniform Embedding Revisited
% -------------------------------------------------------------------------

function [colorrhoP1, colorrhoM1] = UERD_Cost(cover_path, update_choice)
    colorrhoP1 = [];
    colorrhoM1 = [];
    wetConst = 10^13;
    img = jpeg_read(cover_path);
    [~, channel_num] = size(img.coef_arrays);

    [covermat_path, cover_name, ~] = fileparts(cover_path);

    if channel_num == 3
        covermat_path = [covermat_path, 'UERDColorrho', '/', cover_name, '.mat'];
    elseif channel_num == 1
        covermat_path = [covermat_path, 'UERDGrayrho', '/', cover_name, '.mat'];
    else
        covermat_path = [covermat_path, 'UERD', num2str(channel_num), 'rho', '/', cover_name, '.mat'];
    end

    if ~exist(covermat_path, 'file') | update_choice

        for channel = 1:channel_num
            dct_coef = img.coef_arrays{channel};
            [X, Y] = size(dct_coef);

            dct_coef2 = dct_coef;
            dct_coef2(1:8:end, 1:8:end) = 0; % remove DC coefs;
            q_tab = img.quant_tables{round((channel + 1) / 2)};
            q_tab(1, 1) = 0.5 * (q_tab(2, 1) + q_tab(1, 2));
            q_matrix = repmat(q_tab, [X / 8 Y / 8]);

            %energy of each block
            dct_coef2 = im2col(q_matrix .* dct_coef2, [8 8], 'distinct');
            J2 = sum(abs(dct_coef2));
            J = ones(64, 1) * J2;
            J = col2im(J, [8 8], [X Y], 'distinct');

            % decide = q_matrix./J; % version 1

            pad_size = 8;
            im2 = padarray(J, [pad_size pad_size], 'symmetric'); % energies of eight - neighbor blocks
            size2 = 2 * pad_size;
            im_l8 = im2(1 + pad_size:end - pad_size, 1:end - size2);
            im_r8 = im2(1 + pad_size:end - pad_size, 1 + size2:end);
            im_u8 = im2(1:end - size2, 1 + pad_size:end - pad_size);
            im_d8 = im2(1 + size2:end, 1 + pad_size:end - pad_size);
            im_l88 = im2(1:end - size2, 1:end - size2);
            im_r88 = im2(1 + size2:end, 1 + size2:end);
            im_u88 = im2(1:end - size2, 1 + size2:end);
            im_d88 = im2(1 + size2:end, 1:end - size2);
            JJ = (J + 0.25 * (im_l8 + im_r8 + im_u8 + im_d8) + 0.25 * (im_l88 + im_r88 + im_u88 + im_d88));
            decide_ini = q_matrix ./ JJ; % version 2
            decide_ini = decide_ini / (min(min(decide_ini)));
            decide_ini(isnan(decide_ini)) = wetConst;
            decide_ini(decide_ini > wetConst) = wetConst;

            decide_ini_Color{channel} = decide_ini;

        end

        % if update_choice == 1
        [save_path, ~, ~] = fileparts(covermat_path);
        if ~exist(save_path, 'dir'); mkdir(save_path); end
        save(covermat_path, 'decide_ini_Color');
        % end

    else
        load(covermat_path);

    end

    for channel = 1:channel_num
        dct_coef = img.coef_arrays{channel};
        rhoP1 = decide_ini_Color{channel};
        rhoM1 = decide_ini_Color{channel};
        rhoP1(dct_coef > 1023) = wetConst;
        rhoM1(dct_coef <- 1023) = wetConst;

        if channel_num == 1
            colorrhoP1 = rhoP1;
            colorrhoM1 = rhoM1;
        else
            colorrhoP1{channel} = rhoP1;
            colorrhoM1{channel} = rhoM1;
        end

    end

end

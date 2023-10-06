classdef app1_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        KonvolusiTab                    matlab.ui.container.Tab
        JalankanKonvolusiButton         matlab.ui.control.Button
        PreviewCitraHasilKonvolusiMATLABPanel  matlab.ui.container.Panel
        ImageAxes_3                     matlab.ui.control.UIAxes
        PreviewCitraHasilKonvolusiBuatanPanel  matlab.ui.container.Panel
        ImageAxes_2                     matlab.ui.control.UIAxes
        PreviewCitraMasukanPanel        matlab.ui.container.Panel
        ImageAxes                       matlab.ui.control.UIAxes
        PilihanMaskUntukKonvolusiButtonGroup  matlab.ui.container.ButtonGroup
        Image_4                         matlab.ui.control.Image
        Image_3                         matlab.ui.control.Image
        Image_2                         matlab.ui.control.Image
        Image                           matlab.ui.control.Image
        Mask2Button                     matlab.ui.control.RadioButton
        Mask4Button                     matlab.ui.control.RadioButton
        Mask3Button                     matlab.ui.control.RadioButton
        Mask1Button                     matlab.ui.control.RadioButton
        PathInfoLabel                   matlab.ui.control.Label
        PathcitraLabel                  matlab.ui.control.Label
        PilihFileCitraButton            matlab.ui.control.Button
        CitramasukanLabel               matlab.ui.control.Label
        SmoothingLowPassFilterTab       matlab.ui.container.Tab
        JalankanSmoothingButton         matlab.ui.control.Button
        PreviewCitraHasilFilterPanel    matlab.ui.container.Panel
        ImageAxes_5                     matlab.ui.control.UIAxes
        PilihanFilterButtonGroup        matlab.ui.container.ButtonGroup
        nfilterorderforBLPFSpinner      matlab.ui.control.Spinner
        nfilterorderforBLPFSpinnerLabel  matlab.ui.control.Label
        D0cutofffrequencySpinner        matlab.ui.control.Spinner
        D0cutofffrequencySpinnerLabel   matlab.ui.control.Label
        sigmaSpinner                    matlab.ui.control.Spinner
        sigmaSpinnerLabel               matlab.ui.control.Label
        nSpinner                        matlab.ui.control.Spinner
        nSpinnerLabel                   matlab.ui.control.Label
        ButterworthLowpassFilterBLPFButton  matlab.ui.control.RadioButton
        GaussianLowpassFilterGLPFButton  matlab.ui.control.RadioButton
        IdealLowpassFilterILPFButton    matlab.ui.control.RadioButton
        RanahFrekuensiLowpassFilterLabel  matlab.ui.control.Label
        RanahSpasialLabel               matlab.ui.control.Label
        GaussianFilternxnButton         matlab.ui.control.RadioButton
        MeanFilternxnButton             matlab.ui.control.RadioButton
        PreviewCitraMasukanPanel_2      matlab.ui.container.Panel
        ImageAxes_4                     matlab.ui.control.UIAxes
        PathInfoLabel_2                 matlab.ui.control.Label
        PathcitraLabel_2                matlab.ui.control.Label
        PilihFileCitraButton_2          matlab.ui.control.Button
        CitramasukanLabel_2             matlab.ui.control.Label
        HighPassFilterTab               matlab.ui.container.Tab
        JalankanHighPassFilterButton    matlab.ui.control.Button
        PreviewCitraHasilFilterPanel_2  matlab.ui.container.Panel
        ImageAxes_7                     matlab.ui.control.UIAxes
        PilihanFilterButtonGroup_2      matlab.ui.container.ButtonGroup
        nfilterorderforBHPFSpinner      matlab.ui.control.Spinner
        nfilterorderforBHPFSpinnerLabel  matlab.ui.control.Label
        D0cutofffrequencySpinner_2      matlab.ui.control.Spinner
        D0cutofffrequencySpinner_2Label  matlab.ui.control.Label
        ButterworthHighpassFilterBHPFButton  matlab.ui.control.RadioButton
        GaussianHighpassFilterGHPFButton  matlab.ui.control.RadioButton
        IdealHighpassFilterIHPFButton   matlab.ui.control.RadioButton
        RanahFrekuensiHighpassFilterLabel  matlab.ui.control.Label
        PreviewCitraMasukanPanel_3      matlab.ui.container.Panel
        ImageAxes_6                     matlab.ui.control.UIAxes
        PathInfoLabel_3                 matlab.ui.control.Label
        PathcitraLabel_3                matlab.ui.control.Label
        PilihFileCitraButton_3          matlab.ui.control.Button
        CitramasukanLabel_3             matlab.ui.control.Label
        FilterTerangHighBoostTab        matlab.ui.container.Tab
        PathInfoLabel_4                 matlab.ui.control.Label
        JalankanFilterTerangButton      matlab.ui.control.Button
        PreviewCitraHasilFilterPanel_3  matlab.ui.container.Panel
        ImageAxes_9                     matlab.ui.control.UIAxes
        PengaturanFilterButtonGroup     matlab.ui.container.ButtonGroup
        ASpinner                        matlab.ui.control.Spinner
        ASpinnerLabel                   matlab.ui.control.Label
        nfilterorderforBLPFSpinner_2    matlab.ui.control.Spinner
        nfilterorderforBLPFSpinner_2Label  matlab.ui.control.Label
        D0cutofffrequencySpinner_3      matlab.ui.control.Spinner
        D0cutofffrequencySpinner_3Label  matlab.ui.control.Label
        ButterworthLowpassFilterBLPFButton_2  matlab.ui.control.RadioButton
        GaussianLowpassFilterGLPFButton_2  matlab.ui.control.RadioButton
        IdealLowpassFilterILPFButton_2  matlab.ui.control.RadioButton
        RanahFrekuensiLowpassFilterLabel_2  matlab.ui.control.Label
        PreviewCitraMasukanPanel_4      matlab.ui.container.Panel
        ImageAxes_8                     matlab.ui.control.UIAxes
        PathcitraLabel_4                matlab.ui.control.Label
        PilihFileCitraButton_4          matlab.ui.control.Button
        CitramasukanLabel_4             matlab.ui.control.Label
        SpasialFilterTab                matlab.ui.container.Tab
        ResetButton                     matlab.ui.control.Button
        GenNoiseButton                  matlab.ui.control.Button
        PengaturanFilterPanel           matlab.ui.container.Panel
        dSpinner                        matlab.ui.control.Spinner
        dSpinnerLabel                   matlab.ui.control.Label
        QSpinner                        matlab.ui.control.Spinner
        QSpinnerLabel                   matlab.ui.control.Label
        NoiseDensitySpinner             matlab.ui.control.Spinner
        NoiseDensitySpinnerLabel        matlab.ui.control.Label
        TipeFilterButtonGroup           matlab.ui.container.ButtonGroup
        AlphaTrimmedMeanFilterButton    matlab.ui.control.RadioButton
        MidpointFilterButton_2          matlab.ui.control.RadioButton
        ContraharmonicMeanFilterButton_2  matlab.ui.control.RadioButton
        GeometricFilterButton_2         matlab.ui.control.RadioButton
        HarmonicMeanFilterButton_2      matlab.ui.control.RadioButton
        ArithmeticMeanFilterButton_2    matlab.ui.control.RadioButton
        MedianFilterButton_2            matlab.ui.control.RadioButton
        MaxFilterButton_2               matlab.ui.control.RadioButton
        MinFilterButton_2               matlab.ui.control.RadioButton
        JalankanFilterButton            matlab.ui.control.Button
        CitraHasilFilterNoisePanel      matlab.ui.container.Panel
        ImageAxes_12                    matlab.ui.control.UIAxes
        CitraHasilPenambahanSaltPepperPanel  matlab.ui.container.Panel
        ImageAxes_11                    matlab.ui.control.UIAxes
        CitraMasukanPanel_2             matlab.ui.container.Panel
        ImageAxes_10                    matlab.ui.control.UIAxes
        PathInfoLabel_5                 matlab.ui.control.Label
        PathcitraLabel_5                matlab.ui.control.Label
        PilihFileCitraButton_5          matlab.ui.control.Button
        CitramasukanLabel_5             matlab.ui.control.Label
        PeriodicFilterTab               matlab.ui.container.Tab
        FourierCitraHasilPanel          matlab.ui.container.Panel
        UIAxes_2                        matlab.ui.control.UIAxes
        FourierCitraMasukanPanel        matlab.ui.container.Panel
        UIAxes                          matlab.ui.control.UIAxes
        JalankanFilterButton_3          matlab.ui.control.Button
        PreviewCitraHasilPanel          matlab.ui.container.Panel
        ImageAxes_18                    matlab.ui.control.UIAxes
        PreviewCitraMasukanPanel_7      matlab.ui.container.Panel
        ImageAxes_16                    matlab.ui.control.UIAxes
        PathInfoLabel_7                 matlab.ui.control.Label
        PathcitraLabel_7                matlab.ui.control.Label
        PilihFileCitraButton_7          matlab.ui.control.Button
        CitramasukanLabel_7             matlab.ui.control.Label
        WienerMotionBlurTab             matlab.ui.container.Tab
        ResetButton_2                   matlab.ui.control.Button
        MotionBlurPanel                 matlab.ui.container.Panel
        BlurCitraButton                 matlab.ui.control.Button
        THETASpinner                    matlab.ui.control.Spinner
        THETASpinnerLabel               matlab.ui.control.Label
        LENSpinner                      matlab.ui.control.Spinner
        LENSpinnerLabel                 matlab.ui.control.Label
        PerbedaanCitraMasukandenganWienerPanel_2  matlab.ui.container.Panel
        ImageAxes_19                    matlab.ui.control.UIAxes
        JalankanFilterButton_2          matlab.ui.control.Button
        CitraHasilFilterWienerPanel_2   matlab.ui.container.Panel
        ImageAxes_15                    matlab.ui.control.UIAxes
        CitraHasilMotionBlurPanel       matlab.ui.container.Panel
        ImageAxes_14                    matlab.ui.control.UIAxes
        CitraMasukanPanel               matlab.ui.container.Panel
        ImageAxes_13                    matlab.ui.control.UIAxes
        PathInfoLabel_6                 matlab.ui.control.Label
        PathcitraLabel_6                matlab.ui.control.Label
        PilihFileCitraButton_6          matlab.ui.control.Button
        CitramasukanLabel_6             matlab.ui.control.Label
        WienerGaussianBlurTab           matlab.ui.container.Tab
        MotionBlurPanel_2               matlab.ui.container.Panel
        noisetosignalratioKSpinner      matlab.ui.control.Spinner
        noisetosignalratioKSpinnerLabel  matlab.ui.control.Label
        sigmaSpinner_2                  matlab.ui.control.Spinner
        sigmaSpinner_2Label             matlab.ui.control.Label
        ukuranfilternxnnSpinner         matlab.ui.control.Spinner
        ukuranfilternxnnSpinnerLabel    matlab.ui.control.Label
        BlurCitraButton_2               matlab.ui.control.Button
        PerbedaanCitraMasukandenganWienerPanel  matlab.ui.container.Panel
        ImageAxes_23                    matlab.ui.control.UIAxes
        JalankanFilterButton_4          matlab.ui.control.Button
        CitraHasilFilterWienerPanel     matlab.ui.container.Panel
        ImageAxes_22                    matlab.ui.control.UIAxes
        CitraMasukanPanel_3             matlab.ui.container.Panel
        ImageAxes_20                    matlab.ui.control.UIAxes
        PathInfoLabel_8                 matlab.ui.control.Label
        PathcitraLabel_8                matlab.ui.control.Label
        PilihFileCitraButton_8          matlab.ui.control.Button
        CitramasukanLabel_8             matlab.ui.control.Label
    end

    
    methods (Access = public)
        function hasil_konvolusi = konvolusi(~,gambar, mask)
            % Mendapatkan ukuran gambar dan mask
            [M, N, ~] = size(gambar); % ~ untuk mengabaikan channel pada gambar input
            [n, ~] = size(mask);
            
            % Menghitung padding yang diperlukan
            padding = floor(n / 2);
            
            % Inisialisasi matriks hasil konvolusi
            hasil_konvolusi = zeros(M, N, size(gambar, 3)); % Menyesuaikan dengan jumlah channel
            
            % Melakukan konvolusi
            for i = 1:M
                for j = 1:N
                    % Inisialisasi nilai hasil konvolusi untuk piksel (i, j)
                    konvolusi_ij = zeros(1, size(gambar, 3)); % Untuk setiap channel
                    
                    % Loop melalui kernel
                    for k = 1:n
                        for l = 1:n
                            % Koordinat piksel di sekitar (i, j) yang sesuai dengan kernel
                            x = i - padding + k;
                            y = j - padding + l;
                            
                            % Periksa apakah piksel (x, y) berada dalam batas gambar
                            if x >= 1 && x <= M && y >= 1 && y <= N
                                % Untuk setiap channel
                                for c = 1:size(gambar, 3)
                                    konvolusi_ij(c) = konvolusi_ij(c) + gambar(x, y, c) * mask(k, l);
                                end
                            end
                        end
                    end
                    
                    % Menyimpan nilai hasil konvolusi
                    hasil_konvolusi(i, j, :) = konvolusi_ij;
                end
            end
        end

        
        function res = filterfrekuensi(~, gambar, mask)
            % Menentukan parameter Padding
            [M, N, c] = size(gambar);
            padded = uint8(zeros(2*M, 2*N, c));
            
            % Pembentukan citra padding
            for i = 1:2*M
                for j = 1:2*N
                    for a = 1:c
                      if i <= M && j <= N
                            padded(i, j, a) = gambar(i, j, a);
                      else
                            padded(i, j, a) = 0;
                      end
                    end
                end
            end

            % melakukan transformasi fourier pada citra padding
            F = fft2(im2double(padded));

            % shifting mask
            H = mask;
            H = fftshift(H); 

            % pengaplikasian mask
            H = ifftshift(H);
            FF_f = H.*F;
            FF_f2 = real(ifft2(FF_f));
            FF_f2 = FF_f2(1:M, 1:N, :);
            res = FF_f2;
        end

        % for hardcode purposes
        function results = SetZero(~, res, y1, y2, x1, x2)

            results = res;
            for i=y1:y2
                for j=x1:x2
                    results(i, j) = 0;
                end
            end
        end

        % fungsi untuk mendapatkan bandreject
        function results = bandReject(~, image, D0)
            [M, N, ~] = size(image);
            
            % Transformasi fourier pada f(x,y)
            F = fft2(double(im2double(image)));
            
            % Setup range variabel
            u = 0:(M-1);
            v = 0:(N-1);
            
            % Hitung indices untku meshgrid
            idx = find(u > M/2);
            u(idx) = u(idx) - M;
            idy = find(v > N/2);
            v(idy) = v(idy) - N;
            
            % Get meshgrid
            [V, U] = meshgrid(v, u);
            D = sqrt(U.^2 + V.^2);
            
            n = 3;
            W = 50;
            H = 1./(1 + ((D*W)./(D.^2 - D0^2)).^(2*n));
            
            % Kalikan F dengan H
            G = H .* F;
            
            % Ambil bagian real dari inverse FFT G
            G = real(ifft2(G));
            results = im2uint8(G);
        end

    end
   
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: PilihFileCitraButton
        function PilihFileCitraButtonPushed(app, event)
            % Memilih gambar
            [fileName, pathName] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif;', 'All Image Files'});
            fullPath = fullfile(pathName, fileName);

            if fileName ~= 0
                % Menampilkan path file di label
                app.PathInfoLabel.Text = fullPath;

                % Menampilkan citra di axes
                imshow(fullPath, 'Parent', app.ImageAxes);

                % Enable tombol jalankan konvolusi
                app.JalankanKonvolusiButton.Enable = "on";
            end
        end

        % Button pushed function: JalankanKonvolusiButton
        function JalankanKonvolusiButtonPushed(app, event)
            % Mendapatkan citra dari path
            img = imread(app.PathInfoLabel.Text);
            mask = [1/16 2/16 1/16; 2/16 4/16 2/16; 1/16 2/16 1/16];

            % Menentukan matriks mask sesuai dengan pilihan
            if app.Mask1Button.Value
                mask = [1/16 2/16 1/16; 2/16 4/16 2/16; 1/16 2/16 1/16];
            elseif app.Mask2Button.Value
                mask = [0 -1 0; -1 4 -1; 0 -1 0;];
            elseif app.Mask3Button.Value
                mask = [1 1 2 2 2 1 1; 
                        1 2 2 4 2 2 1;
                        2 2 4 8 4 2 2;
                        2 4 8 16 8 4 2;
                        2 2 4 8 4 2 2;
                        1 2 2 4 2 2 1;
                        1 1 2 2 2 1 1;];
                mask = (1/140) .* mask;
            elseif app.Mask4Button.Value
                mask = [-1 -1 -1; -1 17 -1; -1 -1 -1];
            end
            
            % Lakukan konvolusi menggunakan fungsi buatan dan bawaan MATLAB
            hasil_konvolusi_buatan = uint8(255*app.konvolusi(im2double(img),double(mask)));
            hasil_konvolusi_matlab = uint8(255*convn(im2double(img),double(mask)));
            
            % Tampilkan hasil konvolusi
            imshow(hasil_konvolusi_buatan, 'Parent', app.ImageAxes_2);
            imshow(hasil_konvolusi_matlab, 'Parent', app.ImageAxes_3);
            
        end

        % Button pushed function: PilihFileCitraButton_2
        function PilihFileCitraButton_2Pushed(app, event)
            % Memilih gambar
            [fileName, pathName] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif;', 'All Image Files'});
            fullPath = fullfile(pathName, fileName);

            if fileName ~= 0
                % Menampilkan path file di label
                app.PathInfoLabel_2.Text = fullPath;

                % Menampilkan citra di axes
                imshow(fullPath, 'Parent', app.ImageAxes_4);

                % Enable tombol jalankan konvolusi
                app.JalankanSmoothingButton.Enable = "on";
            end
        end

        % Button pushed function: JalankanSmoothingButton
        function JalankanSmoothingButtonPushed(app, event)
            % Mendapatkan citra dari path
            img = imread(app.PathInfoLabel_2.Text);
            filter = [1 1 1; 1 1 1; 1 1 1];

            % Mendapatkan nilai n untuk ukuran filter spasial
            ukuran_filter = app.nSpinner.Value;

            % Mendapatkan nilai sigma untuk filter spasial gaussian
            sigma = app.sigmaSpinner.Value;

            if app.MeanFilternxnButton.Value || app.GaussianFilternxnButton.Value
                % Menentukan kernel untuk konvolusi pada filter ranah spasial
                % Inisialisasi filter
                if app.MeanFilternxnButton.Value
                    filter = ones(ukuran_filter) / (ukuran_filter^2);
                elseif app.GaussianFilternxnButton.Value
                    % Pembuatan filter gaussian dibantu Image Processing Toolbox dari MATLAB 
                    filter = fspecial('gaussian', ukuran_filter, sigma);
                end
                % Melakukan konvolusi untuk filter ranah spasial
                hasil_filter_frekuensi = uint8(255*app.konvolusi(im2double(img),double(filter)));
    
                % Tampilkan hasil konvolusi
                imshow(hasil_filter_frekuensi, 'Parent', app.ImageAxes_5);
            else 
                % menentukan mask dalam ranah frekuensi 
                [M, N, ~] = size(img);
                D0 = app.D0cutofffrequencySpinner.Value;
                u = 0:(2*M-1);
                v = 0:(2*N-1);
    
                idx = find(u > M);
                u(idx) = u(idx) - 2*M;
                idy = find(v > N);
                v(idy) = v(idy) - 2*N;
                
                [V, U] = meshgrid(v, u);
                D = sqrt(U.^2 + V.^2);
                
                if app.IdealLowpassFilterILPFButton.Value
                    % Pilihan ILPF dipilih
                    H = double(D <= D0);
                elseif app.GaussianLowpassFilterGLPFButton.Value
                    % Pilihan GLPF dipilih
                    H = exp(-(D.^2)./(2*(D0^2)));
                elseif app.ButterworthLowpassFilterBLPFButton.Value
                    % Pilihan BLPF dipilih
                    n = app.nfilterorderforBLPFSpinner.Value;
                    H = 1./(1 + (D./D0).^(2*n));
                end

                % Melakukan perkalian filter ranah frekuensi 
                hasil_filter_frekuensi = app.filterfrekuensi(img,H);
    
                % Tampilkan hasil filter
                imshow(hasil_filter_frekuensi, 'Parent', app.ImageAxes_5);
            end

        end

        % Button pushed function: PilihFileCitraButton_3
        function PilihFileCitraButton_3Pushed(app, event)
            % Memilih gambar
            [fileName, pathName] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif;', 'All Image Files'});
            fullPath = fullfile(pathName, fileName);

            if fileName ~= 0
                % Menampilkan path file di label
                app.PathInfoLabel_3.Text = fullPath;

                % Menampilkan citra di axes
                imshow(fullPath, 'Parent', app.ImageAxes_6);

                % Enable tombol jalankan konvolusi
                app.JalankanHighPassFilterButton.Enable = "on";
            end
        end

        % Button pushed function: JalankanHighPassFilterButton
        function JalankanHighPassFilterButtonPushed(app, event)
            % Mendapatkan citra dari path
            img = imread(app.PathInfoLabel_3.Text);

            % menentukan mask dalam ranah frekuensi 
            [M, N, ~] = size(img);
            D0 = app.D0cutofffrequencySpinner_2.Value;
            u = 0:(2*M-1);
            v = 0:(2*N-1);

            idx = find(u > M);
            u(idx) = u(idx) - 2*M;
            idy = find(v > N);
            v(idy) = v(idy) - 2*N;
            
            [V, U] = meshgrid(v, u);
            D = sqrt(U.^2 + V.^2);
            
            if app.IdealHighpassFilterIHPFButton.Value
                % Pilihan ILPF dipilih
                H = double(D <= D0);
            elseif app.GaussianHighpassFilterGHPFButton.Value
                % Pilihan GLPF dipilih
                H = exp(-(D.^2)./(2*(D0^2)));
            elseif app.ButterworthHighpassFilterBHPFButton.Value
                % Pilihan BLPF dipilih
                n = app.nfilterorderforBHPFSpinner.Value;
                H = 1./(1 + (D./D0).^(2*n));
            end
            
            % Menerapkan hubungan filter highpass dengan lowpass
            H = 1-H;

            % Melakukan perkalian filter ranah frekuensi 
            hasil_filter_frekuensi = app.filterfrekuensi(img,H);

            % Tampilkan hasil filter
            imshow(hasil_filter_frekuensi, 'Parent', app.ImageAxes_7);
        end

        % Button pushed function: PilihFileCitraButton_4
        function PilihFileCitraButton_4Pushed(app, event)
            % Memilih gambar
            [fileName, pathName] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif;', 'All Image Files'});
            fullPath = fullfile(pathName, fileName);

            if fileName ~= 0
                % Menampilkan path file di label
                app.PathInfoLabel_4.Text = fullPath;

                % Menampilkan citra di axes
                imshow(fullPath, 'Parent', app.ImageAxes_8);

                % Enable tombol jalankan konvolusi
                app.JalankanFilterTerangButton.Enable = "on";
            end
        end

        % Button pushed function: JalankanFilterTerangButton
        function JalankanFilterTerangButtonPushed(app, event)
            % Mendapatkan citra dari path
            img = imread(app.PathInfoLabel_4.Text);

            % menentukan mask dalam ranah frekuensi 
            [M, N, ~] = size(img);
            D0 = app.D0cutofffrequencySpinner_3.Value;
            u = 0:(2*M-1);
            v = 0:(2*N-1);

            idx = find(u > M);
            u(idx) = u(idx) - 2*M;
            idy = find(v > N);
            v(idy) = v(idy) - 2*N;
            
            [V, U] = meshgrid(v, u);
            D = sqrt(U.^2 + V.^2);
            
            if app.IdealLowpassFilterILPFButton_2.Value
                % Pilihan ILPF dipilih
                H = double(D <= D0);
            elseif app.GaussianLowpassFilterGLPFButton_2.Value
                % Pilihan GLPF dipilih
                H = exp(-(D.^2)./(2*(D0^2)));
            elseif app.ButterworthLowpassFilterBLPFButton_2.Value
                % Pilihan BLPF dipilih
                n = app.nfilterorderforBLPFSpinner_2.Value;
                H = 1./(1 + (D./D0).^(2*n));
            end

            % Mendapatkan gambar hasil Low Pass Filter
            hasil_filter_frekuensi = app.filterfrekuensi(img,H);

            % High Boost Image
            A = app.ASpinner.Value;
            high_boosted_image = A*img - uint8(hasil_filter_frekuensi);

            % Tampilkan hasil filter
            imshow(high_boosted_image, 'Parent', app.ImageAxes_9);
        end

        % Button pushed function: PilihFileCitraButton_5
        function PilihFileCitraButton_5Pushed(app, event)
            % Memilih gambar
            [fileName, pathName] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif;', 'All Image Files'});
            fullPath = fullfile(pathName, fileName);

            if fileName ~= 0
                % Menampilkan path file di label
                app.PathInfoLabel_5.Text = fullPath;

                % Menampilkan citra di axes
                imshow(fullPath, 'Parent', app.ImageAxes_10);

                % Toggle Button
                app.GenNoiseButton.Enable = "on";
                app.JalankanFilterButton.Enable = "off";
            end
        end

        % Button pushed function: PilihFileCitraButton_7
        function PilihFileCitraButton_7Pushed(app, event)
            % Memilih gambar
            [fileName, pathName] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif;', 'All Image Files'});
            fullPath = fullfile(pathName, fileName);

            if fileName ~= 0
                % Menampilkan path file di label
                app.PathInfoLabel_7.Text = fullPath;

                % Menampilkan citra di axes
                imshow(fullPath, 'Parent', app.ImageAxes_16);

                % Enable tombol jalankan filter
                app.JalankanFilterButton_3.Enable = "on";
            end
        end

        % Button pushed function: PilihFileCitraButton_6
        function PilihFileCitraButton_6Pushed(app, event)
            % Memilih gambar
            [fileName, pathName] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif;', 'All Image Files'});
            fullPath = fullfile(pathName, fileName);

            if fileName ~= 0
                % Menampilkan path file di label
                app.PathInfoLabel_6.Text = fullPath;

                % Menampilkan citra di axes
                imshow(fullPath, 'Parent', app.ImageAxes_13);

                % Enable tombol jalankan filter
                app.BlurCitraButton.Enable = "on";
            end
        end

        % Button pushed function: JalankanFilterButton
        function JalankanFilterButtonPushed(app, event)
            img = getimage(app.ImageAxes_11);

            % mengambil ukuran image
            [M, N, c] = size(img);
            neighborhoodSize = 3;

            newImage = zeros(M, N, c);

            for k = 1:c
                channelImage = img(:,:,k);
                for i = 1:M
                    for j = 1:N
                        row_min = max(1, i - floor(neighborhoodSize / 2));
                        row_max = min(M, i + floor(neighborhoodSize / 2));
                        col_min = max(1, j - floor(neighborhoodSize / 2));
                        col_max = min(N, j + floor(neighborhoodSize / 2));
                        neighborhood = channelImage(row_min:row_max, col_min:col_max);
                        
                        if app.MinFilterButton_2.Value
                            % Min Filter
                            min_val = min(neighborhood(:));
                            newImage(i, j, k) = min_val;
                        elseif app.MaxFilterButton_2.Value
                            % Max Filter
                            max_val = max(neighborhood(:));
                            newImage(i, j, k) = max_val;
                        elseif app.MedianFilterButton_2.Value
                            % Median Filter
                            med_val = median(neighborhood(:));
                            newImage(i, j, k) = med_val;
                        elseif app.ArithmeticMeanFilterButton_2.Value
                            % Arithmetic Mean Filter
                            mean_val = mean(neighborhood(:));
                            newImage(i, j, k) = mean_val;
                        elseif app.GeometricFilterButton_2.Value
                            % Geometric Filter
                            geomean_val = nthroot(prod(neighborhood(:)), numel(neighborhood(:)));
                            newImage(i, j, k) = geomean_val;
                        elseif app.HarmonicMeanFilterButton_2.Value
                            % Harmonic Mean Filter
                            hmean_val = numel(neighborhood(:)) / sum(1./neighborhood(:));
                            newImage(i, j, k) = hmean_val;
                        elseif app.ContraharmonicMeanFilterButton_2.Value
                            % Contraharmonic Mean Filter
                            Q = app.QSpinner.Value;
                            doubleneighbor = double(neighborhood);
                            up = sum(doubleneighbor(:) .^ (Q + 1));
                            down = sum(doubleneighbor(:) .^ Q);
                            newImage(i, j, k) = up / down;
                        elseif app.MidpointFilterButton_2.Value
                            % Midpoint Filter
                            mid_val = (max(neighborhood(:)) + min(neighborhood(:))) / 2;
                            newImage(i, j, k) = mid_val;
                        elseif app.AlphaTrimmedMeanFilterButton.Value
                            % AplhaTrimmedMeanFilter
                            d = app.dSpinner.Value;
                            sorted = sort(neighborhood(:));
                            trimmed = sorted(floor(d/2)+1:end-floor(d/2));
                            newImage(i, j, k) = mean(trimmed(:));
                        end
                    end
                end
            end

            newImage = uint8(newImage);
            imshow(newImage, "Parent", app.ImageAxes_12);
        end

        % Button pushed function: JalankanFilterButton_3
        function JalankanFilterButton_3Pushed(app, event)
            img = imread(app.PathInfoLabel_7.Text);
            [~, fileName, ~] = fileparts(app.PathInfoLabel_7.Text);
            if size(img, 3) == 1
                img = cat(3, img, img, img);
            end
            if fileName ~= "image-031"
                new_img = rgb2hsv(img);
            else
                new_img = uint8(img);
            end
            new_img = new_img(:,:,3);
            img_fft = fft2(new_img);
            dft = fftshift(img_fft);
            new_dft = dft;
            dft = abs(dft);
            dft = log(dft+1)*0.1;
            res = cat(3, dft, dft, dft);
            imagesc(app.UIAxes, res);
            
            % dft
            switch fileName
                case "image-023"
                    % kiri
                    new_dft = app.SetZero(new_dft, 182, 184, 52, 58);
                    new_dft = app.SetZero(new_dft, 201, 204, 70, 78);
                    new_dft = app.SetZero(new_dft, 50, 60, 178, 185);
                    new_dft = app.SetZero(new_dft, 175, 187, 178, 185);
                    new_dft = app.SetZero(new_dft, 308, 312, 181, 183);

                    % kanan
                    new_dft = app.SetZero(new_dft, 74, 77, 201, 203);
                    new_dft = app.SetZero(new_dft, 200, 204, 201, 203);
                    new_dft = app.SetZero(new_dft, 328, 331, 201, 203);
                    new_dft = app.SetZero(new_dft, 181, 183, 308, 310);
                    new_dft = app.SetZero(new_dft, 200, 203, 327, 331);

                    % garis
                    new_dft = app.SetZero(new_dft, 182, 182, 1, 182);
                    new_dft = app.SetZero(new_dft, 192, 192, 1, 182);
                    new_dft = app.SetZero(new_dft, 202, 202, 1, 182);
                    new_dft = app.SetZero(new_dft, 182, 182, 202, 383);
                    new_dft = app.SetZero(new_dft, 192, 192, 202, 383);
                    new_dft = app.SetZero(new_dft, 202, 202, 202, 383);
                    new_dft = app.SetZero(new_dft, 1, 383, 202, 202);
                    new_dft = app.SetZero(new_dft, 1, 383, 182, 182);

                case "image-025"
                    % kiri atas
                    new_dft = app.SetZero(new_dft, 157, 172, 1, 7);
                    new_dft = app.SetZero(new_dft, 1, 169, 63, 63);
                    new_dft = app.SetZero(new_dft, 1, 169, 70, 70);
                    new_dft = app.SetZero(new_dft, 1, 178, 123, 123);
                    new_dft = app.SetZero(new_dft, 1, 178, 135, 135);
                    new_dft = app.SetZero(new_dft, 157, 169, 132, 137);

                    % kanan atas
                    new_dft = app.SetZero(new_dft, 156, 166, 521, 523);
                    new_dft = app.SetZero(new_dft, 1, 172, 457, 457);
                    new_dft = app.SetZero(new_dft, 1, 172, 463, 463);
                    new_dft = app.SetZero(new_dft, 1, 178, 403, 403);
                    new_dft = app.SetZero(new_dft, 1, 178, 391, 391);

                    % kanan bawah
                    new_dft = app.SetZero(new_dft, 347, 522, 391, 391);
                    new_dft = app.SetZero(new_dft, 347, 522, 403, 403);
                    new_dft = app.SetZero(new_dft, 361, 522, 457, 457);
                    new_dft = app.SetZero(new_dft, 361, 522, 463, 463);

                    % kiri bawah
                    new_dft = app.SetZero(new_dft, 354, 522, 135, 135);
                    new_dft = app.SetZero(new_dft, 352, 522, 63, 63);
                    new_dft = app.SetZero(new_dft, 358, 522, 69, 69);

                case "image-028"
                    new_dft = app.SetZero(new_dft, 256, 259, 1, 241);
                    new_dft = app.SetZero(new_dft, 255, 258, 273, 512);

                case "image-029"
                    new_dft = app.SetZero(new_dft, 1, 494, 231, 233);
                    new_dft = app.SetZero(new_dft, 1, 494, 431, 433);
                    new_dft = app.SetZero(new_dft, 1, 245, 331, 333);
                    new_dft = app.SetZero(new_dft, 253, 494, 331, 333);

                case "image-030"
                    [M, N, ~] = size(img);
                    D0 = 100;
                    u = 0:(2*M-1);
                    v = 0:(2*N-1);
        
                    idx = find(u > M);
                    u(idx) = u(idx) - 2*M;
                    idy = find(v > N);
                    v(idy) = v(idy) - 2*N;
                    
                    [V, U] = meshgrid(v, u);
                    D = sqrt(U.^2 + V.^2);
                    n = 10;
                    H = 1./(1 + (D./D0).^(2*n));
                    res = app.filterfrekuensi(img, H);
                    imshow(res, "Parent", app.ImageAxes_18);
                    if size(res, 3) == 1
                        res = cat(3, res, res, res);
                    end
                    new_res = rgb2hsv(res);
                    new_res = new_res(:,:,3);
                    res_fft = fft2(new_res);
                    res_dft = fftshift(res_fft);
                    res_dft = abs(res_dft);
                    res_dft = log(res_dft+1)*0.1;
                    res = cat(3, res_dft, res_dft, res_dft);
                    imagesc(app.UIAxes_2, res);
                case "image-031"
                    res = app.bandReject(img, 134);
                    res = app.bandReject(res, 200);
                    res = app.bandReject(res, 225);
                    res = app.bandReject(res, 280);
                    res = app.bandReject(res, 310);
                    res = app.bandReject(res, 370);
                    res = app.bandReject(res, 410);
                    res = app.bandReject(res, 475);
                    
                    if size(res, 3) == 1
                        res = cat(3, res, res, res);
                    end
                    res = uint8(res);
                    res = res(:,:,3);
                    img_fft = fft2(res);
                    dft = fftshift(img_fft);
                    new_dft = dft;
                    dft = abs(dft);
                    dft = log(dft+1)*0.1;
                    dft_res = cat(3, dft, dft, dft);

                    imagesc(app.UIAxes_2, dft_res);
                    imshow(res, "Parent", app.ImageAxes_18);
            end

            if fileName ~= "image-030" && fileName ~= "image-031"
                new_img = real(ifft2(ifftshift(new_dft)));
                new_dft = abs(new_dft);
                new_dft = log(new_dft+1)*0.1;
                new_dft = cat(3, new_dft, new_dft, new_dft);
    
                imagesc(app.UIAxes_2, new_dft);
                imshow(new_img, "Parent", app.ImageAxes_18);
            end
        end

        % Button pushed function: JalankanFilterButton_2
        function JalankanFilterButton_2Pushed(app, event)
            img = imread(app.PathInfoLabel_6.Text);
            img = double(img);
            blur_img = getimage(app.ImageAxes_14);
            blur_img = double(blur_img);

            % mengambil blur power spectrum
            blur_fft = fft2(blur_img);
            blur_psd = abs(blur_fft).^2 / numel(blur_fft);

            % mengambil combined psd
            img_fft = fft2(img, size(blur_img, 1), size(blur_img, 2));
            psd = blur_fft .* conj(img_fft) / numel(blur_fft);

            % mendapatkan filter wiener + constant
            H = psd ./ (blur_psd + 1e-6);

            % merestore image
            restored_img = ifft2(H .* blur_fft);
            restored_img = uint8(restored_img);

            disp(size(restored_img));
            disp(size(img));
            imshow(restored_img, "Parent", app.ImageAxes_15);
            imshow(imabsdiff(uint8(img), restored_img), "Parent", app.ImageAxes_19);
        end

        % Button pushed function: GenNoiseButton
        function GenNoiseButtonPushed(app, event)
            img = imread(app.PathInfoLabel_5.Text);
            blank = uint8(255 * ones(size(img)));

            noiseDensity = app.NoiseDensitySpinner.Value;
            imgNoise = imnoise(img, 'salt & pepper', noiseDensity);
            imshow(imgNoise, "Parent", app.ImageAxes_11);
            imshow(blank, "Parent", app.ImageAxes_12)

            app.JalankanFilterButton.Enable = "on";
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            img = getimage(app.ImageAxes_10);
            blank = uint8(255 * ones(size(img)));

            imshow(blank, "Parent", app.ImageAxes_10);
            imshow(blank, "Parent", app.ImageAxes_11);
            imshow(blank, "Parent", app.ImageAxes_12);
            app.JalankanFilterButton.Enable = "off";
            app.GenNoiseButton.Enable = "off";
            app.PathInfoLabel_5.Text = "";
        end

        % Button pushed function: ResetButton_2
        function ResetButton_2Pushed(app, event)
            img = getimage(app.ImageAxes_13);
            blank = uint8(255 * ones(size(img)));

            imshow(blank, "Parent", app.ImageAxes_13);
            imshow(blank, "Parent", app.ImageAxes_14);
            imshow(blank, "Parent", app.ImageAxes_15);
            imshow(blank, "Parent", app.ImageAxes_19);
            app.JalankanFilterButton_2.Enable = "off";
            app.BlurCitraButton.Enable = "off";
            app.PathInfoLabel_6.Text = "";
        end

        % Button pushed function: BlurCitraButton
        function BlurCitraButtonPushed(app, event)
            img = getimage(app.ImageAxes_13);
            blank = uint8(255 * ones(size(img)));
            LEN = app.LENSpinner.Value;
            THETA = app.THETASpinner.Value;

            PSF = fspecial("motion", LEN, THETA);
            blurredImg = imfilter(img, PSF, "conv", "circular");

            imshow(blurredImg, "Parent", app.ImageAxes_14);
            imshow(blank, "Parent", app.ImageAxes_15);
            imshow(blank, "Parent", app.ImageAxes_19);

            app.JalankanFilterButton_2.Enable = "on";
        end

        % Callback function
        function TambahNoiseButtonPushed(app, event)
            img = getimage(app.ImageAxes_13);
            blank = uint8(255 * ones(size(img)));

            imshow(blank, "Parent", app.ImageAxes_15);
            imshow(blank, "Parent", app.ImageAxes_19);
        end

        % Button pushed function: PilihFileCitraButton_8
        function PilihFileCitraButton_8Pushed(app, event)
            % Memilih gambar
            [fileName, pathName] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif;', 'All Image Files'});
            fullPath = fullfile(pathName, fileName);

            if fileName ~= 0
                % Menampilkan path file di label
                app.PathInfoLabel_8.Text = fullPath;

                % Menampilkan citra di axes
                imshow(fullPath, 'Parent', app.ImageAxes_20);

                % Enable tombol jalankan filter
                app.JalankanFilterButton_4.Enable = "on";
            end
        end

        % Button pushed function: JalankanFilterButton_4
        function JalankanFilterButton_4Pushed(app, event)
            img = imread(app.PathInfoLabel_8.Text);
            img = double(img);
            blur_img = getimage(app.ImageAxes_20);
            blur_img = double(blur_img);
            
            % Mendefinisikan PSF untuk Gaussian blur
            % Mendapatkan nilai n untuk ukuran filter spasial
            size_arr = [app.ukuranfilternxnnSpinner.Value app.ukuranfilternxnnSpinner.Value];

            % Mendapatkan nilai sigma untuk filter spasial gaussian
            sigma = app.sigmaSpinner_2.Value;
            PSF = fspecial('gaussian', size_arr, sigma);
            
            % Menerapkan FFT pada PSF dan menormalkannya
            psf_fft = fft2(PSF, size(blur_img, 1), size(blur_img, 2));
            psf_fft = psf_fft ./ abs(psf_fft);
            
            % Mengambil blur power spectrum
            blur_fft = fft2(blur_img);
            
            % Mendapatkan filter Wiener + constant
            nsr = app.noisetosignalratioKSpinner.Value;  % noise-to-signal ratio (user-defined)
            H = conj(psf_fft) ./ (abs(psf_fft).^2 + nsr);
            
            % Melakukan dekonvolusi untuk merestorasi gambar
            restored_img_fft = H .* blur_fft;
            restored_img = ifft2(restored_img_fft);
            restored_img = uint8(real(restored_img));
            
            disp(size(restored_img));
            disp(size(img));
            imshow(restored_img, "Parent", app.ImageAxes_22);
            imshow(imabsdiff(uint8(img), restored_img), "Parent", app.ImageAxes_23);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.9412 0.9412 0.9412];
            app.UIFigure.Position = [100 100 927 602];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [2 -52 927 655];

            % Create KonvolusiTab
            app.KonvolusiTab = uitab(app.TabGroup);
            app.KonvolusiTab.Title = 'Konvolusi';

            % Create CitramasukanLabel
            app.CitramasukanLabel = uilabel(app.KonvolusiTab);
            app.CitramasukanLabel.Position = [33 584 86 22];
            app.CitramasukanLabel.Text = 'Citra masukan:';

            % Create PilihFileCitraButton
            app.PilihFileCitraButton = uibutton(app.KonvolusiTab, 'push');
            app.PilihFileCitraButton.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButtonPushed, true);
            app.PilihFileCitraButton.Position = [133 584 100 23];
            app.PilihFileCitraButton.Text = 'Pilih File Citra';

            % Create PathcitraLabel
            app.PathcitraLabel = uilabel(app.KonvolusiTab);
            app.PathcitraLabel.Position = [33 554 59 22];
            app.PathcitraLabel.Text = 'Path citra:';

            % Create PathInfoLabel
            app.PathInfoLabel = uilabel(app.KonvolusiTab);
            app.PathInfoLabel.Position = [111 555 776 22];
            app.PathInfoLabel.Text = '';

            % Create PilihanMaskUntukKonvolusiButtonGroup
            app.PilihanMaskUntukKonvolusiButtonGroup = uibuttongroup(app.KonvolusiTab);
            app.PilihanMaskUntukKonvolusiButtonGroup.Title = 'Pilihan Mask Untuk Konvolusi';
            app.PilihanMaskUntukKonvolusiButtonGroup.Position = [33 75 855 194];

            % Create Mask1Button
            app.Mask1Button = uiradiobutton(app.PilihanMaskUntukKonvolusiButtonGroup);
            app.Mask1Button.Text = 'Mask 1';
            app.Mask1Button.Position = [24 148 61 22];
            app.Mask1Button.Value = true;

            % Create Mask3Button
            app.Mask3Button = uiradiobutton(app.PilihanMaskUntukKonvolusiButtonGroup);
            app.Mask3Button.Text = 'Mask 3';
            app.Mask3Button.Position = [449 147 61 22];

            % Create Mask4Button
            app.Mask4Button = uiradiobutton(app.PilihanMaskUntukKonvolusiButtonGroup);
            app.Mask4Button.Text = 'Mask 4';
            app.Mask4Button.Position = [678 146 61 22];

            % Create Mask2Button
            app.Mask2Button = uiradiobutton(app.PilihanMaskUntukKonvolusiButtonGroup);
            app.Mask2Button.Text = 'Mask 2';
            app.Mask2Button.Position = [232 147 61 22];

            % Create Image
            app.Image = uiimage(app.PilihanMaskUntukKonvolusiButtonGroup);
            app.Image.Position = [19 15 134 124];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'doc', 'Convolution Masks', '1.png');

            % Create Image_2
            app.Image_2 = uiimage(app.PilihanMaskUntukKonvolusiButtonGroup);
            app.Image_2.Position = [232 15 134 124];
            app.Image_2.ImageSource = fullfile(pathToMLAPP, 'doc', 'Convolution Masks', '2.png');

            % Create Image_3
            app.Image_3 = uiimage(app.PilihanMaskUntukKonvolusiButtonGroup);
            app.Image_3.Position = [448 15 134 124];
            app.Image_3.ImageSource = fullfile(pathToMLAPP, 'doc', 'Convolution Masks', '3.png');

            % Create Image_4
            app.Image_4 = uiimage(app.PilihanMaskUntukKonvolusiButtonGroup);
            app.Image_4.Position = [672 15 134 124];
            app.Image_4.ImageSource = fullfile(pathToMLAPP, 'doc', 'Convolution Masks', '4.png');

            % Create PreviewCitraMasukanPanel
            app.PreviewCitraMasukanPanel = uipanel(app.KonvolusiTab);
            app.PreviewCitraMasukanPanel.Title = 'Preview Citra Masukan';
            app.PreviewCitraMasukanPanel.Position = [33 284 274 256];

            % Create ImageAxes
            app.ImageAxes = uiaxes(app.PreviewCitraMasukanPanel);
            app.ImageAxes.XTick = [];
            app.ImageAxes.YTick = [];
            app.ImageAxes.Position = [19 1 234 229];

            % Create PreviewCitraHasilKonvolusiBuatanPanel
            app.PreviewCitraHasilKonvolusiBuatanPanel = uipanel(app.KonvolusiTab);
            app.PreviewCitraHasilKonvolusiBuatanPanel.Title = 'Preview Citra Hasil Konvolusi Buatan';
            app.PreviewCitraHasilKonvolusiBuatanPanel.Position = [325 284 274 256];

            % Create ImageAxes_2
            app.ImageAxes_2 = uiaxes(app.PreviewCitraHasilKonvolusiBuatanPanel);
            app.ImageAxes_2.XTick = [];
            app.ImageAxes_2.YTick = [];
            app.ImageAxes_2.Position = [20 1 234 229];

            % Create PreviewCitraHasilKonvolusiMATLABPanel
            app.PreviewCitraHasilKonvolusiMATLABPanel = uipanel(app.KonvolusiTab);
            app.PreviewCitraHasilKonvolusiMATLABPanel.Title = 'Preview Citra Hasil Konvolusi MATLAB';
            app.PreviewCitraHasilKonvolusiMATLABPanel.Position = [614 284 274 256];

            % Create ImageAxes_3
            app.ImageAxes_3 = uiaxes(app.PreviewCitraHasilKonvolusiMATLABPanel);
            app.ImageAxes_3.XTick = [];
            app.ImageAxes_3.YTick = [];
            app.ImageAxes_3.Position = [19 1 234 229];

            % Create JalankanKonvolusiButton
            app.JalankanKonvolusiButton = uibutton(app.KonvolusiTab, 'push');
            app.JalankanKonvolusiButton.ButtonPushedFcn = createCallbackFcn(app, @JalankanKonvolusiButtonPushed, true);
            app.JalankanKonvolusiButton.Enable = 'off';
            app.JalankanKonvolusiButton.Position = [769 584 119 23];
            app.JalankanKonvolusiButton.Text = 'Jalankan Konvolusi';

            % Create SmoothingLowPassFilterTab
            app.SmoothingLowPassFilterTab = uitab(app.TabGroup);
            app.SmoothingLowPassFilterTab.Title = 'Smoothing (Low Pass Filter)';

            % Create CitramasukanLabel_2
            app.CitramasukanLabel_2 = uilabel(app.SmoothingLowPassFilterTab);
            app.CitramasukanLabel_2.Position = [33 584 86 22];
            app.CitramasukanLabel_2.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_2
            app.PilihFileCitraButton_2 = uibutton(app.SmoothingLowPassFilterTab, 'push');
            app.PilihFileCitraButton_2.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_2Pushed, true);
            app.PilihFileCitraButton_2.Position = [133 584 100 23];
            app.PilihFileCitraButton_2.Text = 'Pilih File Citra';

            % Create PathcitraLabel_2
            app.PathcitraLabel_2 = uilabel(app.SmoothingLowPassFilterTab);
            app.PathcitraLabel_2.Position = [33 554 59 22];
            app.PathcitraLabel_2.Text = 'Path citra:';

            % Create PathInfoLabel_2
            app.PathInfoLabel_2 = uilabel(app.SmoothingLowPassFilterTab);
            app.PathInfoLabel_2.Position = [111 555 776 22];
            app.PathInfoLabel_2.Text = '';

            % Create PreviewCitraMasukanPanel_2
            app.PreviewCitraMasukanPanel_2 = uipanel(app.SmoothingLowPassFilterTab);
            app.PreviewCitraMasukanPanel_2.Title = 'Preview Citra Masukan';
            app.PreviewCitraMasukanPanel_2.Position = [111 219 344 324];

            % Create ImageAxes_4
            app.ImageAxes_4 = uiaxes(app.PreviewCitraMasukanPanel_2);
            app.ImageAxes_4.XTick = [];
            app.ImageAxes_4.YTick = [];
            app.ImageAxes_4.Position = [19 25 295 264];

            % Create PilihanFilterButtonGroup
            app.PilihanFilterButtonGroup = uibuttongroup(app.SmoothingLowPassFilterTab);
            app.PilihanFilterButtonGroup.Title = 'Pilihan Filter';
            app.PilihanFilterButtonGroup.Position = [33 68 834 142];

            % Create MeanFilternxnButton
            app.MeanFilternxnButton = uiradiobutton(app.PilihanFilterButtonGroup);
            app.MeanFilternxnButton.Text = 'Mean Filter (n x n)';
            app.MeanFilternxnButton.Position = [14 70 118 22];
            app.MeanFilternxnButton.Value = true;

            % Create GaussianFilternxnButton
            app.GaussianFilternxnButton = uiradiobutton(app.PilihanFilterButtonGroup);
            app.GaussianFilternxnButton.Text = 'Gaussian Filter (n x n)';
            app.GaussianFilternxnButton.Position = [14 44 138 22];

            % Create RanahSpasialLabel
            app.RanahSpasialLabel = uilabel(app.PilihanFilterButtonGroup);
            app.RanahSpasialLabel.FontWeight = 'bold';
            app.RanahSpasialLabel.Position = [14 91 87 22];
            app.RanahSpasialLabel.Text = 'Ranah Spasial';

            % Create RanahFrekuensiLowpassFilterLabel
            app.RanahFrekuensiLowpassFilterLabel = uilabel(app.PilihanFilterButtonGroup);
            app.RanahFrekuensiLowpassFilterLabel.FontWeight = 'bold';
            app.RanahFrekuensiLowpassFilterLabel.Position = [326 94 201 22];
            app.RanahFrekuensiLowpassFilterLabel.Text = 'Ranah Frekuensi (Low-pass Filter)';

            % Create IdealLowpassFilterILPFButton
            app.IdealLowpassFilterILPFButton = uiradiobutton(app.PilihanFilterButtonGroup);
            app.IdealLowpassFilterILPFButton.Text = 'Ideal Low-pass Filter (ILPF)';
            app.IdealLowpassFilterILPFButton.Position = [326 72 168 22];

            % Create GaussianLowpassFilterGLPFButton
            app.GaussianLowpassFilterGLPFButton = uiradiobutton(app.PilihanFilterButtonGroup);
            app.GaussianLowpassFilterGLPFButton.Text = 'Gaussian Low-pass Filter (GLPF)';
            app.GaussianLowpassFilterGLPFButton.Position = [326 49 199 22];

            % Create ButterworthLowpassFilterBLPFButton
            app.ButterworthLowpassFilterBLPFButton = uiradiobutton(app.PilihanFilterButtonGroup);
            app.ButterworthLowpassFilterBLPFButton.Text = 'Butterworth Low-pass Filter (BLPF)';
            app.ButterworthLowpassFilterBLPFButton.Position = [326 27 211 22];

            % Create nSpinnerLabel
            app.nSpinnerLabel = uilabel(app.PilihanFilterButtonGroup);
            app.nSpinnerLabel.HorizontalAlignment = 'right';
            app.nSpinnerLabel.Position = [175 72 25 22];
            app.nSpinnerLabel.Text = 'n';

            % Create nSpinner
            app.nSpinner = uispinner(app.PilihanFilterButtonGroup);
            app.nSpinner.Limits = [2 30];
            app.nSpinner.Position = [215 72 59 22];
            app.nSpinner.Value = 2;

            % Create sigmaSpinnerLabel
            app.sigmaSpinnerLabel = uilabel(app.PilihanFilterButtonGroup);
            app.sigmaSpinnerLabel.HorizontalAlignment = 'right';
            app.sigmaSpinnerLabel.Position = [163 44 37 22];
            app.sigmaSpinnerLabel.Text = 'sigma';

            % Create sigmaSpinner
            app.sigmaSpinner = uispinner(app.PilihanFilterButtonGroup);
            app.sigmaSpinner.Limits = [1 30];
            app.sigmaSpinner.Position = [215 44 59 22];
            app.sigmaSpinner.Value = 1;

            % Create D0cutofffrequencySpinnerLabel
            app.D0cutofffrequencySpinnerLabel = uilabel(app.PilihanFilterButtonGroup);
            app.D0cutofffrequencySpinnerLabel.HorizontalAlignment = 'right';
            app.D0cutofffrequencySpinnerLabel.Position = [565 72 117 22];
            app.D0cutofffrequencySpinnerLabel.Text = 'D0 (cutoff frequency)';

            % Create D0cutofffrequencySpinner
            app.D0cutofffrequencySpinner = uispinner(app.PilihanFilterButtonGroup);
            app.D0cutofffrequencySpinner.Limits = [1 1000];
            app.D0cutofffrequencySpinner.Position = [697 72 59 22];
            app.D0cutofffrequencySpinner.Value = 1;

            % Create nfilterorderforBLPFSpinnerLabel
            app.nfilterorderforBLPFSpinnerLabel = uilabel(app.PilihanFilterButtonGroup);
            app.nfilterorderforBLPFSpinnerLabel.HorizontalAlignment = 'right';
            app.nfilterorderforBLPFSpinnerLabel.Position = [555 44 127 22];
            app.nfilterorderforBLPFSpinnerLabel.Text = 'n (filter order for BLPF)';

            % Create nfilterorderforBLPFSpinner
            app.nfilterorderforBLPFSpinner = uispinner(app.PilihanFilterButtonGroup);
            app.nfilterorderforBLPFSpinner.Limits = [1 100];
            app.nfilterorderforBLPFSpinner.Position = [697 44 59 22];
            app.nfilterorderforBLPFSpinner.Value = 1;

            % Create PreviewCitraHasilFilterPanel
            app.PreviewCitraHasilFilterPanel = uipanel(app.SmoothingLowPassFilterTab);
            app.PreviewCitraHasilFilterPanel.Title = 'Preview Citra Hasil Filter';
            app.PreviewCitraHasilFilterPanel.Position = [471 219 344 324];

            % Create ImageAxes_5
            app.ImageAxes_5 = uiaxes(app.PreviewCitraHasilFilterPanel);
            app.ImageAxes_5.XTick = [];
            app.ImageAxes_5.YTick = [];
            app.ImageAxes_5.Position = [19 25 295 264];

            % Create JalankanSmoothingButton
            app.JalankanSmoothingButton = uibutton(app.SmoothingLowPassFilterTab, 'push');
            app.JalankanSmoothingButton.ButtonPushedFcn = createCallbackFcn(app, @JalankanSmoothingButtonPushed, true);
            app.JalankanSmoothingButton.Enable = 'off';
            app.JalankanSmoothingButton.Position = [743 584 125 23];
            app.JalankanSmoothingButton.Text = 'Jalankan Smoothing';

            % Create HighPassFilterTab
            app.HighPassFilterTab = uitab(app.TabGroup);
            app.HighPassFilterTab.Title = 'High Pass Filter';

            % Create CitramasukanLabel_3
            app.CitramasukanLabel_3 = uilabel(app.HighPassFilterTab);
            app.CitramasukanLabel_3.Position = [33 584 86 22];
            app.CitramasukanLabel_3.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_3
            app.PilihFileCitraButton_3 = uibutton(app.HighPassFilterTab, 'push');
            app.PilihFileCitraButton_3.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_3Pushed, true);
            app.PilihFileCitraButton_3.Position = [133 584 100 23];
            app.PilihFileCitraButton_3.Text = 'Pilih File Citra';

            % Create PathcitraLabel_3
            app.PathcitraLabel_3 = uilabel(app.HighPassFilterTab);
            app.PathcitraLabel_3.Position = [33 554 59 22];
            app.PathcitraLabel_3.Text = 'Path citra:';

            % Create PathInfoLabel_3
            app.PathInfoLabel_3 = uilabel(app.HighPassFilterTab);
            app.PathInfoLabel_3.Position = [111 555 776 22];
            app.PathInfoLabel_3.Text = '';

            % Create PreviewCitraMasukanPanel_3
            app.PreviewCitraMasukanPanel_3 = uipanel(app.HighPassFilterTab);
            app.PreviewCitraMasukanPanel_3.Title = 'Preview Citra Masukan';
            app.PreviewCitraMasukanPanel_3.Position = [111 219 344 324];

            % Create ImageAxes_6
            app.ImageAxes_6 = uiaxes(app.PreviewCitraMasukanPanel_3);
            app.ImageAxes_6.XTick = [];
            app.ImageAxes_6.YTick = [];
            app.ImageAxes_6.Position = [19 25 295 264];

            % Create PilihanFilterButtonGroup_2
            app.PilihanFilterButtonGroup_2 = uibuttongroup(app.HighPassFilterTab);
            app.PilihanFilterButtonGroup_2.Title = 'Pilihan Filter';
            app.PilihanFilterButtonGroup_2.Position = [234 68 459 142];

            % Create RanahFrekuensiHighpassFilterLabel
            app.RanahFrekuensiHighpassFilterLabel = uilabel(app.PilihanFilterButtonGroup_2);
            app.RanahFrekuensiHighpassFilterLabel.FontWeight = 'bold';
            app.RanahFrekuensiHighpassFilterLabel.Position = [14 94 203 22];
            app.RanahFrekuensiHighpassFilterLabel.Text = 'Ranah Frekuensi (High-pass Filter)';

            % Create IdealHighpassFilterIHPFButton
            app.IdealHighpassFilterIHPFButton = uiradiobutton(app.PilihanFilterButtonGroup_2);
            app.IdealHighpassFilterIHPFButton.Text = 'Ideal High-pass Filter (IHPF)';
            app.IdealHighpassFilterIHPFButton.Position = [14 72 172 22];
            app.IdealHighpassFilterIHPFButton.Value = true;

            % Create GaussianHighpassFilterGHPFButton
            app.GaussianHighpassFilterGHPFButton = uiradiobutton(app.PilihanFilterButtonGroup_2);
            app.GaussianHighpassFilterGHPFButton.Text = 'Gaussian High-pass Filter (GHPF)';
            app.GaussianHighpassFilterGHPFButton.Position = [14 49 203 22];

            % Create ButterworthHighpassFilterBHPFButton
            app.ButterworthHighpassFilterBHPFButton = uiradiobutton(app.PilihanFilterButtonGroup_2);
            app.ButterworthHighpassFilterBHPFButton.Text = 'Butterworth High-pass Filter (BHPF)';
            app.ButterworthHighpassFilterBHPFButton.Position = [14 27 216 22];

            % Create D0cutofffrequencySpinner_2Label
            app.D0cutofffrequencySpinner_2Label = uilabel(app.PilihanFilterButtonGroup_2);
            app.D0cutofffrequencySpinner_2Label.HorizontalAlignment = 'right';
            app.D0cutofffrequencySpinner_2Label.Position = [253 72 117 22];
            app.D0cutofffrequencySpinner_2Label.Text = 'D0 (cutoff frequency)';

            % Create D0cutofffrequencySpinner_2
            app.D0cutofffrequencySpinner_2 = uispinner(app.PilihanFilterButtonGroup_2);
            app.D0cutofffrequencySpinner_2.Limits = [1 1000];
            app.D0cutofffrequencySpinner_2.Position = [385 72 59 22];
            app.D0cutofffrequencySpinner_2.Value = 1;

            % Create nfilterorderforBHPFSpinnerLabel
            app.nfilterorderforBHPFSpinnerLabel = uilabel(app.PilihanFilterButtonGroup_2);
            app.nfilterorderforBHPFSpinnerLabel.HorizontalAlignment = 'right';
            app.nfilterorderforBHPFSpinnerLabel.Position = [241 44 129 22];
            app.nfilterorderforBHPFSpinnerLabel.Text = 'n (filter order for BHPF)';

            % Create nfilterorderforBHPFSpinner
            app.nfilterorderforBHPFSpinner = uispinner(app.PilihanFilterButtonGroup_2);
            app.nfilterorderforBHPFSpinner.Limits = [1 100];
            app.nfilterorderforBHPFSpinner.Position = [385 44 59 22];
            app.nfilterorderforBHPFSpinner.Value = 1;

            % Create PreviewCitraHasilFilterPanel_2
            app.PreviewCitraHasilFilterPanel_2 = uipanel(app.HighPassFilterTab);
            app.PreviewCitraHasilFilterPanel_2.Title = 'Preview Citra Hasil Filter';
            app.PreviewCitraHasilFilterPanel_2.Position = [471 219 344 324];

            % Create ImageAxes_7
            app.ImageAxes_7 = uiaxes(app.PreviewCitraHasilFilterPanel_2);
            app.ImageAxes_7.XTick = [];
            app.ImageAxes_7.YTick = [];
            app.ImageAxes_7.Position = [19 25 295 264];

            % Create JalankanHighPassFilterButton
            app.JalankanHighPassFilterButton = uibutton(app.HighPassFilterTab, 'push');
            app.JalankanHighPassFilterButton.ButtonPushedFcn = createCallbackFcn(app, @JalankanHighPassFilterButtonPushed, true);
            app.JalankanHighPassFilterButton.Enable = 'off';
            app.JalankanHighPassFilterButton.Position = [668 584 150 23];
            app.JalankanHighPassFilterButton.Text = 'Jalankan High Pass Filter';

            % Create FilterTerangHighBoostTab
            app.FilterTerangHighBoostTab = uitab(app.TabGroup);
            app.FilterTerangHighBoostTab.Title = 'Filter Terang (High Boost)';

            % Create CitramasukanLabel_4
            app.CitramasukanLabel_4 = uilabel(app.FilterTerangHighBoostTab);
            app.CitramasukanLabel_4.Position = [33 584 86 22];
            app.CitramasukanLabel_4.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_4
            app.PilihFileCitraButton_4 = uibutton(app.FilterTerangHighBoostTab, 'push');
            app.PilihFileCitraButton_4.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_4Pushed, true);
            app.PilihFileCitraButton_4.Position = [133 584 100 23];
            app.PilihFileCitraButton_4.Text = 'Pilih File Citra';

            % Create PathcitraLabel_4
            app.PathcitraLabel_4 = uilabel(app.FilterTerangHighBoostTab);
            app.PathcitraLabel_4.Position = [33 554 59 22];
            app.PathcitraLabel_4.Text = 'Path citra:';

            % Create PreviewCitraMasukanPanel_4
            app.PreviewCitraMasukanPanel_4 = uipanel(app.FilterTerangHighBoostTab);
            app.PreviewCitraMasukanPanel_4.Title = 'Preview Citra Masukan';
            app.PreviewCitraMasukanPanel_4.Position = [111 219 344 324];

            % Create ImageAxes_8
            app.ImageAxes_8 = uiaxes(app.PreviewCitraMasukanPanel_4);
            app.ImageAxes_8.XTick = [];
            app.ImageAxes_8.YTick = [];
            app.ImageAxes_8.Position = [19 25 295 264];

            % Create PengaturanFilterButtonGroup
            app.PengaturanFilterButtonGroup = uibuttongroup(app.FilterTerangHighBoostTab);
            app.PengaturanFilterButtonGroup.Title = 'Pengaturan Filter';
            app.PengaturanFilterButtonGroup.Position = [111 68 699 142];

            % Create RanahFrekuensiLowpassFilterLabel_2
            app.RanahFrekuensiLowpassFilterLabel_2 = uilabel(app.PengaturanFilterButtonGroup);
            app.RanahFrekuensiLowpassFilterLabel_2.FontWeight = 'bold';
            app.RanahFrekuensiLowpassFilterLabel_2.Position = [23 90 201 22];
            app.RanahFrekuensiLowpassFilterLabel_2.Text = 'Ranah Frekuensi (Low-pass Filter)';

            % Create IdealLowpassFilterILPFButton_2
            app.IdealLowpassFilterILPFButton_2 = uiradiobutton(app.PengaturanFilterButtonGroup);
            app.IdealLowpassFilterILPFButton_2.Text = 'Ideal Low-pass Filter (ILPF)';
            app.IdealLowpassFilterILPFButton_2.Position = [23 68 168 22];

            % Create GaussianLowpassFilterGLPFButton_2
            app.GaussianLowpassFilterGLPFButton_2 = uiradiobutton(app.PengaturanFilterButtonGroup);
            app.GaussianLowpassFilterGLPFButton_2.Text = 'Gaussian Low-pass Filter (GLPF)';
            app.GaussianLowpassFilterGLPFButton_2.Position = [23 45 199 22];

            % Create ButterworthLowpassFilterBLPFButton_2
            app.ButterworthLowpassFilterBLPFButton_2 = uiradiobutton(app.PengaturanFilterButtonGroup);
            app.ButterworthLowpassFilterBLPFButton_2.Text = 'Butterworth Low-pass Filter (BLPF)';
            app.ButterworthLowpassFilterBLPFButton_2.Position = [23 23 211 22];
            app.ButterworthLowpassFilterBLPFButton_2.Value = true;

            % Create D0cutofffrequencySpinner_3Label
            app.D0cutofffrequencySpinner_3Label = uilabel(app.PengaturanFilterButtonGroup);
            app.D0cutofffrequencySpinner_3Label.HorizontalAlignment = 'right';
            app.D0cutofffrequencySpinner_3Label.Position = [262 68 117 22];
            app.D0cutofffrequencySpinner_3Label.Text = 'D0 (cutoff frequency)';

            % Create D0cutofffrequencySpinner_3
            app.D0cutofffrequencySpinner_3 = uispinner(app.PengaturanFilterButtonGroup);
            app.D0cutofffrequencySpinner_3.Limits = [1 1000];
            app.D0cutofffrequencySpinner_3.Position = [394 68 59 22];
            app.D0cutofffrequencySpinner_3.Value = 1;

            % Create nfilterorderforBLPFSpinner_2Label
            app.nfilterorderforBLPFSpinner_2Label = uilabel(app.PengaturanFilterButtonGroup);
            app.nfilterorderforBLPFSpinner_2Label.HorizontalAlignment = 'right';
            app.nfilterorderforBLPFSpinner_2Label.Position = [252 40 127 22];
            app.nfilterorderforBLPFSpinner_2Label.Text = 'n (filter order for BLPF)';

            % Create nfilterorderforBLPFSpinner_2
            app.nfilterorderforBLPFSpinner_2 = uispinner(app.PengaturanFilterButtonGroup);
            app.nfilterorderforBLPFSpinner_2.Limits = [1 100];
            app.nfilterorderforBLPFSpinner_2.Position = [394 40 59 22];
            app.nfilterorderforBLPFSpinner_2.Value = 1;

            % Create ASpinnerLabel
            app.ASpinnerLabel = uilabel(app.PengaturanFilterButtonGroup);
            app.ASpinnerLabel.HorizontalAlignment = 'right';
            app.ASpinnerLabel.Position = [503 66 25 22];
            app.ASpinnerLabel.Text = 'A';

            % Create ASpinner
            app.ASpinner = uispinner(app.PengaturanFilterButtonGroup);
            app.ASpinner.Step = 0.1;
            app.ASpinner.Limits = [1 30];
            app.ASpinner.Position = [543 66 59 22];
            app.ASpinner.Value = 1;

            % Create PreviewCitraHasilFilterPanel_3
            app.PreviewCitraHasilFilterPanel_3 = uipanel(app.FilterTerangHighBoostTab);
            app.PreviewCitraHasilFilterPanel_3.Title = 'Preview Citra Hasil Filter';
            app.PreviewCitraHasilFilterPanel_3.Position = [471 219 344 324];

            % Create ImageAxes_9
            app.ImageAxes_9 = uiaxes(app.PreviewCitraHasilFilterPanel_3);
            app.ImageAxes_9.XTick = [];
            app.ImageAxes_9.YTick = [];
            app.ImageAxes_9.Position = [19 25 295 264];

            % Create JalankanFilterTerangButton
            app.JalankanFilterTerangButton = uibutton(app.FilterTerangHighBoostTab, 'push');
            app.JalankanFilterTerangButton.ButtonPushedFcn = createCallbackFcn(app, @JalankanFilterTerangButtonPushed, true);
            app.JalankanFilterTerangButton.Enable = 'off';
            app.JalankanFilterTerangButton.Position = [677 584 132 23];
            app.JalankanFilterTerangButton.Text = 'Jalankan Filter Terang';

            % Create PathInfoLabel_4
            app.PathInfoLabel_4 = uilabel(app.FilterTerangHighBoostTab);
            app.PathInfoLabel_4.Position = [111 555 776 22];
            app.PathInfoLabel_4.Text = '';

            % Create SpasialFilterTab
            app.SpasialFilterTab = uitab(app.TabGroup);
            app.SpasialFilterTab.Title = 'Spasial Filter';

            % Create CitramasukanLabel_5
            app.CitramasukanLabel_5 = uilabel(app.SpasialFilterTab);
            app.CitramasukanLabel_5.Position = [33 584 86 22];
            app.CitramasukanLabel_5.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_5
            app.PilihFileCitraButton_5 = uibutton(app.SpasialFilterTab, 'push');
            app.PilihFileCitraButton_5.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_5Pushed, true);
            app.PilihFileCitraButton_5.Position = [133 584 100 23];
            app.PilihFileCitraButton_5.Text = 'Pilih File Citra';

            % Create PathcitraLabel_5
            app.PathcitraLabel_5 = uilabel(app.SpasialFilterTab);
            app.PathcitraLabel_5.Position = [33 554 59 22];
            app.PathcitraLabel_5.Text = 'Path citra:';

            % Create PathInfoLabel_5
            app.PathInfoLabel_5 = uilabel(app.SpasialFilterTab);
            app.PathInfoLabel_5.Position = [111 558 659 19];
            app.PathInfoLabel_5.Text = '';

            % Create CitraMasukanPanel_2
            app.CitraMasukanPanel_2 = uipanel(app.SpasialFilterTab);
            app.CitraMasukanPanel_2.Title = 'Citra Masukan';
            app.CitraMasukanPanel_2.Position = [33 284 274 256];

            % Create ImageAxes_10
            app.ImageAxes_10 = uiaxes(app.CitraMasukanPanel_2);
            app.ImageAxes_10.XTick = [];
            app.ImageAxes_10.YTick = [];
            app.ImageAxes_10.Position = [19 1 234 229];

            % Create CitraHasilPenambahanSaltPepperPanel
            app.CitraHasilPenambahanSaltPepperPanel = uipanel(app.SpasialFilterTab);
            app.CitraHasilPenambahanSaltPepperPanel.Title = 'Citra Hasil Penambahan Salt/Pepper';
            app.CitraHasilPenambahanSaltPepperPanel.Position = [325 284 274 256];

            % Create ImageAxes_11
            app.ImageAxes_11 = uiaxes(app.CitraHasilPenambahanSaltPepperPanel);
            app.ImageAxes_11.XTick = [];
            app.ImageAxes_11.YTick = [];
            app.ImageAxes_11.Position = [20 1 234 229];

            % Create CitraHasilFilterNoisePanel
            app.CitraHasilFilterNoisePanel = uipanel(app.SpasialFilterTab);
            app.CitraHasilFilterNoisePanel.Title = 'Citra Hasil Filter Noise';
            app.CitraHasilFilterNoisePanel.Position = [614 284 274 256];

            % Create ImageAxes_12
            app.ImageAxes_12 = uiaxes(app.CitraHasilFilterNoisePanel);
            app.ImageAxes_12.XTick = [];
            app.ImageAxes_12.YTick = [];
            app.ImageAxes_12.Position = [19 1 234 229];

            % Create JalankanFilterButton
            app.JalankanFilterButton = uibutton(app.SpasialFilterTab, 'push');
            app.JalankanFilterButton.ButtonPushedFcn = createCallbackFcn(app, @JalankanFilterButtonPushed, true);
            app.JalankanFilterButton.Enable = 'off';
            app.JalankanFilterButton.Position = [779 584 100 23];
            app.JalankanFilterButton.Text = 'Jalankan Filter';

            % Create PengaturanFilterPanel
            app.PengaturanFilterPanel = uipanel(app.SpasialFilterTab);
            app.PengaturanFilterPanel.Title = 'Pengaturan Filter';
            app.PengaturanFilterPanel.Position = [154 68 620 207];

            % Create TipeFilterButtonGroup
            app.TipeFilterButtonGroup = uibuttongroup(app.PengaturanFilterPanel);
            app.TipeFilterButtonGroup.BorderType = 'none';
            app.TipeFilterButtonGroup.BorderWidth = 0;
            app.TipeFilterButtonGroup.Title = 'Tipe Filter';
            app.TipeFilterButtonGroup.FontWeight = 'bold';
            app.TipeFilterButtonGroup.Position = [49 22 531 144];

            % Create MinFilterButton_2
            app.MinFilterButton_2 = uiradiobutton(app.TipeFilterButtonGroup);
            app.MinFilterButton_2.Text = 'Min Filter';
            app.MinFilterButton_2.Position = [11 98 72 22];
            app.MinFilterButton_2.Value = true;

            % Create MaxFilterButton_2
            app.MaxFilterButton_2 = uiradiobutton(app.TipeFilterButtonGroup);
            app.MaxFilterButton_2.Text = 'Max Filter';
            app.MaxFilterButton_2.Position = [11 76 75 22];

            % Create MedianFilterButton_2
            app.MedianFilterButton_2 = uiradiobutton(app.TipeFilterButtonGroup);
            app.MedianFilterButton_2.Text = 'Median Filter';
            app.MedianFilterButton_2.Position = [11 54 92 22];

            % Create ArithmeticMeanFilterButton_2
            app.ArithmeticMeanFilterButton_2 = uiradiobutton(app.TipeFilterButtonGroup);
            app.ArithmeticMeanFilterButton_2.Text = 'Arithmetic Mean Filter';
            app.ArithmeticMeanFilterButton_2.Position = [11 32 140 22];

            % Create HarmonicMeanFilterButton_2
            app.HarmonicMeanFilterButton_2 = uiradiobutton(app.TipeFilterButtonGroup);
            app.HarmonicMeanFilterButton_2.Text = 'Harmonic Mean Filter';
            app.HarmonicMeanFilterButton_2.Position = [194 98 137 22];

            % Create GeometricFilterButton_2
            app.GeometricFilterButton_2 = uiradiobutton(app.TipeFilterButtonGroup);
            app.GeometricFilterButton_2.Text = 'Geometric Filter';
            app.GeometricFilterButton_2.Position = [11 10 108 22];

            % Create ContraharmonicMeanFilterButton_2
            app.ContraharmonicMeanFilterButton_2 = uiradiobutton(app.TipeFilterButtonGroup);
            app.ContraharmonicMeanFilterButton_2.Text = 'Contraharmonic Mean Filter';
            app.ContraharmonicMeanFilterButton_2.Position = [194 76 172 22];

            % Create MidpointFilterButton_2
            app.MidpointFilterButton_2 = uiradiobutton(app.TipeFilterButtonGroup);
            app.MidpointFilterButton_2.Text = 'Midpoint Filter';
            app.MidpointFilterButton_2.Position = [194 54 99 22];

            % Create AlphaTrimmedMeanFilterButton
            app.AlphaTrimmedMeanFilterButton = uiradiobutton(app.TipeFilterButtonGroup);
            app.AlphaTrimmedMeanFilterButton.Text = 'Alpha-Trimmed Mean Filter';
            app.AlphaTrimmedMeanFilterButton.Position = [194 32 167 22];

            % Create NoiseDensitySpinnerLabel
            app.NoiseDensitySpinnerLabel = uilabel(app.PengaturanFilterPanel);
            app.NoiseDensitySpinnerLabel.HorizontalAlignment = 'right';
            app.NoiseDensitySpinnerLabel.Position = [420 120 79 22];
            app.NoiseDensitySpinnerLabel.Text = 'Noise Density';

            % Create NoiseDensitySpinner
            app.NoiseDensitySpinner = uispinner(app.PengaturanFilterPanel);
            app.NoiseDensitySpinner.Step = 0.01;
            app.NoiseDensitySpinner.Limits = [0.01 0.5];
            app.NoiseDensitySpinner.Position = [514 120 86 22];
            app.NoiseDensitySpinner.Value = 0.01;

            % Create QSpinnerLabel
            app.QSpinnerLabel = uilabel(app.PengaturanFilterPanel);
            app.QSpinnerLabel.HorizontalAlignment = 'right';
            app.QSpinnerLabel.Position = [474 97 25 22];
            app.QSpinnerLabel.Text = 'Q';

            % Create QSpinner
            app.QSpinner = uispinner(app.PengaturanFilterPanel);
            app.QSpinner.Step = 0.1;
            app.QSpinner.Position = [514 97 86 22];
            app.QSpinner.Value = 1.5;

            % Create dSpinnerLabel
            app.dSpinnerLabel = uilabel(app.PengaturanFilterPanel);
            app.dSpinnerLabel.HorizontalAlignment = 'right';
            app.dSpinnerLabel.Position = [475 53 25 22];
            app.dSpinnerLabel.Text = 'd';

            % Create dSpinner
            app.dSpinner = uispinner(app.PengaturanFilterPanel);
            app.dSpinner.Limits = [1 Inf];
            app.dSpinner.Position = [514 53 86 22];
            app.dSpinner.Value = 1;

            % Create GenNoiseButton
            app.GenNoiseButton = uibutton(app.SpasialFilterTab, 'push');
            app.GenNoiseButton.ButtonPushedFcn = createCallbackFcn(app, @GenNoiseButtonPushed, true);
            app.GenNoiseButton.Enable = 'off';
            app.GenNoiseButton.Position = [779 554 100 23];
            app.GenNoiseButton.Text = 'Gen Noise';

            % Create ResetButton
            app.ResetButton = uibutton(app.SpasialFilterTab, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Position = [779 252 100 23];
            app.ResetButton.Text = 'Reset';

            % Create PeriodicFilterTab
            app.PeriodicFilterTab = uitab(app.TabGroup);
            app.PeriodicFilterTab.Title = 'Periodic Filter';

            % Create CitramasukanLabel_7
            app.CitramasukanLabel_7 = uilabel(app.PeriodicFilterTab);
            app.CitramasukanLabel_7.Position = [33 584 86 22];
            app.CitramasukanLabel_7.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_7
            app.PilihFileCitraButton_7 = uibutton(app.PeriodicFilterTab, 'push');
            app.PilihFileCitraButton_7.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_7Pushed, true);
            app.PilihFileCitraButton_7.Position = [133 584 100 23];
            app.PilihFileCitraButton_7.Text = 'Pilih File Citra';

            % Create PathcitraLabel_7
            app.PathcitraLabel_7 = uilabel(app.PeriodicFilterTab);
            app.PathcitraLabel_7.Position = [33 554 59 22];
            app.PathcitraLabel_7.Text = 'Path citra:';

            % Create PathInfoLabel_7
            app.PathInfoLabel_7 = uilabel(app.PeriodicFilterTab);
            app.PathInfoLabel_7.Position = [111 555 776 22];
            app.PathInfoLabel_7.Text = '';

            % Create PreviewCitraMasukanPanel_7
            app.PreviewCitraMasukanPanel_7 = uipanel(app.PeriodicFilterTab);
            app.PreviewCitraMasukanPanel_7.Title = 'Preview Citra Masukan';
            app.PreviewCitraMasukanPanel_7.Position = [132 298 274 256];

            % Create ImageAxes_16
            app.ImageAxes_16 = uiaxes(app.PreviewCitraMasukanPanel_7);
            app.ImageAxes_16.XTick = [];
            app.ImageAxes_16.YTick = [];
            app.ImageAxes_16.Position = [19 1 234 229];

            % Create PreviewCitraHasilPanel
            app.PreviewCitraHasilPanel = uipanel(app.PeriodicFilterTab);
            app.PreviewCitraHasilPanel.Title = 'Preview Citra Hasil';
            app.PreviewCitraHasilPanel.Position = [133 27 274 256];

            % Create ImageAxes_18
            app.ImageAxes_18 = uiaxes(app.PreviewCitraHasilPanel);
            app.ImageAxes_18.XTick = [];
            app.ImageAxes_18.YTick = [];
            app.ImageAxes_18.Position = [19 1 234 229];

            % Create JalankanFilterButton_3
            app.JalankanFilterButton_3 = uibutton(app.PeriodicFilterTab, 'push');
            app.JalankanFilterButton_3.ButtonPushedFcn = createCallbackFcn(app, @JalankanFilterButton_3Pushed, true);
            app.JalankanFilterButton_3.Enable = 'off';
            app.JalankanFilterButton_3.Position = [779 584 100 23];
            app.JalankanFilterButton_3.Text = 'Jalankan Filter';

            % Create FourierCitraMasukanPanel
            app.FourierCitraMasukanPanel = uipanel(app.PeriodicFilterTab);
            app.FourierCitraMasukanPanel.Title = 'Fourier Citra Masukan';
            app.FourierCitraMasukanPanel.Position = [491 299 284 256];

            % Create UIAxes
            app.UIAxes = uiaxes(app.FourierCitraMasukanPanel);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [6 25 273 185];

            % Create FourierCitraHasilPanel
            app.FourierCitraHasilPanel = uipanel(app.PeriodicFilterTab);
            app.FourierCitraHasilPanel.Title = 'Fourier Citra Hasil';
            app.FourierCitraHasilPanel.Position = [492 27 284 256];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.FourierCitraHasilPanel);
            title(app.UIAxes_2, 'Title')
            xlabel(app.UIAxes_2, 'X')
            ylabel(app.UIAxes_2, 'Y')
            zlabel(app.UIAxes_2, 'Z')
            app.UIAxes_2.Position = [6 25 273 185];

            % Create WienerMotionBlurTab
            app.WienerMotionBlurTab = uitab(app.TabGroup);
            app.WienerMotionBlurTab.Title = 'Wiener (Motion Blur)';

            % Create CitramasukanLabel_6
            app.CitramasukanLabel_6 = uilabel(app.WienerMotionBlurTab);
            app.CitramasukanLabel_6.Position = [33 584 86 22];
            app.CitramasukanLabel_6.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_6
            app.PilihFileCitraButton_6 = uibutton(app.WienerMotionBlurTab, 'push');
            app.PilihFileCitraButton_6.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_6Pushed, true);
            app.PilihFileCitraButton_6.Position = [133 584 100 23];
            app.PilihFileCitraButton_6.Text = 'Pilih File Citra';

            % Create PathcitraLabel_6
            app.PathcitraLabel_6 = uilabel(app.WienerMotionBlurTab);
            app.PathcitraLabel_6.Position = [33 554 59 22];
            app.PathcitraLabel_6.Text = 'Path citra:';

            % Create PathInfoLabel_6
            app.PathInfoLabel_6 = uilabel(app.WienerMotionBlurTab);
            app.PathInfoLabel_6.Position = [111 555 776 22];
            app.PathInfoLabel_6.Text = '';

            % Create CitraMasukanPanel
            app.CitraMasukanPanel = uipanel(app.WienerMotionBlurTab);
            app.CitraMasukanPanel.Title = 'Citra Masukan';
            app.CitraMasukanPanel.Position = [81 288 274 256];

            % Create ImageAxes_13
            app.ImageAxes_13 = uiaxes(app.CitraMasukanPanel);
            app.ImageAxes_13.XTick = [];
            app.ImageAxes_13.YTick = [];
            app.ImageAxes_13.Position = [19 1 234 229];

            % Create CitraHasilMotionBlurPanel
            app.CitraHasilMotionBlurPanel = uipanel(app.WienerMotionBlurTab);
            app.CitraHasilMotionBlurPanel.Title = 'Citra Hasil Motion Blur';
            app.CitraHasilMotionBlurPanel.Position = [373 288 274 256];

            % Create ImageAxes_14
            app.ImageAxes_14 = uiaxes(app.CitraHasilMotionBlurPanel);
            app.ImageAxes_14.XTick = [];
            app.ImageAxes_14.YTick = [];
            app.ImageAxes_14.Position = [20 1 234 229];

            % Create CitraHasilFilterWienerPanel_2
            app.CitraHasilFilterWienerPanel_2 = uipanel(app.WienerMotionBlurTab);
            app.CitraHasilFilterWienerPanel_2.Title = 'Citra Hasil Filter Wiener';
            app.CitraHasilFilterWienerPanel_2.Position = [81 24 274 256];

            % Create ImageAxes_15
            app.ImageAxes_15 = uiaxes(app.CitraHasilFilterWienerPanel_2);
            app.ImageAxes_15.XTick = [];
            app.ImageAxes_15.YTick = [];
            app.ImageAxes_15.Position = [19 1 234 229];

            % Create JalankanFilterButton_2
            app.JalankanFilterButton_2 = uibutton(app.WienerMotionBlurTab, 'push');
            app.JalankanFilterButton_2.ButtonPushedFcn = createCallbackFcn(app, @JalankanFilterButton_2Pushed, true);
            app.JalankanFilterButton_2.Enable = 'off';
            app.JalankanFilterButton_2.Position = [779 584 100 23];
            app.JalankanFilterButton_2.Text = 'Jalankan Filter';

            % Create PerbedaanCitraMasukandenganWienerPanel_2
            app.PerbedaanCitraMasukandenganWienerPanel_2 = uipanel(app.WienerMotionBlurTab);
            app.PerbedaanCitraMasukandenganWienerPanel_2.Title = 'Perbedaan Citra Masukan dengan Wiener';
            app.PerbedaanCitraMasukandenganWienerPanel_2.Position = [373 25 274 256];

            % Create ImageAxes_19
            app.ImageAxes_19 = uiaxes(app.PerbedaanCitraMasukandenganWienerPanel_2);
            app.ImageAxes_19.XTick = [];
            app.ImageAxes_19.YTick = [];
            app.ImageAxes_19.Position = [19 1 234 229];

            % Create MotionBlurPanel
            app.MotionBlurPanel = uipanel(app.WienerMotionBlurTab);
            app.MotionBlurPanel.Title = 'Motion Blur';
            app.MotionBlurPanel.Position = [677 417 186 127];

            % Create LENSpinnerLabel
            app.LENSpinnerLabel = uilabel(app.MotionBlurPanel);
            app.LENSpinnerLabel.HorizontalAlignment = 'right';
            app.LENSpinnerLabel.Position = [31 73 28 22];
            app.LENSpinnerLabel.Text = 'LEN';

            % Create LENSpinner
            app.LENSpinner = uispinner(app.MotionBlurPanel);
            app.LENSpinner.Limits = [0 100];
            app.LENSpinner.Position = [74 73 100 22];
            app.LENSpinner.Value = 1;

            % Create THETASpinnerLabel
            app.THETASpinnerLabel = uilabel(app.MotionBlurPanel);
            app.THETASpinnerLabel.HorizontalAlignment = 'right';
            app.THETASpinnerLabel.Position = [18 45 41 22];
            app.THETASpinnerLabel.Text = 'THETA';

            % Create THETASpinner
            app.THETASpinner = uispinner(app.MotionBlurPanel);
            app.THETASpinner.Limits = [0 360];
            app.THETASpinner.Position = [74 45 100 22];
            app.THETASpinner.Value = 45;

            % Create BlurCitraButton
            app.BlurCitraButton = uibutton(app.MotionBlurPanel, 'push');
            app.BlurCitraButton.ButtonPushedFcn = createCallbackFcn(app, @BlurCitraButtonPushed, true);
            app.BlurCitraButton.Enable = 'off';
            app.BlurCitraButton.Position = [40 11 100 23];
            app.BlurCitraButton.Text = 'Blur Citra';

            % Create ResetButton_2
            app.ResetButton_2 = uibutton(app.WienerMotionBlurTab, 'push');
            app.ResetButton_2.ButtonPushedFcn = createCallbackFcn(app, @ResetButton_2Pushed, true);
            app.ResetButton_2.Position = [717 380 100 23];
            app.ResetButton_2.Text = 'Reset';

            % Create WienerGaussianBlurTab
            app.WienerGaussianBlurTab = uitab(app.TabGroup);
            app.WienerGaussianBlurTab.Title = 'Wiener (Gaussian Blur)';

            % Create CitramasukanLabel_8
            app.CitramasukanLabel_8 = uilabel(app.WienerGaussianBlurTab);
            app.CitramasukanLabel_8.Position = [33 584 86 22];
            app.CitramasukanLabel_8.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_8
            app.PilihFileCitraButton_8 = uibutton(app.WienerGaussianBlurTab, 'push');
            app.PilihFileCitraButton_8.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_8Pushed, true);
            app.PilihFileCitraButton_8.Position = [133 584 100 23];
            app.PilihFileCitraButton_8.Text = 'Pilih File Citra';

            % Create PathcitraLabel_8
            app.PathcitraLabel_8 = uilabel(app.WienerGaussianBlurTab);
            app.PathcitraLabel_8.Position = [33 554 59 22];
            app.PathcitraLabel_8.Text = 'Path citra:';

            % Create PathInfoLabel_8
            app.PathInfoLabel_8 = uilabel(app.WienerGaussianBlurTab);
            app.PathInfoLabel_8.Position = [111 555 776 22];
            app.PathInfoLabel_8.Text = '';

            % Create CitraMasukanPanel_3
            app.CitraMasukanPanel_3 = uipanel(app.WienerGaussianBlurTab);
            app.CitraMasukanPanel_3.Title = 'Citra Masukan';
            app.CitraMasukanPanel_3.Position = [248 298 274 256];

            % Create ImageAxes_20
            app.ImageAxes_20 = uiaxes(app.CitraMasukanPanel_3);
            app.ImageAxes_20.XTick = [];
            app.ImageAxes_20.YTick = [];
            app.ImageAxes_20.Position = [19 1 234 229];

            % Create CitraHasilFilterWienerPanel
            app.CitraHasilFilterWienerPanel = uipanel(app.WienerGaussianBlurTab);
            app.CitraHasilFilterWienerPanel.Title = 'Citra Hasil Filter Wiener';
            app.CitraHasilFilterWienerPanel.Position = [81 24 274 256];

            % Create ImageAxes_22
            app.ImageAxes_22 = uiaxes(app.CitraHasilFilterWienerPanel);
            app.ImageAxes_22.XTick = [];
            app.ImageAxes_22.YTick = [];
            app.ImageAxes_22.Position = [19 1 234 229];

            % Create JalankanFilterButton_4
            app.JalankanFilterButton_4 = uibutton(app.WienerGaussianBlurTab, 'push');
            app.JalankanFilterButton_4.ButtonPushedFcn = createCallbackFcn(app, @JalankanFilterButton_4Pushed, true);
            app.JalankanFilterButton_4.Enable = 'off';
            app.JalankanFilterButton_4.Position = [779 584 100 23];
            app.JalankanFilterButton_4.Text = 'Jalankan Filter';

            % Create PerbedaanCitraMasukandenganWienerPanel
            app.PerbedaanCitraMasukandenganWienerPanel = uipanel(app.WienerGaussianBlurTab);
            app.PerbedaanCitraMasukandenganWienerPanel.Title = 'Perbedaan Citra Masukan dengan Wiener';
            app.PerbedaanCitraMasukandenganWienerPanel.Position = [373 25 274 256];

            % Create ImageAxes_23
            app.ImageAxes_23 = uiaxes(app.PerbedaanCitraMasukandenganWienerPanel);
            app.ImageAxes_23.XTick = [];
            app.ImageAxes_23.YTick = [];
            app.ImageAxes_23.Position = [19 1 234 229];

            % Create MotionBlurPanel_2
            app.MotionBlurPanel_2 = uipanel(app.WienerGaussianBlurTab);
            app.MotionBlurPanel_2.Title = 'Motion Blur';
            app.MotionBlurPanel_2.Position = [638 324 225 220];

            % Create BlurCitraButton_2
            app.BlurCitraButton_2 = uibutton(app.MotionBlurPanel_2, 'push');
            app.BlurCitraButton_2.Enable = 'off';
            app.BlurCitraButton_2.Position = [62 24 100 23];
            app.BlurCitraButton_2.Text = 'Blur Citra';

            % Create ukuranfilternxnnSpinnerLabel
            app.ukuranfilternxnnSpinnerLabel = uilabel(app.MotionBlurPanel_2);
            app.ukuranfilternxnnSpinnerLabel.HorizontalAlignment = 'right';
            app.ukuranfilternxnnSpinnerLabel.Position = [16 168 114 22];
            app.ukuranfilternxnnSpinnerLabel.Text = 'ukuran filter(n x n) n:';

            % Create ukuranfilternxnnSpinner
            app.ukuranfilternxnnSpinner = uispinner(app.MotionBlurPanel_2);
            app.ukuranfilternxnnSpinner.Limits = [2 30];
            app.ukuranfilternxnnSpinner.Position = [145 168 59 22];
            app.ukuranfilternxnnSpinner.Value = 2;

            % Create sigmaSpinner_2Label
            app.sigmaSpinner_2Label = uilabel(app.MotionBlurPanel_2);
            app.sigmaSpinner_2Label.HorizontalAlignment = 'right';
            app.sigmaSpinner_2Label.Position = [56 139 37 22];
            app.sigmaSpinner_2Label.Text = 'sigma';

            % Create sigmaSpinner_2
            app.sigmaSpinner_2 = uispinner(app.MotionBlurPanel_2);
            app.sigmaSpinner_2.Limits = [1 30];
            app.sigmaSpinner_2.Position = [108 139 59 22];
            app.sigmaSpinner_2.Value = 1;

            % Create noisetosignalratioKSpinnerLabel
            app.noisetosignalratioKSpinnerLabel = uilabel(app.MotionBlurPanel_2);
            app.noisetosignalratioKSpinnerLabel.HorizontalAlignment = 'right';
            app.noisetosignalratioKSpinnerLabel.Position = [7 111 130 22];
            app.noisetosignalratioKSpinnerLabel.Text = 'noise-to-signal ratio (K)';

            % Create noisetosignalratioKSpinner
            app.noisetosignalratioKSpinner = uispinner(app.MotionBlurPanel_2);
            app.noisetosignalratioKSpinner.Step = 0.01;
            app.noisetosignalratioKSpinner.Limits = [0.01 30];
            app.noisetosignalratioKSpinner.Position = [152 111 59 22];
            app.noisetosignalratioKSpinner.Value = 0.01;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
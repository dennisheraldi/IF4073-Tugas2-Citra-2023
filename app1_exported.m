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
            app.TabGroup.Position = [2 1 927 602];

            % Create KonvolusiTab
            app.KonvolusiTab = uitab(app.TabGroup);
            app.KonvolusiTab.Title = 'Konvolusi';

            % Create CitramasukanLabel
            app.CitramasukanLabel = uilabel(app.KonvolusiTab);
            app.CitramasukanLabel.Position = [33 531 86 22];
            app.CitramasukanLabel.Text = 'Citra masukan:';

            % Create PilihFileCitraButton
            app.PilihFileCitraButton = uibutton(app.KonvolusiTab, 'push');
            app.PilihFileCitraButton.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButtonPushed, true);
            app.PilihFileCitraButton.Position = [133 531 100 23];
            app.PilihFileCitraButton.Text = 'Pilih File Citra';

            % Create PathcitraLabel
            app.PathcitraLabel = uilabel(app.KonvolusiTab);
            app.PathcitraLabel.Position = [33 501 59 22];
            app.PathcitraLabel.Text = 'Path citra:';

            % Create PathInfoLabel
            app.PathInfoLabel = uilabel(app.KonvolusiTab);
            app.PathInfoLabel.Position = [111 502 776 22];
            app.PathInfoLabel.Text = '';

            % Create PilihanMaskUntukKonvolusiButtonGroup
            app.PilihanMaskUntukKonvolusiButtonGroup = uibuttongroup(app.KonvolusiTab);
            app.PilihanMaskUntukKonvolusiButtonGroup.Title = 'Pilihan Mask Untuk Konvolusi';
            app.PilihanMaskUntukKonvolusiButtonGroup.Position = [33 22 855 194];

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
            app.PreviewCitraMasukanPanel.Position = [33 231 274 256];

            % Create ImageAxes
            app.ImageAxes = uiaxes(app.PreviewCitraMasukanPanel);
            app.ImageAxes.XTick = [];
            app.ImageAxes.YTick = [];
            app.ImageAxes.Position = [19 1 234 229];

            % Create PreviewCitraHasilKonvolusiBuatanPanel
            app.PreviewCitraHasilKonvolusiBuatanPanel = uipanel(app.KonvolusiTab);
            app.PreviewCitraHasilKonvolusiBuatanPanel.Title = 'Preview Citra Hasil Konvolusi Buatan';
            app.PreviewCitraHasilKonvolusiBuatanPanel.Position = [325 231 274 256];

            % Create ImageAxes_2
            app.ImageAxes_2 = uiaxes(app.PreviewCitraHasilKonvolusiBuatanPanel);
            app.ImageAxes_2.XTick = [];
            app.ImageAxes_2.YTick = [];
            app.ImageAxes_2.Position = [20 1 234 229];

            % Create PreviewCitraHasilKonvolusiMATLABPanel
            app.PreviewCitraHasilKonvolusiMATLABPanel = uipanel(app.KonvolusiTab);
            app.PreviewCitraHasilKonvolusiMATLABPanel.Title = 'Preview Citra Hasil Konvolusi MATLAB';
            app.PreviewCitraHasilKonvolusiMATLABPanel.Position = [614 231 274 256];

            % Create ImageAxes_3
            app.ImageAxes_3 = uiaxes(app.PreviewCitraHasilKonvolusiMATLABPanel);
            app.ImageAxes_3.XTick = [];
            app.ImageAxes_3.YTick = [];
            app.ImageAxes_3.Position = [19 1 234 229];

            % Create JalankanKonvolusiButton
            app.JalankanKonvolusiButton = uibutton(app.KonvolusiTab, 'push');
            app.JalankanKonvolusiButton.ButtonPushedFcn = createCallbackFcn(app, @JalankanKonvolusiButtonPushed, true);
            app.JalankanKonvolusiButton.Enable = 'off';
            app.JalankanKonvolusiButton.Position = [769 531 119 23];
            app.JalankanKonvolusiButton.Text = 'Jalankan Konvolusi';

            % Create SmoothingLowPassFilterTab
            app.SmoothingLowPassFilterTab = uitab(app.TabGroup);
            app.SmoothingLowPassFilterTab.Title = 'Smoothing (Low Pass Filter)';

            % Create CitramasukanLabel_2
            app.CitramasukanLabel_2 = uilabel(app.SmoothingLowPassFilterTab);
            app.CitramasukanLabel_2.Position = [33 531 86 22];
            app.CitramasukanLabel_2.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_2
            app.PilihFileCitraButton_2 = uibutton(app.SmoothingLowPassFilterTab, 'push');
            app.PilihFileCitraButton_2.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_2Pushed, true);
            app.PilihFileCitraButton_2.Position = [133 531 100 23];
            app.PilihFileCitraButton_2.Text = 'Pilih File Citra';

            % Create PathcitraLabel_2
            app.PathcitraLabel_2 = uilabel(app.SmoothingLowPassFilterTab);
            app.PathcitraLabel_2.Position = [33 501 59 22];
            app.PathcitraLabel_2.Text = 'Path citra:';

            % Create PathInfoLabel_2
            app.PathInfoLabel_2 = uilabel(app.SmoothingLowPassFilterTab);
            app.PathInfoLabel_2.Position = [111 502 776 22];
            app.PathInfoLabel_2.Text = '';

            % Create PreviewCitraMasukanPanel_2
            app.PreviewCitraMasukanPanel_2 = uipanel(app.SmoothingLowPassFilterTab);
            app.PreviewCitraMasukanPanel_2.Title = 'Preview Citra Masukan';
            app.PreviewCitraMasukanPanel_2.Position = [111 166 344 324];

            % Create ImageAxes_4
            app.ImageAxes_4 = uiaxes(app.PreviewCitraMasukanPanel_2);
            app.ImageAxes_4.XTick = [];
            app.ImageAxes_4.YTick = [];
            app.ImageAxes_4.Position = [19 25 295 264];

            % Create PilihanFilterButtonGroup
            app.PilihanFilterButtonGroup = uibuttongroup(app.SmoothingLowPassFilterTab);
            app.PilihanFilterButtonGroup.Title = 'Pilihan Filter';
            app.PilihanFilterButtonGroup.Position = [33 15 834 142];

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
            app.nfilterorderforBLPFSpinner.Limits = [1 30];
            app.nfilterorderforBLPFSpinner.Position = [697 44 59 22];
            app.nfilterorderforBLPFSpinner.Value = 1;

            % Create PreviewCitraHasilFilterPanel
            app.PreviewCitraHasilFilterPanel = uipanel(app.SmoothingLowPassFilterTab);
            app.PreviewCitraHasilFilterPanel.Title = 'Preview Citra Hasil Filter';
            app.PreviewCitraHasilFilterPanel.Position = [471 166 344 324];

            % Create ImageAxes_5
            app.ImageAxes_5 = uiaxes(app.PreviewCitraHasilFilterPanel);
            app.ImageAxes_5.XTick = [];
            app.ImageAxes_5.YTick = [];
            app.ImageAxes_5.Position = [19 25 295 264];

            % Create JalankanSmoothingButton
            app.JalankanSmoothingButton = uibutton(app.SmoothingLowPassFilterTab, 'push');
            app.JalankanSmoothingButton.ButtonPushedFcn = createCallbackFcn(app, @JalankanSmoothingButtonPushed, true);
            app.JalankanSmoothingButton.Enable = 'off';
            app.JalankanSmoothingButton.Position = [743 531 125 23];
            app.JalankanSmoothingButton.Text = 'Jalankan Smoothing';

            % Create HighPassFilterTab
            app.HighPassFilterTab = uitab(app.TabGroup);
            app.HighPassFilterTab.Title = 'High Pass Filter';

            % Create CitramasukanLabel_3
            app.CitramasukanLabel_3 = uilabel(app.HighPassFilterTab);
            app.CitramasukanLabel_3.Position = [33 531 86 22];
            app.CitramasukanLabel_3.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_3
            app.PilihFileCitraButton_3 = uibutton(app.HighPassFilterTab, 'push');
            app.PilihFileCitraButton_3.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_3Pushed, true);
            app.PilihFileCitraButton_3.Position = [133 531 100 23];
            app.PilihFileCitraButton_3.Text = 'Pilih File Citra';

            % Create PathcitraLabel_3
            app.PathcitraLabel_3 = uilabel(app.HighPassFilterTab);
            app.PathcitraLabel_3.Position = [33 501 59 22];
            app.PathcitraLabel_3.Text = 'Path citra:';

            % Create PathInfoLabel_3
            app.PathInfoLabel_3 = uilabel(app.HighPassFilterTab);
            app.PathInfoLabel_3.Position = [111 502 776 22];
            app.PathInfoLabel_3.Text = '';

            % Create PreviewCitraMasukanPanel_3
            app.PreviewCitraMasukanPanel_3 = uipanel(app.HighPassFilterTab);
            app.PreviewCitraMasukanPanel_3.Title = 'Preview Citra Masukan';
            app.PreviewCitraMasukanPanel_3.Position = [111 166 344 324];

            % Create ImageAxes_6
            app.ImageAxes_6 = uiaxes(app.PreviewCitraMasukanPanel_3);
            app.ImageAxes_6.XTick = [];
            app.ImageAxes_6.YTick = [];
            app.ImageAxes_6.Position = [19 25 295 264];

            % Create PilihanFilterButtonGroup_2
            app.PilihanFilterButtonGroup_2 = uibuttongroup(app.HighPassFilterTab);
            app.PilihanFilterButtonGroup_2.Title = 'Pilihan Filter';
            app.PilihanFilterButtonGroup_2.Position = [234 15 459 142];

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
            app.nfilterorderforBHPFSpinner.Limits = [1 30];
            app.nfilterorderforBHPFSpinner.Position = [385 44 59 22];
            app.nfilterorderforBHPFSpinner.Value = 1;

            % Create PreviewCitraHasilFilterPanel_2
            app.PreviewCitraHasilFilterPanel_2 = uipanel(app.HighPassFilterTab);
            app.PreviewCitraHasilFilterPanel_2.Title = 'Preview Citra Hasil Filter';
            app.PreviewCitraHasilFilterPanel_2.Position = [471 166 344 324];

            % Create ImageAxes_7
            app.ImageAxes_7 = uiaxes(app.PreviewCitraHasilFilterPanel_2);
            app.ImageAxes_7.XTick = [];
            app.ImageAxes_7.YTick = [];
            app.ImageAxes_7.Position = [19 25 295 264];

            % Create JalankanHighPassFilterButton
            app.JalankanHighPassFilterButton = uibutton(app.HighPassFilterTab, 'push');
            app.JalankanHighPassFilterButton.ButtonPushedFcn = createCallbackFcn(app, @JalankanHighPassFilterButtonPushed, true);
            app.JalankanHighPassFilterButton.Enable = 'off';
            app.JalankanHighPassFilterButton.Position = [668 531 150 23];
            app.JalankanHighPassFilterButton.Text = 'Jalankan High Pass Filter';

            % Create FilterTerangHighBoostTab
            app.FilterTerangHighBoostTab = uitab(app.TabGroup);
            app.FilterTerangHighBoostTab.Title = 'Filter Terang (High Boost)';

            % Create CitramasukanLabel_4
            app.CitramasukanLabel_4 = uilabel(app.FilterTerangHighBoostTab);
            app.CitramasukanLabel_4.Position = [33 531 86 22];
            app.CitramasukanLabel_4.Text = 'Citra masukan:';

            % Create PilihFileCitraButton_4
            app.PilihFileCitraButton_4 = uibutton(app.FilterTerangHighBoostTab, 'push');
            app.PilihFileCitraButton_4.ButtonPushedFcn = createCallbackFcn(app, @PilihFileCitraButton_4Pushed, true);
            app.PilihFileCitraButton_4.Position = [133 531 100 23];
            app.PilihFileCitraButton_4.Text = 'Pilih File Citra';

            % Create PathcitraLabel_4
            app.PathcitraLabel_4 = uilabel(app.FilterTerangHighBoostTab);
            app.PathcitraLabel_4.Position = [33 501 59 22];
            app.PathcitraLabel_4.Text = 'Path citra:';

            % Create PreviewCitraMasukanPanel_4
            app.PreviewCitraMasukanPanel_4 = uipanel(app.FilterTerangHighBoostTab);
            app.PreviewCitraMasukanPanel_4.Title = 'Preview Citra Masukan';
            app.PreviewCitraMasukanPanel_4.Position = [111 166 344 324];

            % Create ImageAxes_8
            app.ImageAxes_8 = uiaxes(app.PreviewCitraMasukanPanel_4);
            app.ImageAxes_8.XTick = [];
            app.ImageAxes_8.YTick = [];
            app.ImageAxes_8.Position = [19 25 295 264];

            % Create PengaturanFilterButtonGroup
            app.PengaturanFilterButtonGroup = uibuttongroup(app.FilterTerangHighBoostTab);
            app.PengaturanFilterButtonGroup.Title = 'Pengaturan Filter';
            app.PengaturanFilterButtonGroup.Position = [111 15 699 142];

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
            app.nfilterorderforBLPFSpinner_2.Limits = [1 30];
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
            app.PreviewCitraHasilFilterPanel_3.Position = [471 166 344 324];

            % Create ImageAxes_9
            app.ImageAxes_9 = uiaxes(app.PreviewCitraHasilFilterPanel_3);
            app.ImageAxes_9.XTick = [];
            app.ImageAxes_9.YTick = [];
            app.ImageAxes_9.Position = [19 25 295 264];

            % Create JalankanFilterTerangButton
            app.JalankanFilterTerangButton = uibutton(app.FilterTerangHighBoostTab, 'push');
            app.JalankanFilterTerangButton.ButtonPushedFcn = createCallbackFcn(app, @JalankanFilterTerangButtonPushed, true);
            app.JalankanFilterTerangButton.Enable = 'off';
            app.JalankanFilterTerangButton.Position = [677 531 132 23];
            app.JalankanFilterTerangButton.Text = 'Jalankan Filter Terang';

            % Create PathInfoLabel_4
            app.PathInfoLabel_4 = uilabel(app.FilterTerangHighBoostTab);
            app.PathInfoLabel_4.Position = [111 502 776 22];
            app.PathInfoLabel_4.Text = '';

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
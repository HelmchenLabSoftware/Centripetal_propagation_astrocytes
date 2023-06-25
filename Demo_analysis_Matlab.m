
%% Compute delay maps from astrocytic in vivo calcium imaging data

% Code written by Peter Rupprecht (rupprecht@hifo.uzh.ch) in 2023

% Code is deposited in this repository: https://github.com/HelmchenLabSoftware/Centripetal_propagation_astrocytes
% See this repository for an explanation of the code

% Please cite this paper when using the code:
% Rupprecht, Peter, Christopher M. Lewis, and Fritjof Helmchen.
% "Centripetal integration of past events by hippocampal astrocytes."
% https://www.biorxiv.org/content/10.1101/2022.08.16.504030v1


% generic pattern of input files; if multiple files are recognized by this
% pattern, separate delay maps will be computed for each file, but the
% final result will be obtained by averaging across these delay maps
filenames = dir('Raw calcium imaging data/FOV_excerpt_recording*.tif');

% specific for this recording; change to make the scaling of the delay map correct
framerate = 4.419;



%% Go through each recording segment; delay maps will be computed for each segment and averaged afterwards
clear Corr_maps_all
for kkk = 1:numel(filenames)
    
    filename = fullfile(filenames(kkk).folder,filenames(kkk).name);
    disp(['Computing delay map for',32,filenames(kkk).name],'.');

    
    %% Read raw data from hard disk
    L = imfinfo(filename);
    TifLink = Tiff(filename, 'r');
    nb_frames = numel(L);
    movie = zeros(L(1).Height,L(1).Width,nb_frames,'uint16');
    for i = 1:nb_frames
        TifLink.setDirectory(i);
        movie(:,:,i) = TifLink.read();
    end 
    movie = double(movie);
    
    %% Compute reference (mean across FOV)
    mean_activity = squeeze(mean(mean(movie,1),2));

    %% Drift correction (align segments with respect to each other, if necessary)
    
    if kkk == 1
        
        % compute alignment template from first segment
        ref_excerpt = mean(movie,3);
        
    else
        
        % compute shift to alignment template using very simple movement correction
        this_excerpt = mean(movie,3);

        result_convC = fftshift(fftshift(real(ifft2(conj(fft2(ref_excerpt)).*fft2(ref_excerpt))), 1), 2);
        [x0,y0] = find(result_convC == max(result_convC(:)));

        result_convB = fftshift(fftshift(real(ifft2(conj(fft2(ref_excerpt)).*fft2(this_excerpt))), 1), 2);
        [x1,y1] = find(result_convB == max(result_convB(:)));

        dx = x1(1)-x0(1);
        dy = y1(1)-y0(1);

        movie = circshift(movie,-[dx dy 0]);
        
    end
    
    %% Compute delay map based on correlation function peaks
    max_delay = 30; % in #frames; do increase value if the frame rate is higher
    
    Corr_map = zeros(size(movie,1),size(movie,2));
    for j = 1:size(movie,2)
        % to speed up the program, replace the following line with this code: 
        % parfor k = 1:size(movie,1) (requires parallel computing toolbox)
        for k = 1:size(movie,1)
            % extract time trace of this pixel
            trace = squeeze( mean(mean(movie(k,j,:),1),2));
            % initialize cross_correlation vector
            cross_correlation = zeros(2*max_delay+1,1);
            for kk = -max_delay:max_delay
                % cross-correlate X and Y
                X = mean_activity;
                Y = circshift(trace,[kk 0]);
                if kk >= 0
                    cross_correlation(kk+31) = corr(X(kk+1:end),Y(kk+1:end));
                else
                    cross_correlation(kk+31) = corr(X(1:end+kk),Y(1:end+kk));
                end
            end
            % find peak of the cross-correlation
            [ix,xi] = max(cross_correlation);
            delay = xi - max_delay;
            % assign peak delay to the delay map pixel
            Corr_map(j,k) = delay;
        end
    end
    % save in one data structure
    Corr_maps_all{kkk} = Corr_map;
    
end
disp('All delay maps are computed.');

%% Average delay maps across segments
% Concatenate delay maps
Corr_maps_all_concatenated = [];
for kkk = 1:numel(Corr_maps_all)
    Corr_maps_all_concatenated = cat(3,Corr_maps_all_concatenated,Corr_maps_all{kkk});
end
% Average concatenated delay maps
delay_map = mean(Corr_maps_all_concatenated,3)/framerate;


%% Visualize results
figure(77), imagesc(delay_map); colormap(parula); caxis([-2 2]); colorbar; axis equal off
title('Map of delays (s)')

figure(78), imagesc(ref_excerpt); colormap(gray); caxis([quantile(ref_excerpt(:),0.0) quantile(ref_excerpt(:),0.95)]); axis equal off
title('Anatomical reference')



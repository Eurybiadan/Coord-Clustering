clear;
% close all;


SCALING_FACTOR = 10;

[pname1]=uigetdir(pwd, 'Select the first coordinate set. This folder needs to contain the reference images.');
[pname2]=uigetdir(pwd, 'Select the second coordinate set.');
[pname3]=uigetdir(pwd, 'Select the third coordinate set.');
[pname4]=uigetdir(pwd, 'Select the fourth coordinate set.');

[fNames] = read_folder_contents(pname1,'tif');
% delete('Summary.csv');
%% Load shit
for i=1:length(fNames)

    imName = fNames{i};
    coordName = [fNames{i}(1:end-4) '_coords.csv'];
    filepath1 = fullfile(pname1,coordName);
    filepath2 = fullfile(pname2,coordName);
    filepath3 = fullfile(pname3,coordName);
    filepath4 = fullfile(pname4,coordName);
        
    im = imread(fullfile(pname1,imName));
    coords{1} = dlmread(filepath1);
    coords{1} = coords{1}(round(coords{1}(:,1))>0 & round(coords{1}(:,1))<size(im,2), :);
    coords{1} = coords{1}(round(coords{1}(:,2))>0 & round(coords{1}(:,2))<size(im,1), :);
    
    coords{2} = dlmread(filepath2);
    coords{2} = coords{2}(round(coords{2}(:,1))>0 & round(coords{2}(:,1))<size(im,2), :);
    coords{2} = coords{2}(round(coords{2}(:,2))>0 & round(coords{2}(:,2))<size(im,1), :);
    
    coords{3} = dlmread(filepath3);
    coords{3} = coords{3}(round(coords{3}(:,1))>0 & round(coords{3}(:,1))<size(im,2), :);
    coords{3} = coords{3}(round(coords{3}(:,2))>0 & round(coords{3}(:,2))<size(im,1), :);
    
    coords{4} = dlmread(filepath4);
    coords{4} = coords{4}(round(coords{4}(:,1))>0 & round(coords{4}(:,1))<size(im,2), :);
    coords{4} = coords{4}(round(coords{4}(:,2))>0 & round(coords{4}(:,2))<size(im,1), :);
    
    allcoords = [coords{1}; coords{2}; coords{3}; coords{4}];
    coordlen = [length(coords{1}); length(coords{2}); length(coords{3}); length(coords{4})]';
    coordbounds = cumsum([1 coordlen]);
    samecoord = [];
    
    [smdists, distind] = pdist2(allcoords, allcoords,'euclidean','Smallest',length(coords));
    
    threshdist = smdists(2:end,:);
    threshold = mean(threshdist(threshdist~=0))*1; %median(threshdist(:))*1; %+std(threshdist(:));

    thresholddisk = strel('disk',round(threshold*SCALING_FACTOR)-1,0);
    
    simple_constellation = zeros(size(im)*SCALING_FACTOR,'uint8');
    simple_constellation = repmat(simple_constellation,[1 1 length(coords)]);
    
    labelled_constellation = zeros(size(im)*SCALING_FACTOR,'double');
    labelled_constellation = repmat(labelled_constellation,[1 1 length(coords)]); 
    
    for c=1:length(coords)
        
        RnS_coords{c} = round(coords{c}*SCALING_FACTOR);
        
        inds = sub2ind(size(im)*SCALING_FACTOR,RnS_coords{c}(:,1),RnS_coords{c}(:,2));
        
        this_constellation = false(size(im)*10);
        this_constellation( inds ) = true;
        this_constellation = imdilate(this_constellation, thresholddisk);
        
        simple_constellation(:,:,c) = uint8(this_constellation);
        
        this_constellation = zeros(size(im)*10);
        for j=1:size(RnS_coords{c},1)
            this_constellation( RnS_coords{c}(j,2), RnS_coords{c}(j,1) ) = j;
        end
        this_constellation = conv2(this_constellation, thresholddisk.Neighborhood,'same');
        labelled_constellation(:,:,c) = double(this_constellation');
    end
    
    flattened_constellations = uint8(sum(simple_constellation,3));
    
    figure(1); clf; imagesc( flattened_constellations); axis image; hold on;
    theconncomps = bwconncomp(flattened_constellations);
    
    
    clustercoords=[];
    
    %Go through our connected components and separate ones that have two
    % "peak" value areas; these (concievably) are not real.
    for c=1:length(theconncomps.PixelIdxList)
    
        mask = zeros(size(im)*SCALING_FACTOR,'uint8');
        
        mask(theconncomps.PixelIdxList{c}) = 1;
        masked_constellation = flattened_constellations.*mask;
        maximus = imregionalmax(masked_constellation);
        figure(2); clf; imagesc(masked_constellation); axis image;
        subconncomp = bwconncomp(maximus);
        
        if subconncomp.NumObjects > 1
                        
            subregionstats = regionprops(subconncomp, maximus,'MinorAxisLength','Area', 'Centroid' );
            
            for o=1:subconncomp.NumObjects
                if subregionstats(o).Area > 5
                    clustercoords = [clustercoords; subregionstats(o).Centroid]; 
                end
            end
            
            figure(3); clf; imagesc(flattened_constellations.*uint8(maximus)); 
            axis image;
        else
            subregionstats = regionprops(subconncomp, 'Centroid' );
            
            clustercoords=[clustercoords; subregionstats.Centroid];
            
        end
    end
    
    clustercoords=round(clustercoords);
    
    clusterlist = nan( size(clustercoords,1), length(coords) );
    
    for c=1:size(clustercoords,1)
        clusterlist(c,:)= labelled_constellation(clustercoords(c,2), clustercoords(c,1),:);
    end
    %% Reconcile multiple clusters overlapping
    for m=1:length(coords)
        for c=1:max(clusterlist(:,m))
            
            if size(clusterlist(clusterlist(:,m)==c,:),1) > 1

                overlapped_inds = find(clusterlist(:,m)==c);
                overlapped_clusters = clusterlist(overlapped_inds,:);
                overlapped_clusters(overlapped_clusters==0) = NaN;
                
                uniqueinds = unique(overlapped_clusters(overlapped_clusters~=0));
                
                if length(uniqueinds)>length(coords)

                    allcombos = fliplr(combvec_matrix(1:size(overlapped_clusters,1), length(coords))');
                    
                    overlapamt = zeros(size(allcombos,1),1);
                    for j=1:size(allcombos,1)
                        
                        overlapmask = zeros( size(labelled_constellation) );                        
                        for k=1:size(allcombos,2)

                            overlap_cluster = clusterlist(overlapped_clusters(allcombos(j,k),k), m);
                            overlapmask = overlapmask | (labelled_constellation(:,:,k)== overlap_cluster );
                        end
                        
                         imagesc(sum(overlapmask,3)); axis image; 
                         pause(0.1)
                        overlapamt(j)=sum(overlapmask(:));
                    end

                    [overlapmax,overlapwinner] = max(overlapamt);

                    winnerinds = sub2ind(size(overlapped_clusters),allcombos(overlapwinner,:), 1:length(coords));
                    overlapped_clusters( winnerinds );
                    clusterlist(overlapped_inds(1),:) = overlapped_clusters( winnerinds );
                    
%                     overlapped_clusters( winnerinds ) = [];
                    clusterlist(overlapped_inds(2:end),:)=0;                    
                
                else
                    clusterlist(overlapped_inds(1), 1:length(uniqueinds)) = uniqueinds;
                    
                    clusterlist(overlapped_inds(2:end),:)=0;
                end
                                                
%                 pause;
            end

        end
    end
    %% Display clusters
    figure(10); clf; hold on; imagesc( flattened_constellations' ); colormap gray; axis image;
    overallmask = zeros( size(flattened_constellations));
    for c=1:size(clusterlist,1)
        if all(clusterlist(c,:)==0)
            continue;
        else
            showableinds = clusterlist(c,:)~=0;
        end
        
        overlapmask = zeros( [size(labelled_constellation,1),size(labelled_constellation,2)] );
        
        for k=1:length(coords)
            if clusterlist(c,k)~=0
                overlapmask = overlapmask + (labelled_constellation(:,:,k)==clusterlist(c,k));
                overallmask = overallmask | (labelled_constellation(:,:,k)==clusterlist(c,k));
            end
        end                
        overlapprops = regionprops(bwconncomp(overlapmask),'Centroid');
        
        clustercenter(c,:) = [overlapprops.Centroid(2) overlapprops.Centroid(1)]; %overlapprops.Centroid;

        % % if max(overlapmask(:)) == 4
        % %     plot(thiscenter(:,1),thiscenter(:,2),'g*');
        % if sum(showableinds) >= 3
        %     plot(clustercenter(c,1),clustercenter(c,2),'g*');
        % elseif sum(showableinds) == 2
        % 
        %     if all(clusterlist(c, 2:3)>0) % If JM and RC marked it but GH didn't, mark a GH colored X, and so on.
        %         plot(clustercenter(c,1),clustercenter(c,2),'mx');
        %     elseif all(clusterlist(c, [1 3])>0)  % GH and RC did, JM didn't
        %         plot(clustercenter(c,1),clustercenter(c,2),'yx');
        %     elseif all(clusterlist(c, 1:2)>0)
        %         plot(clustercenter(c,1),clustercenter(c,2),'bx'); % GH and JM did, RC didn't
        %     else                
        %         disp('You should not see this!')
        %     end
        % elseif sum(showableinds) == 1
        %     [~, whichobs]=max(clusterlist(c,:));
        %     switch(whichobs)
        %         case 1
        %             color='m.';
        %         case 2
        %             color='y.';
        %         case 3
        %             color='b.';
        %     end
        %     plot(clustercenter(c,1),clustercenter(c,2),color,'MarkerSize',12)
        % 
        % end
        if sum(showableinds) == 4
             plot(clustercenter(c,1),clustercenter(c,2),'b*');
        elseif sum(showableinds) == 3
            plot(clustercenter(c,1),clustercenter(c,2),'g*');
        elseif sum(showableinds) == 2            
            plot(clustercenter(c,1),clustercenter(c,2),'r*');
        elseif sum(showableinds) == 1
            plot(clustercenter(c,1),clustercenter(c,2),'y*');            
        end

    end
    %%
    title('Constellation');


    saveas(gcf, fullfile(pname1, [fNames{i}(1:end-4) '_overlap.png']));
%     saveas(gcf, fullfile(pname, [fNames{i}(1:end-4) '_agreement.svg']))
    
    % Make our grid of true/false positive and dice coefficient
           
    % theperms = perms(1:3);
    % theperms = theperms(:,1:2);
    % fid=fopen(fullfile(pname1,'Summary.csv'),'a');
    % 
    % for p=1:size(theperms,1)
    % 
    %     goldstdexist = existlist(:,theperms(p,1));
    %     comparisoneixst = existlist(:,theperms(p,2));
    % 
    %     numtruepositive = sum( goldstdexist & comparisoneixst );
    %     numfalsepositive = sum( goldstdexist < comparisoneixst );
    % 
    %     numfalsenegative = sum( goldstdexist > comparisoneixst );
    % 
    %     numauto = numtruepositive+numfalsepositive;
    %     nummanual = numtruepositive+numfalsenegative;
    % 
    %     truepositiverate = numtruepositive/nummanual;
    %     falsediscoveryrate = numfalsepositive/numauto;
    %     dicecoeff = 2*numtruepositive ./ (nummanual+numauto);
    %     tmp = goldstdexist & comparisoneixst;
    % 
    %     fprintf(fid, '%s,%d,%d,%f,%f,%f,\n',imName, theperms(p,1), theperms(p,2), truepositiverate, falsediscoveryrate, dicecoeff);
    % end
    % fclose(fid);
    
end
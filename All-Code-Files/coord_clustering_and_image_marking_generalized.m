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
    coord_lists{1} = dlmread(filepath1);
    coord_lists{1} = coord_lists{1}(round(coord_lists{1}(:,1))>0 & round(coord_lists{1}(:,1))<size(im,2), :);
    coord_lists{1} = coord_lists{1}(round(coord_lists{1}(:,2))>0 & round(coord_lists{1}(:,2))<size(im,1), :);
    
    coord_lists{2} = dlmread(filepath2);
    coord_lists{2} = coord_lists{2}(round(coord_lists{2}(:,1))>0 & round(coord_lists{2}(:,1))<size(im,2), :);
    coord_lists{2} = coord_lists{2}(round(coord_lists{2}(:,2))>0 & round(coord_lists{2}(:,2))<size(im,1), :);
    
    coord_lists{3} = dlmread(filepath3);
    coord_lists{3} = coord_lists{3}(round(coord_lists{3}(:,1))>0 & round(coord_lists{3}(:,1))<size(im,2), :);
    coord_lists{3} = coord_lists{3}(round(coord_lists{3}(:,2))>0 & round(coord_lists{3}(:,2))<size(im,1), :);
    
    coord_lists{4} = dlmread(filepath4);
    coord_lists{4} = coord_lists{4}(round(coord_lists{4}(:,1))>0 & round(coord_lists{4}(:,1))<size(im,2), :);
    coord_lists{4} = coord_lists{4}(round(coord_lists{4}(:,2))>0 & round(coord_lists{4}(:,2))<size(im,1), :);
    
    allcoords = [coord_lists{1}; coord_lists{2}; coord_lists{3}; coord_lists{4}];
    coordlen = [length(coord_lists{1}); length(coord_lists{2}); length(coord_lists{3}); length(coord_lists{4})]';
    coordbounds = cumsum([1 coordlen]);
    samecoord = [];
    
    [smdists, distind] = pdist2(allcoords, allcoords,'euclidean','Smallest',length(coord_lists));
    
    threshdist = smdists(2:end,:);
    threshold = mean(threshdist(threshdist~=0))*1; %median(threshdist(:))*1; %+std(threshdist(:));

    thresholddisk = strel('disk',round(threshold*SCALING_FACTOR)-1,0);
    
    simple_constellation = zeros(size(im)*SCALING_FACTOR,'uint8');
    simple_constellation = repmat(simple_constellation,[1 1 length(coord_lists)]);
    
    labelled_constellation = zeros(size(im)*SCALING_FACTOR,'double');
    labelled_constellation = repmat(labelled_constellation,[1 1 length(coord_lists)]); 
    
    for c=1:length(coord_lists)
        
        RnS_coords{c} = round(coord_lists{c}*SCALING_FACTOR);
        
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
        figure(2); clf; imagesc(masked_constellation); axis image; hold on;
        subconncomp = bwconncomp(maximus);
        
        if subconncomp.NumObjects > 1
                        
            subregionstats = regionprops(subconncomp, maximus,'MinorAxisLength','Area', 'Centroid' );
            
            for o=1:subconncomp.NumObjects
                if subregionstats(o).Area > 5
                    clustercoords = [clustercoords; subregionstats(o).Centroid]; 
                    plot( subregionstats(o).Centroid(1), subregionstats(o).Centroid(2),'*')
                end
            end
            
            figure(3); clf; imagesc(flattened_constellations.*uint8(maximus)); 
            axis image;
        else
            subregionstats = regionprops(subconncomp, 'Centroid' );
            
            clustercoords=[clustercoords; subregionstats.Centroid];
            plot( subregionstats.Centroid(1), subregionstats.Centroid(2),'*')
        end
        
    end
    
    clustercoords=round(clustercoords);
    
    clusterlist = nan( size(clustercoords,1), length(coord_lists) );
    
    for c=1:size(clustercoords,1)
        clusterlist(c,:)= labelled_constellation(clustercoords(c,2), clustercoords(c,1),:);
    end
    %% Reconcile multiple clusters overlapping
    reconciled_clusterlist = [];

% for m=1:length(coord_lists)
        for c=1:size(clusterlist,1)
            
            coi = clusterlist(c,:);

            % Find all rows containing each cluster member.
            roi = [];
            for m=1:length(coord_lists)
                if coi(m) ~= 0
                    roi = [roi; clusterlist(clusterlist(:,m) == coi(m),:)];
                end
            end

            roi = unique(roi, 'rows');

            % Find all instances of the members of cluster c in the coordinate lists-
            % if it exists more than once, then go in here- otherwise, do
            % nothing
            if size(roi,1) > 1

                % Get the rows where this cluster appears more than once
                % and the clusters themselves- if there's a 0, that means
                % there is no cluster- so we could have clusters with m-1,
                % m-2, etc. members.
                overlapped_clusters = roi;
                overlapped_clusters(overlapped_clusters==0) = NaN;
                                
                while size(overlapped_clusters,1) > 0
                                              
                    % Display of each cluster. Uncomment if curious.
                    overlapmask = zeros( size(labelled_constellation,1), size(labelled_constellation,2) );  
                    for k=1:size(overlapped_clusters,1)
                        for cl=1:length(coord_lists)
                            what_a_cluster = labelled_constellation(:,:,cl)==overlapped_clusters(k,cl);
                            overlapmask = overlapmask + what_a_cluster;
                            
                            figure(42); clf; hold on;
                            imagesc(overlapmask); axis image; 
                            % disp(num2str(overlapped_clusters(k,cl)))
                        end
                    end

                    allcombos = fliplr(combvec_matrix(1:size(overlapped_clusters,1), length(coord_lists))');
                    
                    overlapamt = zeros(size(allcombos,1),1);
                    for j=1:size(allcombos,1)
                                                
                        overlapmask = zeros( size(labelled_constellation,1), size(labelled_constellation,2) );  
                        for k=1:size(allcombos,2)

                            what_a_cluster = overlapped_clusters(allcombos(j,k), k);                            
                            overlapmask = overlapmask | (labelled_constellation(:,:,k)== what_a_cluster );
                        end
                        
                        % figure(42); hold on;
                        %  imagesc(sum(overlapmask,3)); axis image; 
                        %  pause(0.1)
                        overlapamt(j)=sum(overlapmask(:));
                    end
                    
                    % Find out which has the greatest overlap- take the
                    % first one that has the greatest overlap.
                    [overlapmax,overlapwinner] = max(overlapamt);

                    winnerinds = sub2ind(size(overlapped_clusters),allcombos(overlapwinner,:), 1:length(coord_lists));
                    winning_cluster = overlapped_clusters( winnerinds );
                    
                    disp('Found greatest overlap. Adding best option, cleaning up cluster list...')

                    % Add them to our reconciled clusterlist
                    reconciled_clusterlist = [reconciled_clusterlist; winning_cluster];       

                    % Remove the winning indexes from the overlapped cluster
                    for j=1:size(overlapped_clusters,1)
                        for k=1:size(overlapped_clusters,2)
                            if overlapped_clusters(j,k) == winning_cluster(k) 
                                overlapped_clusters(j,k) = NaN;
                            end
                        end
                    end

                    % Remove blank rows.
                    overlapped_clusters = overlapped_clusters(~all(isnan(overlapped_clusters),2),:);

                    % If there's only one row, then just add it to our
                    % recoiled list- if there's multiple, but unique rows,
                    % then do the same.
                    if size(overlapped_clusters,1) == 1
                        disp('Adding remaining cluster to list.')
                        reconciled_clusterlist = [reconciled_clusterlist; overlapped_clusters];
                        overlapped_clusters = [];
                    elseif size(overlapped_clusters,1) > 1

                        % Do a check for columnwise uniqueness.
                        we_are_special_snowflakes = true;
                        for k=1:size(overlapped_clusters,2)
                            if length(unique(overlapped_clusters(:,k))) == 1
                                we_are_special_snowflakes=false;
                            end
                        end
                        % If each column is unique, then each remaining row
                        % of the cluster is unique. We then need to check
                        % to see if these remaining clusters are actually
                        % clusters (e.g. they're actually connected). If
                        % they aren't, separate them and add them to the
                        % clusterlist. If they are, just add them.
                        if we_are_special_snowflakes
                            disp('Adding remaining clusters to list.')

                            for j=1:size(overlapped_clusters,1)
                                thiscluster = overlapped_clusters(j,:);
                                
                                if sum(~isnan(thiscluster)) == 1
                                    reconciled_clusterlist = [reconciled_clusterlist; thiscluster];
                                else
                                    overlapmask = zeros( size(labelled_constellation,1), size(labelled_constellation,2) );  
                                    for k=1:size(allcombos,2)            
                                        overlapmask = overlapmask + (labelled_constellation(:,:,k)== thiscluster(k) );
                                    end
                                    overlapconcomp = bwconncomp(overlapmask>0);
                                    overlapprops = regionprops(bwconncomp(overlapmask>0),'Centroid');
                                    
                                    if overlapconcomp.NumObjects > 1
                                        %If we have more than one object, then
                                        %add each new cluster as we did above.                                 
                                        for o=1:overlapconcomp.NumObjects
                                            if subregionstats(o).Area > 5
                                                subreg_coords = subregionstats(o).Centroid; 
                                                
                                                thiscluster = labelled_constellation(subreg_coords (c,2), subreg_coords (c,1),:);
                                    
                                                thiscluster(thiscluster==0)=NaN;
                                                reconciled_clusterlist  = [reconciled_clusterlist ; thiscluster];
                                            end
                                        end
                                       
                                    else
                                        reconciled_clusterlist = [reconciled_clusterlist; thiscluster];
                                    end
                                end
                            end
                            overlapped_clusters = [];                            
                            
                        else
                            disp('We are not special snowflakes.')                                
                        end
                    end
                    
                end
                disp(['Done reconciling index: ' num2str(c)])

            else
                thiscluster = roi;
                thiscluster(thiscluster==0)=NaN;                
                reconciled_clusterlist = [reconciled_clusterlist; thiscluster];
            end

        end
    % end
    
    %% Display clusters
    figure(10); clf; hold on; imagesc( flattened_constellations' ); colormap gray; axis image;
    overallmask = zeros( size(flattened_constellations));

   
    for c=1:size(reconciled_clusterlist,1)
        
        showableinds = ~isnan(reconciled_clusterlist(c,:));
        
        overlapmask = zeros( size(labelled_constellation,1), size(labelled_constellation,2) );  

        for k=1:length(coord_lists)          
                overlapmask = overlapmask + (labelled_constellation(:,:,k)==reconciled_clusterlist(c,k));
                overallmask = overallmask + (labelled_constellation(:,:,k)==reconciled_clusterlist(c,k));            
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
            plot(clustercenter(c,1),clustercenter(c,2),'y*');
        elseif sum(showableinds) == 1
            plot(clustercenter(c,1),clustercenter(c,2),'r*');            
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
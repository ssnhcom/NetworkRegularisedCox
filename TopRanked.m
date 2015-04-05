%% Get Top Ranked Genes List
function[TopRankedGenes] = TopRanked(genes_list)

M = size(genes_list,1); % Number of genes
N = size(genes_list,2); % Number of evaluation

clear prob_list, clear tmp_list
TopRankedGenes = zeros(M,N);
for i = 1:M
    disp(i)
    if(i==1)
        tmp = genes_list(i,:);
    else
        tmp = [tmp, genes_list(i,:)];
    end
    score = tabulate(tmp);
    for j = 1:N
        %fprintf('j %d\n', j);
        for v = 1:i
            %fprintf('v %d\n', v);
            indexScore = strfind(score(:,1), char(genes_list(v,j)));
            index = find(not(cellfun('isempty', indexScore)));
            if(size(index,1) > 1)
                for k = 1:size(index,1)
                    %fprintf('k %d\n', k);
                    if(strcmp(genes_list(v,j),score(index(k,1),1)))
                        index = index(k,1);
                        break;
                    end
                end
            end
            %fprintf('here\n');
            tmp_list(v,j) = cell2mat(score(index,3));
        end
    end
    if( i == 1)
        TopRankedGenes(i,:) = tmp_list;
    else
        TopRankedGenes(i,:) = sum(tmp_list);
    end
end


clear tmp, clear tmp_list, clear i, clear j, clear v, clear k, clear M, clear N, clear indexScore, clear index, clear score
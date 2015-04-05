% nmi: Normalized Mutual Information
% mi; Mutual Information
function mi = MutualInformation( A, B )

total = length(A);
A_ids = unique(A);
A_class = length(A_ids);
B_ids = unique(B);
B_class = length(B_ids);

% Mutual information
idAOccur = double (repmat( A, A_class, 1) == repmat( A_ids', 1, total ));
idBOccur = double (repmat( B, B_class, 1) == repmat( B_ids', 1, total ));
idABOccur = idAOccur * idBOccur';

Px = sum(idAOccur') / total;
Py = sum(idBOccur') / total;
Pxy = idABOccur / total;

MImatrix = Pxy .* log2(Pxy ./(Px' * Py)+eps);
mi = sum(MImatrix(:));

% Entropies
Hx = -sum(Px .* log2(Px + eps),2);
Hy = -sum(Py .* log2(Py + eps),2);


%Normalized Mutual information
% http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
nmi = 2 * mi / (Hx+Hy);

% nmi = MI / sqrt(Hx*Hy); another version of NMI
%http://i11www.iti.uni-karlsruhe.de/extra/publications/ww-cco-06.pdf

end
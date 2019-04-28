ex = csvread('gene_interaction\expression.csv');
ex = transpose(ex);
ex = abs(ex);
norm_ex = ex/max(ex);

listFS = {'co_expression','DO', 'GObp', 'GOcc', 'GOmf', 'HPO','PPI','pubmed'};
listFS = {'joint_interaction'};

for m = 1:length(listFS)
    selection_method = listFS{m}
    disp(selection_method)
    din = strcat('gene_interaction\', selection_method, '.csv')
    W = csvread(din);
    w = sparse(W);
    degrees = sum(W);
    ind = find(degrees == 0);
    degrees(ind) = 1;
    D1 = sparse(diag(1./degrees));
    d = 0 
    while d < 1.05
        ranking = degrees / norm(degrees,1);
        if d<1.0
            A = eye(size(w)) - d*(w'*D1);
            b = (1-d)*norm_ex;
            ranking = A\b;
        end                 
        dout = strcat('gene_interaction\score\', selection_method, '_', num2str(d*100), '.txt')
        fileID = fopen(dout,'w');
        fprintf(fileID,'%f\n',ranking);
        fclose(fileID);
        d = d + .05
    end
end
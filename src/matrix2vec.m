function match_pair = matrix2vec(M)
Map = zeros(2,length(M(:,1)));
for i = 1:length(M(:,1))
    indi = find(M(i,:)==1);
    if isempty(indi)~=1
        Map(1,i) = i;
        Map(2,i) = indi(1);
    end
end
zero_ind = ismember(Map',[0,0],'rows');
match_pair = Map(:,zero_ind==0);
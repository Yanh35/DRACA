function [newF] = DRACA_Demo(nTypes,instanseIdx,Gcell,Rcell)

LD=Rcell{3};
[NL,ND]=size(LD);
max_iter = 20;
alpha = 100000;
G_enum = cell(size(Gcell));
G_denom = cell(size(Gcell));

threshold = 0.00001;


for ii = 1:length(Gcell)
    G_enum{ii}=0;
    G_denom{ii} =0;
end
Scell = cell(size(Gcell));
for iter = 1:max_iter
    %% initialize the S

    mus=zeros(length(Rcell),1);


    for rr = 1:length(instanseIdx)
        i = fix(instanseIdx{rr}/nTypes);
        j = mod(instanseIdx{rr},nTypes)+1;
        if j ==0
                i = i-1;
                j = 6;
        end
        Gmatii = Gcell{i};
        Gmatjj = Gcell{j};

        Rmat = Rcell{rr};

        Smat = (Gmatii'*Gmatii)\Gmatii'*Rmat*Gmatjj/(Gmatjj'*Gmatjj);
        Smat(isnan(Smat))=0;
        Scell{rr}=Smat;

        result = Rcell{rr}-Gmatii*Scell{rr}*Gmatjj';
        R = sum(sum(result.^2));
        mus(rr)=R;
    end

    %% get the optimal weights for relation matrices according to the reconstruction loss
    [bestP,Ws]=getOptimalWeights(mus,alpha);

    %% update G with relation matrices and constraint matrices
    for rr = 1:length(instanseIdx)
        i = fix(instanseIdx{rr}/nTypes);
        j = mod(instanseIdx{rr},nTypes)+1;
        if j ==0
                i = i-1;
                j = 6;
        end
        temp1 = Rcell{rr}*Gcell{j}*Scell{rr}';
        temp1(isnan(temp1))=0;
        t = abs(temp1);

        temp1p = (t+temp1)/2;
        temp1n = (t-temp1)/2;

        temp2 = Scell{rr}*Gcell{j}'*Gcell{j}*Scell{rr}';
        temp2(isnan(temp2))=0;
        t = abs(temp2);

        temp2p = (t+temp2)/2;
        temp2n = (t-temp2)/2;

        temp3 = Rcell{rr}'*Gcell{i}*Scell{rr};
        temp3(isnan(temp3))=0;
        t = abs(temp3);
        t(t<0)=0;
        temp3p = (t+temp3)/2;
        temp3n = (t-temp3)/2;

        temp4 = Scell{rr}'*Gcell{i}'*Gcell{i}*Scell{rr};
        temp4(isnan(temp4))=0;
        t = abs(temp4);
        t(t<0)=0;
        temp4p = (t+temp4)/2;
        temp4n = (t-temp4)/2;
        G_enum{i} = G_enum{i}+Ws(rr).* temp1p+Ws(rr).*Gcell{i}*temp2n;
        G_denom{i}= G_denom{i}+Ws(rr).*temp1n+Ws(rr).*Gcell{i}*temp2p;

        G_enum{j} = G_enum{j}+ Ws(rr).*temp3p+Ws(rr).*Gcell{j}*temp4n;
        G_denom{j}= G_denom{j}+Ws(rr).*temp3n+Ws(rr).*Gcell{j}*temp4p;

    end

    for ii = 1:length(Gcell)
        G_denom{ii}=G_denom{ii}+eps;
        factor = sqrt(G_enum{ii}./G_denom{ii});
        Gcell{ii}=Gcell{ii}.*factor;
        Gcell{ii}(isnan(Gcell{ii}))=0;
        Gcell{ii}(isinf(Gcell{ii}))=0;
    end

    %% compare the target approximation (||R15-G1S15G5'||^2) with threshold 

    result = Rcell{3}-Gcell{1}*Scell{3}*Gcell{4}';
    R = sum(sum(result.^2));
    if R<threshold
        break;
    end
end

    
newF=Gcell{1}*Scell{3}*Gcell{4}';

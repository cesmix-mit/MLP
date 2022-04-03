
elem = ["Ni", "Cu", "Li", "Mo", "Si", "Ge"];

for i = 1:length(elem)
    cd(elem(i));
    filename = "trainerror.txt";
    fileID = fopen(filename,'r');
    err = fscanf(fileID, '%f');
    fclose(fileID);
    err = reshape(err,[6, 6])';
    err(:,1:3) = 1e3*err(:,1:3);
    cd("..");
    if i == 1
        eNi = err;
    elseif i== 2
        eCu = err;
    elseif i== 3
        eLi = err;        
    elseif i== 4
        eMo = err;
    elseif i== 5
        eSi = err;        
    elseif i== 6
        eGe = err;        
    end
end

eACE=[eNi(:,1) eCu(:,1) eLi(:,1) eMo(:,1) eSi(:,1) eGe(:,1)];
i=2;eSNAP=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];
i=3;ePOD=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];
i=4;fACE=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];
i=5;fSNAP=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];
i=6;fPOD=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];

for i = 1:length(elem)
    cd(elem(i));
    filename = "testerror.txt";
    fileID = fopen(filename,'r');
    err = fscanf(fileID, '%f');
    fclose(fileID);
    err = reshape(err,[6, 6])';
    err(:,1:3) = 1e3*err(:,1:3);
    cd("..");
    if i == 1
        eNi = err;
    elseif i== 2
        eCu = err;
    elseif i== 3
        eLi = err;        
    elseif i== 4
        eMo = err;
    elseif i== 5
        eSi = err;        
    elseif i== 6
        eGe = err;        
    end
end

eACEtest=[eNi(:,1) eCu(:,1) eLi(:,1) eMo(:,1) eSi(:,1) eGe(:,1)];
i=2;eSNAPtest=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];
i=3;ePODtest=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];
i=4;fACEtest=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];
i=5;fSNAPtest=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];
i=6;fPODtest=[eNi(:,i) eCu(:,i) eLi(:,i) eMo(:,i) eSi(:,i) eGe(:,i)];

% k = 4;
% [eACE(k,:); fACE(k,:); eSNAP(k,:); fSNAP(k,:); ePOD(k,:); fPOD(k,:)]

for k = 2:2:6
[eACE(k,:); eACEtest(k,:); eSNAP(k,:); eSNAPtest(k,:); ePOD(k,:); ePODtest(k,:)]
end


for k = 2:2:6
[fACE(k,:); fACEtest(k,:); fSNAP(k,:); fSNAPtest(k,:); fPOD(k,:); fPODtest(k,:)]
end




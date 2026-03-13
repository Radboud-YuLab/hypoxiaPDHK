function [subModel, rxnsIdxClosed] = getSubModel(model, subsysToKeep, rxnsToKeep, subsysToClose, rxnsToClose, ...
                                                 keepLipids, keepATP, keepATPox, keepGln)
% reduce a model to a smaller submodel. 
%
%   model           An ihuman-based GEM
%   subsysToKeep    (cell array) subsystems to keep
%   rxnsToKeep      (cell array) reactions to keep
%   subsysToClose   (cell array) subsystems to close
%   rxnsToClose     (cell array) reactions to close
%   keepLipids      flag whether to force non-zero flux through lipids
%   keepATP         flag whether to force non-zero flux through ATP (NGAM)
%   keepATPox       flag whether to force non-zero flux through ATP (oxphos)
%   keepGln         flag whether to force non-zero flux through Gln-Glu-AKG
%   
%   subModel        smaller model with only indicated subsystems
%   rxnsIdxClosed   reactions that have been removed, given as indices of the OG model
%
% Rosemary Yu. 		Last edited: 2026-03-13


% set bounds
model = setHams(model);

mbiomass = find(ismember(model.rxns, 'MAR09932')); %maintenance biomass
ATP = find(ismember(model.rxns, 'MAR03964')); %NGAM
oATP = find(ismember(model.rxns, 'MAR06916')); %oxphos-produced ATP
glc = find(ismember(model.rxns, 'MAR09034')); % glucose exchange

model.ub(mbiomass) = 1000;
model.lb(glc) = -1;
model.ub(ATP) = 1000;
model.ub(oATP) = 1000;

% close rxns in subsysToClose and rxnsToClose
rxns = [];
for idx = 1:numel(model.subSystems)
    i = model.subSystems{idx};
    if any( strcmp(subsysToClose, string(i)) ) 
        rxns{end+1,1} = model.rxns{idx}; 
    end
end
rxns = [rxns; rxnsToClose];

rxnsIdxClosed = getIndexes(model, rxns, 'rxns');
model.ub(rxnsIdxClosed) = 0;
model.lb(rxnsIdxClosed) = 0;

% check model has flux through lipids, ATP, oATP, and Gln
[f1,f2,f3,f4] = checkFlag(model, keepLipids, keepATP, keepATPox, keepGln);
if f1 == 0 || f2 == 0 || f3 == 0 || f4 == 0
    disp(' Too many subsystems removed, please check "subsysToRemove" ')
else
    % find reactions to keep 
    rxnsToKeepIdx = getIndexes(model, rxnsToKeep, 'rxns');

    %find reactions to close
    rxns = [];
    for idx = 1:numel(model.subSystems)
        i = model.subSystems{idx};
        if ~any( strcmp(subsysToKeep, string(i)) ) && ~ismember(idx, rxnsToKeepIdx) && ~ismember(idx, rxnsIdxClosed)
            rxns{end+1,1} = model.rxns{idx}; 
        end
    end
    rxnIdx = getIndexes(model, rxns, 'rxns');

    % close these reactions one at a time    
    for i = 1:numel(rxnIdx) 
        if rem(i,20) == 0
            disp(['closing ' num2str(i) ' of ' num2str(numel(rxnIdx)) ' rxns'])
        end
        % save current ub and lb
        ub = model.ub(rxnIdx(i));
        lb = model.lb(rxnIdx(i));
        % close the rxn
        model.ub(rxnIdx(i)) = 0;
        model.lb(rxnIdx(i)) = 0;
        % check model has flux through lipids, ATP, oATP, and Gln
        [f1,f2,f3,f4] = checkFlag(model, keepLipids, keepATP, keepATPox, keepGln);
        if f1 == 0 || f2 == 0 || f3 == 0 || f4 == 0
            model.ub(rxnIdx(i)) = ub;
            model.lb(rxnIdx(i)) = lb;
        else
            %if there is solution then can be closed
            rxnsIdxClosed(end+1) = rxnIdx(i); 
        end
    end

    model = removeReactionsFull(model, rxnsIdxClosed, true, true, true);
    subModel = removeDeadEnds(model);
end

end %end function getSubModel




% local
function [res_lipids, res_ATP, res_oATP, res_Gln] = checkFlag(model, keepLipids, keepATP, keepATPox, keepGln)

%mbiomass = find(ismember(model.rxns, 'MAR09932')); %maintenance biomass 
ATP = find(ismember(model.rxns, 'MAR03964')); %NGAM
oATP = find(ismember(model.rxns, 'MAR06916')); %oxphos-produced ATP
PDH_1 = find(ismember(model.rxns, 'MAR08746')); %PDH step 1
Gln_a = find(ismember(model.rxns, 'MAR03804')); %Glutamate entry to TCA, NADP/NADPH
Gln_b = find(ismember(model.rxns, 'MAR03802')); %Glutamate entry to TCA, NAD/NADH

FApool = find(ismember(model.rxns, 'MAR09209')); 
PELDpool = find(ismember(model.rxns, 'MAR00656')); 
PSLDpool = find(ismember(model.rxns, 'MAR00661')); 
cholpool = find(ismember(model.rxns, 'MAR09285')); 
lipidpool = [FApool, PELDpool, PSLDpool, cholpool];

res_Gln = 0;
if keepGln == true
    model.c(:) = 0;
    model.c(Gln_a) = 1;
    sol4a = optimizeCbModel(model);
    model.c(:) = 0;
    model.c(Gln_b) = 1;
    sol4b = optimizeCbModel(model);
    if sol4a.f > 0 && sol4b.f > 0
        res_Gln = 1;
    end
else
    res_Gln = 1;
end

res_lipids = 0;
if keepLipids == true
    j = 0;
    for i = 1:numel(lipidpool)
        model.c(:) = 0;
        model.c(lipidpool(i)) = 1;
        sol = optimizeCbModel(model);
        if sol.f > 1e-5
            j = j+1;
        end
    end
    if j == numel(lipidpool)
        res_lipids = 1;
    end
else
    res_lipids = 1;
end

model_glc = setExchangeBounds(model, 'glucose', -1); %only glc allowed
res_ATP = 0;
if keepATP == true
    model_glc.c(:) = 0;
    model_glc.c(ATP) = 1;
    sol3a = optimizeCbModel(model_glc);
    model_glc.c(:) = 0;
    model_glc.c(PDH_1) = 1;
    sol3b = optimizeCbModel(model_glc);
    if sol3a.f > 1.9 && sol3b.f > 0.9 
        % glycolysis only: NGAM = 2
        res_ATP = 1;
    end
else
    res_ATP = 1;
end

model_glc = setExchangeBounds(model, {'glucose', 'O2'}, [-1, -1000]); %only glc and O2 allowed
res_oATP = 0;
if keepATPox == true
    model_glc.c(:) = 0;
    model_glc.c(ATP) = 1;
    sol2 = optimizeCbModel(model_glc);
    if sol2.f > 30 && sol2.x(oATP) > 25 && sol2.x(PDH_1) > 1.9
        % complete glycolysis + oxphos of 1 glucose: NGAM = 31.5, oATP = 27.5 and PHD_1 = 2
        res_oATP = 1;
    end
else
    res_oATP = 1;
end

end %end function checkFlag (local)
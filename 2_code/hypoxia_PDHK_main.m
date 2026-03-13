% main script for GEM simulations and analysis of hypoxia & PDHK1 project.
% requires COBRA toolbox and RAVEN toolbox.
% Rosemary Yu. Last update: 2026-03-13

%% init
initCobraToolbox
changeCobraSolver('gurobi','LP')
setRavenSolver('gurobi');

%% reduce HumanGEM to include only listed subsystems & reactions
% runtime ca. 30 min

clear
clc
load('3_models\Human-GEM.mat');
model = ihuman;

% subsystems and reactions to keep
subsysToKeep = {'Glycolysis / Gluconeogenesis'
                'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism'
                'Oxidative phosphorylation'
                'Pyruvate metabolism'
                }; 

rxnsToKeep = {'MAR09932' % maintenance biomass
              'MAR09034' % glucose exchange
              'MAR09048' % O2 exchange
              'MAR03964' % NGAM
              'MAR06916' % ATP from oxphos
              'MAR06328' % ATP transport between m and c
              'MAR08746' % PDH step 1
              'MAR20069' % PDH step 2
              'MAR06412' % PDH step 3
              'MAR06409' % PDH step 4
              'MAR04894' % oxygen transport between c and r
              'MAR04865' % malate transport from c to m (with Pi)
              'MAR09135' % lactate exchange
              'MAR05998' % lactate transport betwween c and e (with H+)
              'MAR04388' % lactate production from pyruvate
              'MAR09058' % CO2 exchange
              'MAR04919' % CO2 transport between c and e
              'MAR04922' % CO2 transport between c and m
              'MAR09044' % threonine exchange
              'MAR09047' % H2O exchange
              'MAR09067' % glycine exchange
              % keep NADH-NADPH conversion and glycolate-glyoxalate shunt:
              'MAR04271' % H+ [i] + NADH [m] + NADP+ [m] -> H+ [m] + NAD+ [m] + NADPH [m]
              'MAR08779' % H+ [m] + NADPH [m] + glyoxalate [m] -> NADP+ [m] + glycolate [m]
              'MAR09675' % glycolate [m] <=> glycolate [c]
              'MAR07702' % H+ [c] + NADPH [c] + glyoxalate [c] <=> NADP+ [c] + glycolate [c]
              'MAR07708' % glyoxalate [c] <=> glyoxalate [m]
              % keep Glutamine-Glutamate-AKG path:
              'MAR09063' % Glutamine exchange
              'MAR05078' % Glutamine transport between c and e
              'MAR01381' % Glutamine transport between c and m
              'MAR03892' % Glutamine -> glutamate (GLS)
              'MAR03804' % Glutamate entry to TCA (GLUD), NADP/NADPH
              'MAR03802' % Glutamate entry to TCA (GLUD), NAD/NADH
              % keep lipid exchange paths:
              'MAR00223' % H2O [c] + palmitoyl-CoA [c] => CoA [c] + H+ [c] + palmitate [c]
              'MAR00224' % ATP[c] + H2O[c] + palmitate[c] => ADP[c] + H+[c] + palmitate[e] + Pi[c]
              'MAR00611' % palmitate exchange
              'MAR09209' % fatty acid pool exchange
              'MAR00656' % PE-LD exchange
              'MAR00661' % PS-LD exchange
              'MAR09285' % cholesterol exchange
              };

% grab subsystems & reactions to close  
subsysToClose = readtable("1_raw_data\subsys_to_close.xlsx", ReadVariableNames=false);
subsysToClose = subsysToClose.Var1;
rxnsToClose = readtable("1_raw_data\rxns_to_close.xlsx", ReadVariableNames=false); 
rxnsToClose = rxnsToClose.Var1;

% main
tic
sub_model = getSubModel(model, subsysToKeep, rxnsToKeep, subsysToClose, rxnsToClose, ...
                        true, true, true, true);
toc

% save submodel
save('3_models\sub_model.mat','sub_model')
replaceMetNames = strcat(sub_model.metNames, {' ['}, sub_model.comps(sub_model.metComps), {']'});
sub_model.mets = replaceMetNames;
writeCbModel(sub_model,'3_models\sub_model.xlsx')

%% analysis submodel
load('3_models\sub_model.mat');
model = sub_model;

rxnNames = {'glucose'; 'O2'; 'lactate'; 'CO2'; 'glutamine'; 'palmitate';
            'ATPc'; 'ATPm'; 'PDH'};
rxns = {'MAR09034'; 'MAR09048'; 'MAR09135'; 'MAR09058'; 'MAR09063'; 'MAR00611'; 
        'MAR03964';'MAR06916';'MAR08746'};
rxnIndexes = getIndexes(model, rxns, 'rxns');

model = setExchangeBounds(model, {'glucose', 'O2'}, [-10, -60], [0, 0], false);
%model = setParam(model, 'lb', rxnIndexes(ismember(rxnNames, 'PDH')), 20);
model = setParam(model, 'obj', rxnIndexes(ismember(rxnNames, 'ATPc')), 1);
sol = optimizeCbModel(model);

lb = model.lb(rxnIndexes);
ub = model.ub(rxnIndexes);
obj = model.c(rxnIndexes);
fluxes = sol.x(rxnIndexes);

res = table(rxnNames, rxns, lb, ub, obj, fluxes);
disp(res)


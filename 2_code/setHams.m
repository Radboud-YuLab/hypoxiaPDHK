function exchModel = setHams(model)
% setHams
%
% Set a Hams culture medium for a humanGEM based model. Simplified from
% setHamsMedium
%
%   model           An ihuman-based GEM
%   exchModel       (struct) Model with Ham's media constraints
%
% Ivan Domenzain.      Last edited: 2019-11-29
% Rosemary Yu. 		   Last edited: 2026-02-16


%Remove unconstrained field, if available
if isfield(model,'unconstrained')
    model = rmfield(model,'unconstrained');
end

%Ham's media composition
mediaComps ={'glucose'
             'arginine'
             'histidine'
             'lysine'
             'methionine'
             'phenylalanine'
             'tryptophan'
             'tyrosine'
             'alanine'
             'glycine'
             'serine'
             'threonine'
             'aspartate'
             'glutamate'
             'asparagine'
             'glutamine'
             'isoleucine'
             'leucine'
             'proline'
             'valine'
             'cysteine'
             'thiamin'
             'hypoxanthine'
             'folate'
             'biotin'
             'pantothenate'
             'choline'
             'inositol'
             'nicotinamide'
             'pyridoxine'
             'riboflavin'
             'thymidine'
             'aquacob(III)alamin'
             'lipoic acid'
             'sulfate'
             'linoleate'
             'linolenate'
             'O2'
             'H2O'
             'retinoate'
             'Fe2+'
             'Pi'
             'alpha-tocopherol'
             'gamma-tocopherol'}; 
             
%Default flux bounds [LB, UB]
fluxBounds = [-ones(length(mediaComps),1), ones(length(mediaComps),1)]*1000;

%Set uptake fluxes for media mets
[exchModel,~] = setExchangeBounds(model,mediaComps,fluxBounds(:,1),fluxBounds(:,2),true);

%Check if model is feasible
sol = solveLP(exchModel); 
if ~isempty(sol.x)
	disp(['Constrained model is feasible'])
else
    disp(['Constrained model is unfeasible'])
	exchModel = [];
end
end
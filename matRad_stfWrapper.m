function stf = matRad_stfWrapper(ct,cst,Pln)
   stf = [];
   Fields = [];   
   
   %Generate stf for every modality
   for k=1:length(Pln.propStf)

      CurrStf = matRad_generateStf(ct,cst,Pln.OriginalPlans(k));
      Fields = [Fields; fieldnames(CurrStf)];
      stf = [stf,{CurrStf}];
      
   end
   
   %Add fields to the single stf structures when missing in order to make
   %them compatible
   TotFields = unique(Fields);
   for k=1:length(stf)
      IsPlanField = find(~isfield(stf{k},TotFields));
      if any(IsPlanField)
         for m=[IsPlanField]
            stf{k} = setfield(stf{k},{1},TotFields{m},[]);
         end
      end
   end
   stf = [stf{:}];
end

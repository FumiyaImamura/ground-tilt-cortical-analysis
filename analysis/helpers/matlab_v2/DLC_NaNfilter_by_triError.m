function NaNed_DLC=DLC_NaNfilter_by_triError(DLC, rec3d, DLC_triError_thresh_in_mm,options)
arguments
	DLC
	rec3d
	DLC_triError_thresh_in_mm
	options.Output (1,1) string {mustBeMember(options.Output, ["SameAsInput","AllNumeric_concatenatedIn2ndDimension","AllSeparatedInCellArray"])}="SameAsInput"
end


NaNed_DLC=[];
triError=DLC.triError;
flds=string(fields(DLC));
for k=flds(:)'
	if k=="triError"
		continue;
	end
	currDLC=DLC.(k);
	if iscell(currDLC) && length(currDLC)==length(rec3d.name_ind)	% 14 body parts in cell array
		NaNed_DLC.(k)=cell(size(currDLC));
		for i=1:numel(currDLC)
			coord=currDLC{i};
			ind=tri(i);
			currtri=triError(ind{:});
			currtri=repmat(currtri,1,3);
			if size(coord)~=size(currtri)
				error('triangulation array and coordinates array must have matching sizes')
			end
			trilog=currtri>DLC_triError_thresh_in_mm;	% indices of timepoints with high triangulation error
			coord(trilog)=NaN;
			NaNed_DLC.(k){i}=coord;
		end
		if options.Output=="AllNumeric_concatenatedIn2ndDimension"
			NaNed_DLC.(k)=cell2mat(NaNed_DLC.(k));
		end
	elseif isnumeric(currDLC) && size(currDLC,2)==16				% joint angles
		NaNed_DLC.(k)=cell(1,size(currDLC,2));
		allIND=arrayfun(rec3d.ind_name,rec3d.Dotname_ind);
		for i=1:size(allIND,1)
			ind=tri(allIND(i,:));
			currtri=triError(ind{:});
			trilog=currtri<=DLC_triError_thresh_in_mm;	% indices of timepoints with low error
			trilog=logical(prod(trilog,2));			% all 3 body parts should have high accuracy
			
			ind=tri(i);
			coord=currDLC(ind{:});
			if size(coord)~=size(trilog)
				error('triangulation array and trilog array must have matching sizes')
			end
			coord(~trilog)=NaN;
			NaNed_DLC.(k){i}=coord;
		end
		if ~(options.Output=="AllSeparatedInCellArray")
			NaNed_DLC.(k)=cell2mat(NaNed_DLC.(k));
		end
	else
		error(append("The field ",k," has weird shape or structure"))
	end
end


	function IND=tri(ind)
		IND=repmat({':'},1,ndims(triError));
		IND{2}=ind(:)';
	end

end
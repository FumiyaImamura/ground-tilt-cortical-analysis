function NaNed_DLC=selectiveDLC_NaNfilter_by_triError(DLC_cell_or_numeric_2ndDimIsParams,triError, rec3d, DLC_triError_thresh_in_mm,options)
arguments
	DLC_cell_or_numeric_2ndDimIsParams
	triError
	rec3d
	DLC_triError_thresh_in_mm
	options.XYZOrScalar	(1,1) {mustBeMember(options.XYZOrScalar,["XYZ","Scalar"])}
	options.Output (1,1) string {mustBeMember(options.Output, ["AllNumeric_concatenatedIn2ndDimension","AllSeparatedInCellArray"])}="SameAsInput"
end

DLC=DLC_cell_or_numeric_2ndDimIsParams;
NaNed_DLC=[];
if ~isnumeric(triError)
	error("triError must be numeric")
end

hightrilog=generate_high_triError_log();
hightrilog=divideArray2ChunksOfCells_2ndDim(hightrilog,1);
lowtrilog=generate_low_triError_log();
lowtrilog=divideArray2ChunksOfCells_2ndDim(lowtrilog,1);

data=DLC;
switch options.XYZOrScalar
	case "XYZ"
		if iscell(data) && length(data)==length(rec3d.name_ind)
		elseif isnumeric(data) && size(data,2)==3*length(rec3d.name_ind)
			data=divideArray2ChunksOfCells_2ndDim(data,3);
		else
			error("Size must match")
		end

		hightrilog=cellfun(@(a)repmat(a,1,3),hightrilog,UniformOutput=false);
		for k=1:length(data)
			data{k}(hightrilog{k})=NaN;
		end

	case "Scalar"
		if isnumeric(data)
			data=divideArray2ChunksOfCells_2ndDim(data,1);
		end

		if (iscell(data) && length(data)==size(rec3d.Dotname_ind,1))
			hightrilog=[];
			for k=1:length(data)
				currinds=rec3d.ind_name(rec3d.Dotname_ind(k,:));
				currtrilog=lowtrilog(currinds);
				currtrilog=cell2mat(currtrilog);
				currtrilog=prod(currtrilog,2);
				hightrilog{k}=~logical(currtrilog);
			end
			for k=1:length(data)
				data{k}(hightrilog{k})=NaN;
			end
		elseif (iscell(data) && length(data)==length(rec3d.name_ind))
			for k=1:length(data)
				data{k}(hightrilog{k})=NaN;
			end
		else
			error("Size must match")
		end
		
end
NaNed_DLC=data;
if ~(options.Output=="AllSeparatedInCellArray")
	NaNed_DLC=cell2mat(NaNed_DLC);
end


	function hightrilog=generate_high_triError_log()	%
		hightrilog=[];
		if any(isnan(triError));error("There cannot be nan in triError");end
		hightrilog=zeros(size(triError));
		hightrilog(triError>DLC_triError_thresh_in_mm)=1;
		hightrilog=logical(hightrilog);
	end
	function lowtrilog=generate_low_triError_log()	%
		lowtrilog=[];
		if any(isnan(triError));error("There cannot be nan in triError");end
		lowtrilog=zeros(size(triError));
		lowtrilog(triError<=DLC_triError_thresh_in_mm)=1;
		lowtrilog=logical(lowtrilog);
	end

	function celldata=divideArray2ChunksOfCells_2ndDim(numeric,chunksize)
		if size(numeric,2)/chunksize~=round(size(numeric,2))
			error("the size must be divisible by chunk size")
		end
		celldata=[];
		for i=1:size(numeric,2)/chunksize
			IND=repmat({':'},1,ndims(numeric));
			IND{2}=[chunksize*(i-1)+1:chunksize*i];
			celldata{i}=numeric(IND{:});
		end
	end


end
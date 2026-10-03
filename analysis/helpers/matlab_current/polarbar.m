function polarbar(dataVec, options,opts)
    arguments
        dataVec (:,1) double                                    % input angles, in radians
		options.ThetaDir (1,1) string {mustBeMember(options.ThetaDir,["clockwise","counterclockwise"])} = "counterclockwise"
		options.ThetaZeroLocation (1,1)  = "right"
        opts.?matlab.graphics.chart.primitive.Histogram     % all PolarHistogram name–value args
	end

	edges=linspace(0,2*pi,numel(dataVec)+1);
	bincenter=edges(1:end-1)+diff(edges,1)/2;
	bincenter=num2cell(bincenter);
	dataVec=num2cell(dataVec(:)');
	valuelist=cellfun(@(a,b)repmat(b,a,1),dataVec,bincenter,UniformOutput=false);
	valuelist=cat(1,valuelist{:});
	opts=namedargs2cell(opts);
	polarhistogram(valuelist,edges,opts{:})

	ax = gca;
	ax.ThetaDir           = options.ThetaDir;
	ax.ThetaZeroLocation = options.ThetaZeroLocation;
end

function [GaussHotEncoded_Data,DiscretizingFunc_cellVector]=GaussianHotEncode_SpecifedSPparams(tpk,edges_cellVector,sigma,params2discretize_stringArray4jointDiscretization, options)
% GaussHotEncoded_Data is a cell array when Output="CellArray",
% concatenated to a numeric array when Output="NumericArray_concatenateOutputIn2ndDimension"

arguments
    tpk
    edges_cellVector
	sigma
end
arguments (Repeating)
    params2discretize_stringArray4jointDiscretization (1,:) string {mustBeMember(params2discretize_stringArray4jointDiscretization, ...
                                          ["theta","phi","tau","X","Y","Z"])}
end
arguments
    options.Output (1,1) string {mustBeMember(options.Output, ["CellArray","NumericArray_concatenateOutputIn2ndDimension"])}="CellArray"
end

[Discretized_t1k_cellVector, DiscretizingFunc_cellVector ]...
	= discretizeSpecifiedSPparams(tpk, edges_cellVector, params2discretize_stringArray4jointDiscretization{:});

GaussHotEncoded_Data=cell(1,length(Discretized_t1k_cellVector));
for i=1:length(Discretized_t1k_cellVector)
	Data_numeric				=Discretized_t1k_cellVector{i};
	encoding_dim				=1;
	parameter_dim_mustBeSize1	=2;
	sigma						=sigma;
	bincount					=DiscretizingFunc_cellVector{i}.bincount;
	Data_disc_numeric = GaussianHotEncoder_usingDiscretizedData(Data_numeric, encoding_dim, parameter_dim_mustBeSize1, sigma, bincount);
	GaussHotEncoded_Data{i}=Data_disc_numeric;
end

if options.Output=="NumericArray_concatenateOutputIn2ndDimension"
	cell2mat(GaussHotEncoded_Data);
end


end

%%
function testEncodedData()
%% for data that was jointly encoded
GaussHotEncoded_Data;
tpk;

ind=1;	% index of GaussHotEncoded_Data that was encoded jointly
ind2plot=1:2;	% indices of original tpk (second dimension) used for above joint encoding

IND=repmat({':'},1,ndims(tpk));
IND{2}=ind2plot;
orgdata=tpk(IND{:});
gausEnc=GaussHotEncoded_Data{ind};

f=@(a)BlockCat(reshape(a,size(a,1),size(a,2),[]),3,1);
orgdata=f(orgdata);
gausEnc=f(gausEnc);

figure;
CSP.figureTile;
ax1=nexttile;
plot(orgdata);
ax2=nexttile;
plot(gausEnc);
linkaxes([ax1,ax2],'x')
%% testing gaussian hot encoding: plotting the one that was jointly encoded
spname="p8_tpkd";
DLCname="p8_DLC";

A=CSP(info,plane=1);
A.load_specified(1,spname,DLCname);
tpk=A.(spname);


%%%%%%%%%
[edges_cellVector,bincenter_cellVector]=makeSPBins();
GaussHotEncoded_Data=GaussianHotEncode_SpecifedSPparams(tpk,edges_cellVector,1,["theta","phi"],Output="CellArray");

ind=1;	% index of GaussHotEncoded_Data that was encoded jointly
ind2plot=1:2;	% indices of original tpk (second dimension) used for above joint encoding

IND=repmat({':'},1,ndims(tpk));
IND{2}=ind2plot;
orgdata=tpk(IND{:});
gausEnc=GaussHotEncoded_Data{ind};

f=@(a)BlockCat(reshape(a,size(a,1),size(a,2),[]),3,1);
orgdata=f(orgdata);
gausEnc=f(gausEnc);

CSP.figureTile;
ax1=nexttile;
plot(orgdata);
ax2=nexttile;
plot(gausEnc);
linkaxes([ax1,ax2],'x')








%% for data that was individually encoded
ind=2;	% index of GaussHotEncoded_Data that was not encoded jointly
ind2plot=1:2;	% indices of original tpk (second dimension) used for above joint encoding

GaussHotEncoded_Data;
tpk;
paramIndex=2;
edges;		% create this on your own

plotGaussianHotEncodedData(tpk, GaussHotEncoded_Data, edges, paramIndex, options);

%% testing gaussian hot encoding: plot parameters that were not jointly discretized
spname="p8_tpk";
DLCname="p8_DLC";

A=CSP(info,plane=1);
A.load_specified(1,spname,DLCname);
tpk=A.(spname);


%%%%%%%%%
	[edges_cellVector,bincenter_cellVector]=makeSPBins();
GaussHotEncoded_Data=GaussianHotEncode_SpecifedSPparams(tpk,edges_cellVector,1,"theta","phi","tau","X","Y","Z",Output="CellArray");

ind=1;	% index of GaussHotEncoded_Data that to plot
edges=edges_cellVector;

figure;
plotGaussianHotEncodedData(tpk, GaussHotEncoded_Data, edges, ind);


end
function StrChunks=splitString2Chunks(inputStr,delim_StrArray,options)
arguments
	inputStr	string
	delim_StrArray (1,:) string
	options.direction string {mustBeMember(options.direction,["first","last"])}="last"
end
anyToken="(.*)";
delim=cellstr(delim_StrArray);
delim=cellfun(@(a)regexptranslate("escape",a),delim,UniformOutput=false);
Searcher=delim(:)';
switch options.direction
	case "first"
		for i=1:length(delim)
			if i==1
				Searcher{i}=append(NonStrToken(delim{i}),delim{i},Non2StrToken(delim{i},delim{i+1}));
			elseif i==length(delim)
				Searcher{i}=append(delim{i},anyToken);
			else
				Searcher{i}=append(delim{i},Non2StrToken(delim{i},delim{i+1}));
			end
		end
	case "last"
		for i=1:length(delim)
			if i==1
				Searcher{i}=append(anyToken,delim{i});
			elseif i==length(delim)
				Searcher{i}=append(Non2StrToken(delim{i},delim{i-1}),delim{i},NonStrToken(delim{i}));
			else
				Searcher{i}=append(Non2StrToken(delim{i},delim{i-1}),delim{i});
			end
		end
end
Searcher=append(Searcher{:});
StrChunks=regexp(inputStr,Searcher,"tokens","once");
	function newstr=Non2StrToken(str1,str2)
		newstr=NonStrToken(append(str1,str2));
	end
	function newstr=NonStrToken(str)
		newstr=append("([^",str,"])*");
	end
end
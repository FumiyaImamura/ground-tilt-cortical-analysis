function [Trainer, Tester] = splitPredResp(Pred, Resp, groupnum,dim2split)
% splitPredResp splits high-dimensional arrays Pred and Resp into cross-validation groups.
%
%   [Trainer, Tester] = splitPredResp(Pred, Resp, groupnum)
%
%   INPUTS:
%       Pred     - A high-dimensional array where:
%                  Dimension 1: Time.
%                  Dimension 2: Predictors.
%                  Dimensions 3 and higher: Trials and conditions.
%
%       Resp     - A high-dimensional array with the same number of dimensions as Pred.
%                  The sizes in dimensions 1 and 3+ must match Pred, though dimension 2
%                  (responses) may differ.
%
%       groupnum - A positive integer specifying the number of groups (folds) for splitting.
%
%   OUTPUTS:
%       Trainer  - A structure with the following fields:
%                   Trainer.Pred: A cell array (length = groupnum) where each cell contains
%                                 training data from Pred.
%                   Trainer.Resp: A cell array (length = groupnum) where each cell contains
%                                 training data from Resp.
%
%       Tester   - A structure with the following fields:
%                   Tester.Pred: A cell array (length = groupnum) where each cell contains
%                                testing data from Pred.
%                   Tester.Resp: A cell array (length = groupnum) where each cell contains
%                                testing data from Resp.
%
%   The function works as follows:
%       1. It verifies that Pred and Resp have at least 3 dimensions and that their
%          dimensions (except the second) match.
%       2. It reshapes both arrays so that dimensions 4 and higher are collapsed into
%          the third dimension.--> actually, no collapsing
%       3. It splits the third dimension into groupnum groups by taking every groupnum-th
%          element as a test set and the remaining as the training set for that fold.
%
%   Example:
%       [Trainer, Tester] = splitPredResp(Pred, Resp, 5);
%
arguments
	Pred
	Resp
	groupnum
	dim2split
end

    %% Check inputs
    dimsPred = size(Pred);
    dimsResp = size(Resp);
    
    % if numel(dimsPred) < 3
    %     error('Pred must have at least 3 dimensions.');
    % end
    % if numel(dimsResp) < 3
    %     error('Resp must have at least 3 dimensions.');
    % end
    
    % Check that dimensions 1 and 3 onward match between Pred and Resp.
    if ~isequal(dimsPred([1, 3:end]), dimsResp([1, 3:end]))
        error('Pred and Resp dimensions mismatch: dimensions 1 and 3 onward must be the same.');
    end
    
    if ~isscalar(groupnum) || groupnum < 1 || mod(groupnum, 1) ~= 0
        error('groupnum must be a positive integer scalar.');
    end

 %% Split the data into groups.
 for i=1:groupnum
	 ind=[i:groupnum:size(Pred,dim2split)];
	 exind=setdiff(1:size(Pred,dim2split),ind);

	 IND=repmat({':'},1,ndims(Pred));IND{dim2split}=ind;
	 testPred{i}=Pred(IND{:});
	 IND{dim2split}=exind;
	 trainPred{i}=Pred(IND{:});

	 IND=repmat({':'},1,ndims(Pred));IND{dim2split}=ind;
	 testResp{i}=Resp(IND{:});
	 IND{dim2split}=exind;
	 trainResp{i}=Resp(IND{:});
 end

 %% Package outputs into structures.
 Trainer.Pred = trainPred;
 Trainer.Resp = trainResp;

 Tester.Pred  = testPred;
 Tester.Resp  = testResp;

	function collapse_3rd_dim_and_beyond(data)
		data=reshape(pagetranspose(data),size(data,2),[])';
	end

end

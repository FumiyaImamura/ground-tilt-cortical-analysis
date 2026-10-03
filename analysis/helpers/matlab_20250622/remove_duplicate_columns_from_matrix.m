function [X_new, keep_idx] = remove_duplicate_columns_from_matrix(X)
    % This function removes columns that are scalar multiples of one another.
    % For each column, we adjust the sign (so the first nonzero is positive)
    % and normalize the column. Then, we use unique() with a rounding tolerance
    % to account for numerical precision. In groups of columns that are multiples,
    % only one representative is kept.
    
    [n, m] = size(X);
    X_norm = zeros(n, m);
    tol = 1e-10;  % tolerance for floating-point comparisons
    
    for j = 1:m
        col = X(:, j);
        if norm(col) == 0
            % If the column is a zero vector, leave it unchanged.
            X_norm(:, j) = col;
        else
            % Find the first nonzero element and adjust sign so that it is positive.
            idx = find(col, 1, 'first');
            if col(idx) < 0
                col = -col;
            end
            % Normalize the column.
            X_norm(:, j) = col / norm(col);
        end
    end
    
    % To reduce numerical issues, round the normalized columns.
    X_round = round(X_norm/tol) * tol;
    
    % Each column is now a row in X_round'. Use unique to find unique rows.
    [~, keep_idx, ~] = unique(X_round', 'rows', 'stable');
    
    % Extract the columns to be kept.
    X_new = X(:, keep_idx);
end

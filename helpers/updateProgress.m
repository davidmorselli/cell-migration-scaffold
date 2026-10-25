function updateProgress(total, reset)
    persistent completed

    if isempty(completed) || reset
        completed = 0;
    end

    if reset
        return
    end

    completed = completed + 1;
    fprintf('Completed %d/%d (%.1f%%)\n', ...
        completed, total, 100 * completed / total);
end
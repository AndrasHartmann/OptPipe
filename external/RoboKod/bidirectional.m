function [is_bidirectional] = bidirectional(values)
    is_bidirectional = any(values(values < 0)) && any(values(values > 0));
end
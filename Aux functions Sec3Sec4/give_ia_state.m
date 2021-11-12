function [ia_state] = give_ia_state(state)
    ia_state = state(find(triu(ones(size(state,1)))));
end
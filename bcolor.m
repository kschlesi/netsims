function h = bcolor(inmat)
% provides a balanced color plot (no row/cols left out) with no edge lines
    if ~ismatrix(inmat)
        error('input matrix must be two-dimensional'); 
    end
    
    pad = mean(mean(inmat));
    h = pcolor(padarray(inmat,[1 1],pad,'post'));
    set(h, 'EdgeColor', 'none');
    
end
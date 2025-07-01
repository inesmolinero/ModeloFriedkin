% guardarFigura.m
function guardarFigura(fig,outDir,fn)
    try saveas(fig,fullfile(outDir,fn)); catch, set(fig,'Visible','on'); pause(1); end
    close(fig);
end
function s = safeName(tag)
    s = regexprep(tag,'[^a-zA-Z0-9_]','_');
end

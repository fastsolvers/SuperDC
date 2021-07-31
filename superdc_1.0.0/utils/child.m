function ch = child(tr)

ch = {[]};

for i = 1:length(tr)-1
    ch{tr(i)} = [];
end

for i = 1:length(tr)-1
    ch{tr(i)} = [ch{tr(i)} i];
end
% read results file
fid = fopen(strcat(pwd(),'/Traffic_network_framework/results.txt'));
nedges = 76;

%%
thisline = fgetl(fid);
cct = 0; capacity = []; labels = {}; f = zeros(nedges,1);
while thisline>0
    bfind = strfind(thisline,'[');
    thisline = strrep(thisline,']','');
    if ~isempty(bfind)
        cct = cct +1;
        commas = strfind(thisline,',');
        capacity(cct) = str2double(thisline(1:commas(1)-1));
        labels{cct} = strrep(thisline(commas(1)+1:commas(2)-1),' ','');
        thisline = thisline(bfind+1:end);
        f(:,cct) = 0; flowct = 1;
    end
    
    spfind = strfind(thisline,' ');
    start = 1; 
    for i = 1:length(spfind)
        thisnum = str2double(thisline(start:spfind(i)));
        if ~isnan(thisnum)
             f(flowct,cct) = thisnum;
             flowct = flowct+1;
        end
        start = spfind(i);
    end
    thisnum = str2double(thisline(start:end));
    if ~isnan(thisnum)
        f(flowct,cct) = thisnum;
        flowct = flowct+1;
    end
    
    thisline = fgetl(fid);
end

% keep flow solution only
solct = 0;
fsol = zeros(nedges,1);
csol = [];
for i=1:length(labels)
    thislabel = labels{i};
    if contains(thislabel,'upper')
    else
        solct = solct + 1;
        csol(solct) = capacity(i);
        fsol(:,solct) = f(:,i);
    end
end

fold = f;
f = fsol;
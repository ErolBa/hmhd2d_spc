function nameshift(pre,ext,si,ei,offset)
% shift the file name number
% Input parameter :
% pre    : prefix of the filename
% ext    : extension of the filename
% si     : starting index
% ei     : end index
% offset : offset number


if si>ei 
    ss=si;
    si=ei;
    ei=ss;
end

if offset==0 
    return
elseif offset>0
    if (ei+offset)>999 
        return
    else
        ii=ei:-1:si;    % large number first
    end
else
    if (si+offset)<0 
        return
    else
        ii=si:1:ei;     % small number first
    end
end
    
for i=ii
    fname=sprintf([pre,'%03i','.',ext],i);
    fname1=sprintf([pre,'%03i','.',ext],i+offset);
    unix(['mv ',fname,' ',fname1]);
end

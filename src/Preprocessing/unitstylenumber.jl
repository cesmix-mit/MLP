function unitstylenumber(style::String)

if lowercase(style) == "lj" 
    unitstylenum = 0;
elseif lowercase(style) == "real"
    unitstylenum = 1;    
elseif lowercase(style) == "metal"
    unitstylenum = 2;
elseif lowercase(style) == "si"
    unitstylenum = 3;
elseif lowercase(style) == "cgs"
    unitstylenum = 4;
elseif lowercase(style) == "electron"
    unitstylenum = 5;
elseif lowercase(style) == "micro"
    unitstylenum = 6;
elseif lowercase(style) == "nano"
    unitstylenum = 7;
else
    error("Invalid unit style");
end

return unitstylenum

end


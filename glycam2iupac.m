function IUPACout = glycam2iupac(glycam)
IUPACout = glycam;
bondlinks = unique(regexp(glycam,'(a|b)[12]-[0-9|OH|OME|OtBu]+','match'));
for i = 1:length(bondlinks)
    bondlinksnew = ['(',bondlinks{i},')'];
    IUPACout = strrep(IUPACout,bondlinks{i},bondlinksnew);
end
specmsexp = 'DNeup5Ac|DNeup5Gc|KDN|KDO|LGalpNAc|DGalpNAc|LGlcpNAc|DGlcpNAc|LManpNAc|DManpNAc|LGalpA|DGalpA|LGlcpA|DGlcpA|LIdopA|DIdopA';
specms = unique(regexp(glycam,specmsexp,'match'));
specmsori = {'DNeup5Ac','DNeup5Gc','KDN','KDO','LGalpNAc','DGalpNAc','LGlcpNAc','DGlcpNAc','LManpNAc','DManpNAc','LGalpA','DGalpA','LGlcpA','DGlcpA','LIdopA','DIdopA'};
specmsnew = {'Neu5Ac','Neu5Gc','KDN','KDO','GalNAc','GalNAc','GlcNAc','GlcNAc','ManNAc','ManNAc','GalA','GalA','GlcA','GlcA','IdoA','IdoA'};
specms = specmsori(ismember(specmsori,specms));
newspecms = specmsnew(ismember(specmsori,specms));
for i = 1:length(specms)
    IUPACout = strrep(IUPACout,specms{i},newspecms{i});
end
regmsexp = 'Man|Gal|Glc|Ido|All|Alt|Gul|Tal|Xyl|Lyx|Rib|Ara|Fru|Psi|Sor|Tag|Fuc|Rha|Qui';
regms = unique(regexp(glycam,['(L|D)(',regmsexp,')(f|p)'],'match'));
for i = 1:length(regms)
    IUPACout = strrep(IUPACout,regms{i},regms{i}(2:end-1));
end
end
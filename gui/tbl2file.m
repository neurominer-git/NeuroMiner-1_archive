function [ERR, STATUS, fil, typ] = tbl2file(tbl, filename, sheetname )
% Write table data to either a Excel or a text-based file
ERR=[]; STATUS = 0;
try
    if ispc 
        fil = sprintf('%s.xls',filename); typ='xls';
        s_rownames = size(tbl.rownames);
        if s_rownames(1)==1,
            tbl.rownames = tbl.rownames';
        end
        Ix = isnan(tbl.array); tbl_array = num2cell(tbl.array); tbl_array(Ix)={[]};
        tbl = [tbl.colnames; [tbl.rownames tbl_array]];
        STATUS = xlswrite(fil, tbl, sheetname);                 
    else
        fil = sprintf('%s_%s.txt', filename, sheetname); typ='txt';
        fid = fopen(fil,'w');
        %Print header row
        for i=1:size(tbl.colnames,2)
            fprintf(fid,'%s',tbl.colnames{i});
            if i<size(tbl.colnames,2),fprintf(fid,'\t'); end
        end
        %Print table rows
        for i=1:numel(tbl.rownames)
            fprintf(fid,'\n%s', tbl.rownames{i});
            for j=1:size(tbl.array,2)
                fprintf(fid,'\t%g', tbl.array(i,j));
            end
        end
        fclose(fid);
        STATUS = 1;
    end
catch ERR
    errordlg(ERR.message);
end
   
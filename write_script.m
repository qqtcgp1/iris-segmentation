function write_script(filename, variable)
filename1=strcat(filename,'.txt');
filename2=strcat(filename,'.py');


fid=fopen(filename2,'w');
fprintf(fid, '(');

if size(variable,2)==2
    for i=1:length(variable)-1
        fprintf(fid,'(%.10f,%.10f',variable(i,1),variable(i,2));
        fprintf(fid,'),');
    end
    fprintf(fid,'(%.10f,%.10f',variable(end,1),variable(end,2));
    fprintf(fid,'))');
else
    for i=1:length(variable)-1
        fprintf(fid,'(%.10f,%.10f,%.10f',variable(i,1),variable(i,2),variable(i,3));
        fprintf(fid,'),');
    end
    fprintf(fid,'(%.10f,%.10f,%.10f',variable(end,1),variable(end,2),variable(end,3));
    fprintf(fid,'))');
end

fclose(fid);
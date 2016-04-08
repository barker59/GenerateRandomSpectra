function fprintsetup(SaveFileName,Filename)
[myfile,message] = fopen(SaveFileName ,'wt');
if myfile < 0
    disp(message);
else
fprintf(myfile,['Author:,' 'Aaron Barker\n']);
fprintf(myfile,['File Created:,',datestr(now),'\n']);
fprintf(myfile,['Filename:',',',Filename,'\n,\n' ]);
fclose(myfile);
end
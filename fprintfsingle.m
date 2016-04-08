function fprintfsingle(SaveFileName,MonthlyWindows,Probability,Length,ParameterName,ParameterValue,Unit)
DataToWrite=[Probability];
myfile = fopen(SaveFileName ,'a');
fprintf(myfile,['Duration:, ',Length ',','hours\n']);
fprintf(myfile,[ParameterName ',' ,ParameterValue ',' ,Unit '\n']);
fprintf(myfile,['Month:,','January,','February,','March,','April,','May,','June,','July,','August,','September,','October,','Novemeber,','December\n']);
fprintf(myfile,['Probability:,']);
dlmwrite(SaveFileName,DataToWrite, '-append','roffset',0,'delimiter',',');
fprintf(myfile,',\n');
fclose(myfile);

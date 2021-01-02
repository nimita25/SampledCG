filename = ["G1";
    "G2";
    "G3";
    "G4";
    "G5";
    "G14";
    "G15";
    "G16";
    "G17";
    "G22";
    "G23";
    "G24";
    "G25";
    "G26";
    "G35";
    "G36";
    "G37";
    "G38";
    "G43";
    "G44";
    "G45";
    "G46";
    "G47";
    "G48";
    "G49";
    "G50";
    "G51";
    "G52";
    "G53";
    "G54";
    "G55";
    "G58";
    "G60";
    "G63"];

filename = ["G1";
    "G2";
    "G3";
    "G4";
    "G5";
    "G14";
    "G15";
    "G16";
    "G17"];

for i = 1:length(filename)
    profile clear
    profile -memory on
    
    data = strcat('large_instances/',filename(i),'.mat');
    disp(data);
%     [v,x,s,so] = Penalty('G',data,0.1,3600,0.1);
%     fileID = fopen('Output/Memory_Algo2.txt','a');
    
    [f,mc] = Penalty_SeDuMi('G',data);
    fileID = fopen('Output/Memory_SeDuMi.txt','a');
%     
%     [f,mc] = Penalty_SDPT3('G',data,1e-3);
%     fileID = fopen('Output/Memory_SDPT3.txt','a');
%     
%     [f,mc] = Penalty_SDPNAL('G',data,1e-3);
%     fileID = fopen('Output/Memory_SDPNAL.txt','a');
    
%     data = strcat('G/',filename(i));
%     data = convertStringsToChars(data);
%     disp(data);
%     Test_MaxCut_CGAL(data);
%     fileID = fopen('Memory_CGAL10.txt','a');
    
    
    p = profile('info');
    memoryUsed = max([p.FunctionTable.PeakMem]);
    fprintf(fileID,'\n');
    fprintf(fileID,'%s&', filename(i));
    fprintf(fileID,'%.3f', memoryUsed/1024);
    fclose(fileID);
    profile clear

end

    
clear all;
clc;

%% Initial conditions

n_node=155; % Number of nodes

% C1=Calcite        IC=1 BC=0
% C2=CO2            IC=0 BC=1
% C3=C4             IC=0 BC=0
% C4=Ca             IC=0 BC=0

C1(1:n_node,2)=1;
C1(1:n_node,1)=[1:n_node];
C2(1:n_node,2)=0;
C2(1:n_node,1)=[1:n_node];
C3(1:n_node,2)=0;
C3(1:n_node,1)=[1:n_node];
C41(1:n_node,2)=0;
C41(1:n_node,1)=[1:n_node];

fid = fopen( 'Initil_Condition_C1.txt', 'wt' );
fprintf(fid,'%f %f \n', [C1]');
fclose(fid);
fid = fopen( 'Initil_Condition_C2.txt', 'wt' );
fprintf(fid,'%f %f \n', [C2]');
fclose(fid);
fid = fopen( 'Initil_Condition_C3.txt', 'wt' );
fprintf(fid,'%f %f \n', [C3]');
fclose(fid);
fid = fopen( 'Initil_Condition_C4.txt', 'wt' );
fprintf(fid,'%f %f \n', [C41]');
fclose(fid);

% Transfer initial condition files from API folder to ic-ferst folder
system('cp ../phreeqc/Initil_Condition_C1.txt ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional');
system('cp ../phreeqc/Initil_Condition_C2.txt ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional');
system('cp ../phreeqc/Initil_Condition_C3.txt ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional');
system('cp ../phreeqc/Initil_Condition_C4.txt ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional');


T_total=0.1; %Total time
delta_t=10e-3; % dt

steps=T_total/delta_t;

for step=1:steps
    clear readData;
    clear Par;
    %clear all;
    %% RUNNIG IC-FERST
    
    system('../phreeqc ./expcom.sh'); 								%path of folder of matlab file
    ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional 		%path of folder of Ic-ferst file
    system('./RUNICFERST.sh');
    
    
    %% Convert export file from IC-FERST (.VTU) to CSV file
    
    system('cp ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional/porous_compositional_1.vtu ../phreeqc/');
    
    system('sudo python code.py');
    
    
    %% Call CSV file in MATLAB
    
    % Open data file from IC-FERST
    
    fid = fopen('porous_compositional_1.0.csv','rt');
    headerline = fgetl(fid);
    headerfields = regexp(headerline, ',', 'split');
    numcol = length(headerfields);
    %numrow = count_remaining_lines(fid);
    frewind(fid);
    fgetl(fid);    %re-read and discard the headerline
    
    
    % Read data in form CSV file
    readData =textscan(fid,repmat('%f',1,numcol),'HeaderLines',1,'Delimiter',',');
    
    %---------------------------------------------------------------------
    
    N=numcol;
    for i= 1:N
        Par(:,i)=readData{1,i}(:,1);
    end
    % Remove Raw Rows
    ff=size(Par);
    for ii=1:n_node+1
        k=1;
        Row_Rem=[];
        Ax=Par(ii,N-2);
        Ay=Par(ii,N-1);
        ff=size(Par);
        for jj=2:ff(1,1)
            Bx=Par(jj,N-2);
            By=Par(jj,N-1);
            if (Ax==Bx) && (Ay==By)
                Row_Rem(k)=jj;
                k=k+1;
                ff=size(Par);
            end
            ff=size(Row_Rem);
        end
        m=0;
        for ss=1:ff(1,2)-1
            Par(Row_Rem(ss+1)-m,:) = [];
            m=m+1;
        end
    end
    Par(1,61)= 1;
    
    % Make the value vector for PREEQCRM
    
    fid1 = fopen('Mesh.txt','rt');
    headerline = fgetl(fid1);
    headerfields = regexp(headerline, ',', 'split');
    numcol = 4;
    frewind(fid1);
    fgetl(fid1);
    MeshData =textscan(fid1,repmat('%f',1,numcol),'HeaderLines',1,'Delimiter',',');
    Component1MassFraction=[];
    for ic=1:n_node
        for jc=1:n_node
            if (round(Par(ic,N-2),5)==round((MeshData{1,2}(jc,1)),5)) && ((round(Par(ic,N-1),5))==round((MeshData{1,3}(jc,1)),5))
                Phase1density(jc,1)=Par(ic,2)/1000;
                Component1MassFraction(jc,1)=Par(ic,5);
                Component2MassFraction(jc,1)=Par(ic,6);
                Component3MassFraction(jc,1)=Par(ic,8);
                Component4MassFraction(jc,1)=Par(ic,9);
            end
        end
    end
    Component1MassFraction(3,1)=Component1MassFraction(4,1);
    Component2MassFraction(3,1)=Component1MassFraction(4,1);
    Component3MassFraction(3,1)=Component1MassFraction(4,1);
    Component4MassFraction(3,1)=Component1MassFraction(4,1);
    
    %% making the input PhreeqcRM file (advect.pqi) and Run 

    node=1
    time_step=step
    Component1M(1:n_node,1)=1:n_node;
    Component2M(1:n_node,1)=1:n_node;
    Component3M(1:n_node,1)=1:n_node;
    Component4M(1:n_node,1)=1:n_node;
    for i=1:n_node
        
        T=50;
        pH=7;
        pe=4;
        density=Phase1density(i,1);
        CO2=Component2MassFraction(i,1);
        Calcite=Component1MassFraction(i,1);
        C4=Component3MassFraction(i,1) ;
        Ca=Component4MassFraction(i,1) ;
        fid = fopen( 'advect.pqi', 'wt' );
        formatSpec1 = 'SOLUTION 1 Pure water \n temp %8.3f \n pH %8.3f \n pe %8.3f \n redox     pe\n unit     mmol/kgw \n density %8.3f\n C(4)%8.3f \n Ca %8.3f \n -water    1 # kg \n\n EQUILIBRIUM_PHASES 1 \n CO2(g)    0 %8.3f \n Calcite   0 %8.3f \n \n \n END \n \n \n SELECTED_OUTPUT \n -file            advect.sel \n -reset           false \n -totals          C Ca';
        fprintf(fid,formatSpec1,T,pH,pe,density,C4,Ca,CO2,Calcite);
        fclose(fid);
        
        %% RUN PhreeqcRM
        
        system('./expcom.sh');
        time_step= step
        node=node+1
        
        %% Get the data from phreeqcRM
        
        fileID = fopen('Components_cpp.txt','rt');
        formatSpec = '%f';
        PhreeqcRMdata = fscanf(fileID,formatSpec);    
        
        %-------->
        Component1M(i,2)=Component1MassFraction(i,1);
        Component2M(i,2)=Component2MassFraction(i,1);
        Component3M(i,2)=PhreeqcRMdata(5,1);
        Component4M(i,2)=PhreeqcRMdata(6,1);
                fid = fopen('Advect_cpp.chem.txt' );
        S = textscan(fid,'%s','delimiter','\n') ;
        S = S{1}  ;
        fclose(fid) ;
        
        idxXi = strfind(S, 'pH');
        idxXi = find(not(cellfun('isempty',idxXi)));
        
        idxYi = strfind(S, 'Calcite');
        idxYi = find(not(cellfun('isempty',idxYi)));
        
        idxYi2 = strfind(S, 'Calcite');
        idxYi2 = find(not(cellfun('isempty',idxYi2)));
        
        A=size(idxYi);
        pH(i+1,2)= str2num(S{idxXi(3,1)}(9:13)) ;
        Calcite_SI(i+1,2)= str2num(S{idxYi(A(1,1),1)}(18:22)) ;
        Calcite_dissoluted(i+1,2)= str2num(S{idxYi2(A(1,1)-1,1)}(18:28)) ;

        %
    end
    %% Make the IC-FERST input file from PhreeqcRM
    
    fid = fopen( 'Initil_Condition_C1.txt', 'wt' );
    fprintf(fid,'%f %f \n', [Component1M]');
    fclose(fid);
    fid = fopen( 'Initil_Condition_C2.txt', 'wt' );
    fprintf(fid,'%f %f \n', [Component2M]');
    fclose(fid);
    fid = fopen( 'Initil_Condition_C3.txt', 'wt' );
    fprintf(fid,'%f %f \n', [Component3M]');
    fclose(fid);
    fid = fopen( 'Initil_Condition_C4.txt', 'wt' );
    fprintf(fid,'%f %f \n', [Component4M]');
    fclose(fid);
    
    system('cp /phreeqc/Initil_Condition_C1.txt ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional');
    system('cp /phreeqc/Initil_Condition_C2.txt ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional');
    
    system('cp /phreeqc/Initil_Condition_C3.txt ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional');
    system('cp /phreeqc/Initil_Condition_C4.txt ../ICFERST/icferst-master/legacy_reservoir_prototype/tests/Porous_compositional');
    
    %% Save Data
    ICFERST_DATA(:,:,step)=Par;
    C1M_ICF(:,:,step)=Component1MassFraction;
    C2M_ICF(:,:,step)=Component2MassFraction;
    C3M_ICF(:,:,step)=Component3MassFraction;
    C4M_ICF(:,:,step)=Component4MassFraction;
    pH_ICF(:,:,step)=pH;
    Calcite_SI_ICF(:,:,step)=Calcite_SI;
    Calcite_dissoluted_ICF(:,:,step)=Calcite_dissoluted;
    
    
end
TotalC1=[];
TotalC2=[];
TotalC3=[];
TotalC4=[];

TotalC1(:,1)=C1(:,2);
TotalC2(:,1)=C2(:,2);
TotalC3(:,1)=C3(:,2);
TotalC4(:,1)=C41(:,2);
TotalpH(:,1)=pH(:,2);
TotalCalcite_SI(:,1)=Calcite_SI(:,2);
TotalCalcite_dissoluted(:,1)=Calcite_dissoluted(:,2);

for i=1:step
    TotalC1(:,i+1)=C1M_ICF(:,1,i);
    TotalC2(:,i+1)=C2M_ICF(:,1,i);
    TotalC3(:,i+1)=C3M_ICF(:,1,i);
    TotalC4(:,i+1)=C4M_ICF(:,1,i);
    TotalpH(:,i+1)=pH_ICF(:,1,i);
    TotalCalcite_SI(:,i+1)=Calcite_SI_ICF(:,1,i);
    TotalCalcite_dissoluted(:,i+1)=Calcite_dissoluted_ICF(:,1,i);
end
p(:,1)=(MeshData{1,2}(:,1));
p(:,2)=(MeshData{1,3}(:,1));
t=delaunayn(p);

for j=1:step
    writeVTK(C_1(j),t,p,TotalC1(:,j));
    writeVTK(C_2(j),t,p,TotalC2(:,j));
    writeVTK(C_3(j),t,p,TotalC3(:,j));
    writeVTK(C_4(j),t,p,TotalC4(:,j));
end

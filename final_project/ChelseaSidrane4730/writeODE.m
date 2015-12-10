function writeODE(A,b,flag)
% extract parameters
if flag %if DAE
    n=length(b)/5;
else
    n=length(b);
end
% create file
fID = fopen('nPendODE.m','w');
fprintf(fID,'function zdot = nPendODE(t,z,p)\n');
fprintf(fID,'%% Unpack parameters\n');
fprintf(fID,'n=p.n; g=p.g;\n');
% create and print parameter and state assignments
substr = {'m','Ig','l','d'};
str = '';
for j=1:length(substr)
    for i=1:n
        strNew = sprintf('%s%d = p.%s(%d); ',substr{j},i,substr{j},i);
        str = [str,strNew];
    end
    str = [str, '\n'];
end
fprintf(fID,[str,'\n']); % print to file and add newline
%
fprintf(fID,'%% Unpack state\n');
str = '';
for i=1:n
    strNew = sprintf('t%d = z(%d); ',i,i);
    str = [str,strNew];
    strNew = sprintf('td%d = z(n+%d); ',i,i);
    str = [str,strNew];
end
str = [str,'\n'];
fprintf(fID,[str,'\n']); % print to file and add newline
% print the vector of unknowns that was symbolically solved for
fprintf(fID,'%% Assign matrices\n');
if flag % if DAE
    matSize = length(b);
else
    matSize = n;
end
str = '';
for i=1:matSize
    for j=1:matSize
        expr = char(A(i,j));
        strNew = sprintf('A(%d,%d) = %s;     ',i,j,expr);
        str = [str,strNew];
    end
    str = [str,'\n'];
end
fprintf(fID,[str,'\n']); % print to file and add newline
%
str = '';
for i=1:matSize
    expr = char(b(i));
    strNew = sprintf('b(%d,1) = %s;\n',i,expr);
    str = [str,strNew];
end
fprintf(fID,[str,'\n']); % print to file and add newline
% solve system of equations
fprintf(fID,'%% Solve system of equations\n');
if flag==0
    fprintf(fID,'tdd = A\\b; \n');
elseif flag ==1 % we are dealing with the DAE case
    fprintf(fID,'uks = A\\b; \n');
    fprintf(fID,'tdd = uks(2*n+1:2*n+n); \n');
end

fprintf(fID,'\n'); % add newline
% assign change of state vector
fprintf(fID,'%% Assign change of state\n');
fprintf(fID,'zdot = [z(n+1:2*n); tdd];\n');

fclose(fID);

n = 15; %nb des noeuds
d = 300; %nb des itérations
c = 1/(n+1);
a = 200; %longueur de l'espace
b = 150; %largeur de l'espace
r_det = (((a.^2)+(b.^2)).^(1/2))/3; %rayon de détection
r_sec = r_det/3; %rayon de sécurité

x1 = randsample(double(0:5:a),n); %double(0:5:a) to generate a list of values from 0 to a seperated by 5
                                 %randsample is for taking values randomly
                                 %from 0 to n with 5 seperateds
y1 = randsample(double(0:5:b),n);
%dist = ((x1(1)-x1(5)).^2 + (y1(1)-y1(5)).^2).^(1/2);
M = positionToAdjMatrix(x1,y1,r_det);

%M = reshape(M,[n,n]);
P = pMatrix(M,c,n);

lap = Laplacienne(M,n);

W = randStocastic(n)

C = adjacentToCom(M,n)

A = matCom(M,W);

%show_nodes(M,W,x1,y1,d,n)
show_nodes_sec(C,x1,y1,r_sec,d,n)

%show_graph(M,n,d)

function mat = positionToAdjMatrix(x,y,r_det)
    n = length(x);
    mat = zeros(n);
    
    for i=1:n
        for j=1:n
            
            dist_ij = ((x(i)-x(j)).^2 + (y(i)-y(j)).^2).^(1/2);
            if dist_ij <= r_det
                mat(i,j) = 1;
            end
        end
    end
    
    
    
    for i=1:n
        mat(i,i) = 0;
    end

    
end

function L = Laplacienne(M,n)
    D = eye(length(M));
    deg = sum(M,2);
    for i=1:n
       D(i,i) = deg(i);
    end
    
    L = D - M;
end

function C = adjacentToCom(M,n)
    C(1:n,1:n) = M;
    fact = sum(M,2);
    %M = double(M);
    for i=1:n
        if fact(i) == 0
            fact(i) = 1;
        end
        C(i) = C(i)/fact(i);
    end
end

function P = pMatrix(M,c,n)
    I = eye(length(M));
    L = Laplacienne(M,n);
    P = I - c.*L;
end

function m = randStocastic(n)
    matrix = rand (n,n)
    sumo  = sum(matrix,2);
    for i=1:n
        matrix(i,:) = matrix(i,:)/sumo(i,:);
    end
    m = matrix;
end

function A = matCom(M,W)

    I = eye(length(M));
    J = ones(length(M));
    Q = M + I;
    Qbar = J - Q;
    %C = dot(W.*Qbar,J);
    A = Q.*W + ((W.*Qbar)*J).*I;
end

function [x,y]  = mouvement(x,y,A,n)
    for i=1:n 
        for j=1:n
            x(i)= x(i)+(x(j)-x(i))*A(i,j);
            y(i)= y(i)+(y(j)-y(i))*A(i,j);
        end
    end
end

function [x,y] = mouvement_sec(x,y,r_sec,C,n)
    for i=1:n
        xc = x(i);
        yc = y(i);
        
        for j=1:n
            xc = xc + (x(j)-x(i))*C(i,j);
            yc = yc + (y(j)-y(i))*C(i,j);
            
            l= ((xc - x(i)).^2 + (yc - y(i)).^2).^(1/2);
            
            if l > 2*r_sec
                x(i) = x(i) + (xc - x(i))*(r_sec/l);
                y(i) = y(i) + (yc - y(i))*(r_sec/l);

            end
        end
    end
end

function show_graph(adj,n)
    %rows = zeros(1,n*n);
    nodes = {};
    for i=1:n
        nodes{i} = int2str(i);
    end
    G = graph(adj,nodes);
    plot(G);
    %for i=1:d
        %G = graph(adj,nodes);
        %plot(G);
        %pause(1)
    %end
    %G = digraph(rows,cols,weight,nodes); 
    %plot(G);
end

function show_nodes(M,W,x,y,d,n)
    for i=1:d
       A = matCom(M,W);
       show_graph(M,n)
       scatter(x,y,'filled')
       for i=1:n
            text (x(i),y(i),int2str(i))
       end
       [x,y] = mouvement(x,y,A,n);
       %x
       %y
       pause(0.2)
    end
end

function show_nodes_sec(C,x,y,r_sec,d,n)
    for i=1:d
       scatter(x,y,'filled')
       for i=1:n
            text (x(i),y(i),int2str(i))
       end
       [x,y] = mouvement_sec(x,y,r_sec,C,n);
       %x
       %y
       pause(0.02)
    end
end



%show_graph(M,n);
 
 %cols = find(M==1);
    %cols;
    %%%%%%%%%%%
    %j=0;
    %rows=0;
    %for i=1:n
        %row_ones = find(M(i,:)==1);
        %for j=j+1:j+length(row_ones)
            %rows(j) = i;
        %end
    %end
    %rows;
    
    %cols = num2cell(cols);
    %for i=1:length(cols)
        %cols{i} = int2str(cols{i});    
    %end 
    
    %rows = num2cell(rows);
    %for i=1:length(rows)
        %rows{i} = int2str(rows{i});    
    %end   
    
    %weight = randi(100,size(n));
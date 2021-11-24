n = 15; %nombre des noeuds
d = 40; %nombre d'itérations
c = 1/(n+1); %constante
a = 200; %longueur du plan
b = 150; %largeur du plan
r_det = (((a.^2)+(b.^2)).^(1/2))/3; %rayon de détection
r_sec = r_det/3; %rayon de sécurité

x_agent = randsample(double(0:5:a),n); %double(0:5:a) pour generer une liste de valeur de "0" à "a" avec un espace de 5 entre les valeurs
                                 %randsample pour choisir "n" valeurs 
                                 %de 0 a n avec 5 séparés
y_agent = randsample(double(0:5:b),n);

M = positionToAdjMatrix(x_agent,y_agent,r_det);

%M = reshape(M,[n,n]);
P = pMatrix(M,c,n);

lap = Laplacienne(M,n);

W = randStocastic(n)

C = adjacentToCom(M,n)

A = matCom(M,W);

%show_nodes(M,W,x_agent,y_agent,d,n) %en respectant la distence de
%detection entre agents
show_nodes_sec(C,x_agent,y_agent,r_sec,d,n)%affichage du graphe en 

%show_graph(M,n,d)

function mat = positionToAdjMatrix(x,y,r_det) %matrice d'adjacence qui contient :
                                              % "1" sur (i,j) et (j,i) si
                                              % l'agent i est proche d'une distance
                                              % inferieur au rayonde detection de l'agent j
    n = length(x); % nombre d'agent
    mat = zeros(n); %matrice de zero
    
    for i=1:n
        for j=1:n
            
            dist_ij = ((x(i)-x(j)).^2 + (y(i)-y(j)).^2).^(1/2); %calcul de distence entre les agents
            if dist_ij <= r_det
                mat(i,j) = 1; %si distance entre les deux agents est inférieur au rayon de sécurité ajouter 1 en (i,j) et (j,i)
            end
        end
    end
    
    
    
    for i=1:n
        mat(i,i) = 0;%0 sinon
    end

    
end

function Laplacienne = Laplacienne(M,n)
    D = eye(length(M)); %Matrice d'identité
    deg = sum(M,2);
    for i=1:n
       D(i,i) = deg(i); %matrice de degrée (chaque (i,j) contient le nombre de liens entre le neoud i et j)
    end
    
    Laplacienne = D - M; %matrice laplacienne est la difference de la matrice de degree et la matrice d'adjacence
end

function C = adjacentToCom(M,n) %extraire la matrice de communication à partir de la matrice d'adjacence
    C(1:n,1:n) = M; %faire une copie de la amtrice M
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
    matrix = rand (n,n) %generer une matrice de taille nxnaleatoire de valeur entre 0 et 1
    somme_facteur  = sum(matrix,2);
    for i=1:n
        matrix(i,:) = matrix(i,:)/somme_facteur(i,:);%chaque ligne de la matrice est divisé par la somme des elements de la ligne
    end
    m = matrix;
end

function A = matCom(M,W) %matrice de communication à partir de la matrice d'adjacence et la matrice de priorité

    I = eye(length(M)); %matrice d'identite
    J = ones(length(M));%matrice qui contient des "1" partout
    Q = M + I;
    Qbar = J - Q;
    %C = dot(W.*Qbar,J);
    A = Q.*W + ((W.*Qbar)*J).*I; %matrice de communication
end

function [x,y]  = mouvement(x,y,A,n) %couple de mouvements d'un agent
    for i=1:n 
        for j=1:n
            x(i)= x(i)+(x(j)-x(i))*A(i,j); %deplacement elementaire sur l'axe des abscisses
            y(i)= y(i)+(y(j)-y(i))*A(i,j); %deplacement elementaire sur l'axe des cordonnées
        end
    end
end

function [x,y] = mouvement_sec(x,y,r_sec,C,n) %retourne les coordonées de mouvement d'un agent en respectant le rayon de sécurité
    for i=1:n
        xc = x(i);
        yc = y(i);
        
        for j=1:n
            xc = xc + (x(j)-x(i))*C(i,j);
            yc = yc + (y(j)-y(i))*C(i,j);
            
            l = ((xc - x(i)).^2 + (yc - y(i)).^2).^(1/2);
            
            if l > (1.5)*r_sec %on fait un mouvement d'agent si la distance entre deux agents est superieur a 1.5*(rayon de sécurité)
                x(i) = x(i) + (xc - x(i))*(r_sec/l);
                y(i) = y(i) + (yc - y(i))*(r_sec/l);

            end
        end
    end
end

function show_graph(adj,n)%affichage du graphe en utilisant la matrice d'adjacence
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

function show_nodes(M,W,x,y,d,n) %affichage du graphique qui contient l'animation des noeuds 
    for i=1:d
       A = matCom(M,W);
       show_graph(M,n)
       scatter(x,y,'filled')
       for i=1:n
            text(x(i),y(i),int2str(i))
       end
       [x,y] = mouvement(x,y,A,n);
       %x
       %y
       pause(0.2) % 0.2secondes entre chqaue iteration
    end
end

function show_nodes_sec(C,x,y,r_sec,d,n) %affichage du graphique qui contient l'animation des noeuds en utilisant le rayon de sécurité
    for i=1:d
        %cmap = colormap(jet(size(C,2)));
        %cmap = cmap(randperm(length(cmap)),:)
        %Set colororder and plot
        %ax = axes('colororder',cmap)
        %plot(x,y)

       scatter(x,y,'filled','MarkerFaceColor',[0 .2 .7])
       for i=1:n
           circle(x(i),y(i),3);
       end
       for i=1:n
            text (x(i),y(i),int2str(i))
       end
       [x,y] = mouvement_sec(x,y,r_sec,C,n);
       %x
       %y
       pause(0.2)
    end
end

function h = circle(x,y,r)%fonction pour dessiner un cercle avec une couleur specifié
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit, 'r');
hold off
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
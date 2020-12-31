function [] = mesh(ref_area, subspace)

    % import global parameters needed for domain construction
    global vertices;
    global Ndof;
    global boundaries;
    global inputs;

    Domain.InputVertex = vertices;


    % ---------------------------------------------
    % Definizione del dominio a partire dai Vertici
    % ---------------------------------------------

    % Dichiaro le variabili per delimitare il dominio
    Domain.Boundary.Values = 1:4;
    % lato di bordo 1 dal nodo 1 al nodo 2
    % lato di bordo 2 dal nodo 2 al nodo 3
    % lato di bordo 3 dal nodo 3 al nodo 4
    % lato di bordo 4 dal nodo 4 al nodo 1

    Domain.Holes.Hole = [];       % non ci sono buchi nel dominio
    Domain.Segments.Segment = []; % non ci sono lati forzati nel dominio

    % --------------------------------------------------
    % Definizione delle condizioni al contorno a partire
    % dai Vertici e dai lati di bordo
    % --------------------------------------------------

    % numerical (constant) values of BCs (useless since BCs are set in expand function)
    BC.Values = [0.0 12.0 0.0 14.0 0.0 16.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];

    % assign BCs markers to borders
    BC.Boundary.Values = boundaries;
    % assign BCs markers to (input) vertices
    BC.InputVertexValues = inputs;

    BC.Holes.Hole = [];
    BC.Segments.Segment = [];



    % --------------------------------------------
    % Inserimento dei parametri di triangolazione
    % --------------------------------------------

    RefiningOptions.CheckArea  = 'Y';
    RefiningOptions.CheckAngle = 'N';
    RefiningOptions.AreaValue  = ref_area;
    RefiningOptions.AngleValue = [];
    RefiningOptions.Subregions = [];    

    % Construct triangulation

    global geom;
    [geom] = bbtr30(Domain,BC,RefiningOptions);

    % --------------------------------------------------
    % --------------------------------------------------

    geom.elements.coordinates = geom.elements.coordinates(...
                    1:geom.nelements.nVertexes,:);
    geom.elements.triangles = geom.elements.triangles(...
                    1:geom.nelements.nTriangles,:);
    geom.elements.borders = geom.elements.borders(...
                    1:geom.nelements.nBorders,:);
    geom.elements.neighbourhood = geom.elements.neighbourhood(...
                    1:geom.nelements.nTriangles,:);

    % --------------------------------------------------

    j  = 1;
    Dj = 1;
    for i=1:size(geom.pivot.nodelist)
         if geom.pivot.nodelist(i)==0
            geom.pivot.pivot(i)=j;
            j = j+1;
         else
            geom.pivot.pivot(i)=-Dj;
            Dj = Dj + 1;
         end
    end

    % --------------------------------------------------

    geom.pivot.pivot = transpose(geom.pivot.pivot);

    % --------------------------------------------------

    % geom.pivot.Di dopo le operazioni seguenti contiene l`indice dei nodi
    % di Dirichlet e il corrispondente marker

    [X,I] = sort(geom.pivot.Di(:,1));
    geom.pivot.Di = geom.pivot.Di(I,:);
    
    if subspace == 'P2'
        
        P2;
        Ndof = 6;
        
    else
        
        Ndof = 3;
        
    end

    clear X I;
    
end
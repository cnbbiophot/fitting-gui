function J=ajusteGlobalGral_jacobiano(dimJacobiano, params, indParams, data, varGlobales, numCurvas, numPtosxCurva, varFix, funcJacobianoHandle)
                                   
% Programa general de c�lculo del jacobiano (para cualquier funci�n).
% Primero calcula las derivadas parciales de todos los datos con respecto a las variables de ajuste (cada columna corresponde a una variable de ajuste). 
% Luego distribuye las derivadas en la columna correspondiente (segun varGlobales).
% 
% funcJacobianoHandle - Nombre de la funci�n que calcula las derivadas parciales del ajuste

% Unai, 18May2015
% agv, 24-Oct-2018

J=zeros(dimJacobiano); % Inicializa el Jacobiano
colJ=1; % Columna del Jacobiano
numVarFit=numel(varGlobales); % N�mero de variables por curva
derivParciales=funcJacobianoHandle(params, indParams, data); % Calcula una matriz de derivadas parciales de (n�total de puntos, numVarFit)
indexVarFix=1;  % Indice para construir la variable esFija (para descartar del Jacobiano las derivadas parciales de las variables fijas que no se van a ajustar)

for varFit=1:numVarFit % Calcula el jacobiano, variable por variable

    esGlobal=varGlobales(varFit);
    switch esGlobal     % Construye la variable esFija
        case 0 % Si no es global, esFija es un vector con 0s (params variables) y 1s (params fijos)
            esFija = varFix(indexVarFix : indexVarFix + numCurvas-1);
            indexVarFix = indexVarFix + numCurvas;
        case 1 % Si es global, esFija es 0 (variable) o 1 (fija)
            esFija = varFix(indexVarFix);
            indexVarFix = indexVarFix + 1;
    end        
    jParcial=derivParciales(:,varFit);  % Vector de las derivadas parciales de la variable n� varFit

    switch esGlobal % Distribuye la derivada parcial, dependiendo si la variable es global o independiente, seg�n corresponda.
        
        case 0 % si es variable independiente pone las derivadas parciales de cada curva en diferentes columnas (ceros en la parte de la columna de la otra curva)
            
            if numel(find(esFija==0)) == 0 % Si todas las variables son fijas, cuntinue 
                continue 
            end 
            
            colJdesde=colJ; % Columna del par�metro X en la curva y1
            colJhasta = colJdesde + numel(find(esFija==0)) - 1;  % Columna del par�metro X en la curva yn=n�curvas
            iteracion=1;    % Iteraci�n para la posici�n en la columna
            filaJhasta=0;   % Fila que indica d�nde acaba cada curva
            
            posCurvaParam = find(esFija==0); % N� curva en cada par�metro variable
            
            for colJnoGlobal=colJdesde:colJhasta     
                if posCurvaParam(iteracion)-1 == 0 % Cada curva empieza en una fila
                    filaJdesde = 1;
                else
                    filaJdesde = sum(numPtosxCurva(1:posCurvaParam(iteracion)-1)) + 1;
                end
                filaJhasta=sum(numPtosxCurva(1:posCurvaParam(iteracion)));
                J(filaJdesde:filaJhasta,colJnoGlobal)=jParcial(filaJdesde:filaJhasta);
                iteracion=iteracion+1;
            end
            colJ=colJhasta+1;
        
        case 1 % si es variable global pone las derivadas parciales directamente en toda la columna
        
            switch esFija
                case 0
                    J(:,colJ)=jParcial;
                    colJ=colJ+1;
                case 1  % Si es fija, no entra en el Jacobiano
                    continue
            end
    end
end
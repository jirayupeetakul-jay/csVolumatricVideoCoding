function [x, iter] = AStarOMP(y,V,N,optionsStruct)

% function [x, iter] = runAStarOMP(y,V,N,optionsStruct)
%
% function to perform A*OMP for reconstruction of a sparse signal x from 
% the random measurements y = Vx where V is the observation matrix. It 
% searchs for the sparsest approximation of y using at maximum 'Kmax' 
% vectors among columns of V and outputs the reconstructed sparse signal 
% x such that Vx = y. Orthogonalization is performed using '\' (mldivide) 
% operator of Matlab.
%
% AStarOMP may be configured using the optionsStruct. If run without the 
% optionsStruct, default values are used. Note that these default values
% are only meaningful for the noise-free case. It is strongly recommended 
% to adjust the parameters using the optionsStruct if the measurements are 
% noisy. See below for explanation of the optionsStruct.
% 
% For details of A*OMP, see 
% [1] Karahanoglu and Erdogan, "A* Orthogonal Matching Pursuit: Best-First 
% Search for Compressed Sensing Signal Recovery," Digital Signal
% Processing, 2012
% [2] Karahanoglu and Erdogan, "A Comparison of termination criteria for 
% A*OMP," in EUSIPCO'2012
% 
% Columns of V should be normalized.
%
% INPUT PARAMETERS
%   y:  observed vector (Mx1 where M is the number of observations)
%   V:  dictionary (MxN)
%   N:  length of the returned sparse signal x
%   optionsStruct: struct that contains A*OMP options (may be omitted)
%
% A*OMP OPTIONS:
% Following fields can be set in the optionsStruct:
%   I:  number of initial paths in the search stack (I>=1, integer)   
%   B:  number of branches added to the search stack at each iteration 
%       (B>=1, integer)  
%   P:  number of maximum paths in the search stack (B>=1, integer)
%   alpha: parameter for multiplicative auxiliary function (0<alpha<=1)
%   beta: parameter for adaptive auxiliary function (beta>=1)
%   AuxMode: auxiliary function mode ('Adap' for adaptive, 'Mul' for 
%       multiplicative, 'Add' for additive cost functions [1], 'AdapMul' for
%       adaptive-multiplicative cost function[2])
%   Control:  enable/disable check if input variables are correct
%   Display:  enable/disable display output ('Off': only errors and 
%       warnings are displayed. 'On': Parameter values and additional 
%       messages are also displayed.)
%   Kmax:  maximum allowed sparsity level (integer, maximum number of nonzero 
%       entries in x) (see [2])
%   eps: termination threshold for the residual power (terminates when 
%       ||r||_2/||y||_2 < eps, where r is the residue of the observation
%       vector y.) (see [2])
%
%   default options:
%       I = 3;
%       B = 2;
%       P = 200;
%       alpha = 0.97;
%       beta = 1.2;
%       AuxMode = 'AdapMul';
%       Control = 'On';
%       Display = 'On';
%       Kmax = M/2; (M: observation length)
%       eps = 1e-7;
%
% OUTPUT:
%   x:  (Nx1)reconstructed K-sparse signal
%   iter: number of total A*OMP iterations
%
% TERMINATION OF THE SEARCH
%   In [1] and [2], two different termination criteria are discussed:
%   
%   First one in [1] is based on the sparsity, that is terminate when the 
%   best path has length Kmax. To use this termination criterion, set 
%   Kmax equal to the sparsity level, and eps = 0.
%   
%   In [2], we evaluate termination based on the residual power, which
%   stops when ||r||_2/||y||_2 < eps, where r is the residue of the 
%   observation vector y. To use this criterion, set Kmax > the actual
%   sparsity level and eps wrt. the noise level in the problem (and very 
%   small for the noise-free case). It is recommended to set Kmax at least
%   bit larger than the sparsity level of interest. Performance may degrade
%   if Kmax is selected very close to the actual sparsity level.
%
%   Practically, the residue-based termination criteria can be used with 
%   Kmax = M/2. It is recommended to set eps ~ 1e-7 for noise-free case.
%
%   Note: It is recommended to avoid 
%       - "Add" cost model for residue-based termination, and 
%       - "AdapMul" cost model for sparsity-based termination
%   as they may lead to suboptimal results.
%   
%   COST MODELS
%   [1] and [2] discuss different cost models for A*OMP. We recommend 
%   AdapMul cost model with residue-based termination, as this combination 
%   returns returns the fastest and most accurate results for most cases. 
%   (see [2]).
% 
% AStarOMP.m version 1.1
% Nazim Burak Karahanoglu, Sabancý University, 09.2012
%   email: karahanoglu at sabanciuniv.edu

%-------------------------------------------------------------------------%

%% default parameters:
I = 3;
B = 2;
P = 200;
alpha = 0.97;
beta = 1.2;
AuxMode = 'AdapMul';
AuxFunc = inline('pathErr.* (AuxFuncParam*pathErr/oldErr)^(Kmax-PL);','pathErr','AuxFuncParam', 'Kmax', 'PL','oldErr');
AuxFuncParam = alpha;
NC = 0;
ND = 0;
M = length(y);
Kmax = M/2;
eps = 1e-7;
noVectors = size(V,2);

%%  Parameters in optionsStruct 
if exist('optionsStruct')
    if( isfield(optionsStruct,'Control') )
        if strcmp(optionsStruct.Control,'Off')
            NC = 1;
        end
    end
    
    if( isfield(optionsStruct,'Display') )
        if strcmp(optionsStruct.Display,'Off')
            ND = 1;
        end
    end
    
    if isfield(optionsStruct,'I')
        I = optionsStruct.I;
    else
        if ~ND
            display('optionsStruct.I (# initial paths) not specified.');
        end
    end;
    
    if isfield(optionsStruct,'B')
        B = optionsStruct.B;
    else
        if ~ND
            display('optionsStruct.B (# branches per iteration) not specified.');
        end
    end;
    
    if isfield(optionsStruct,'P')
        P = optionsStruct.P;
        
    else
        if ~ND
            display('optionsStruct.P (# maximum search paths in stack) not specified.');
        end
    end;
    
    if isfield(optionsStruct,'AuxMode')
        AuxMode = optionsStruct.AuxMode;
    else
        if ~ND
            display('optionsStruct.AuxMode (auxiliary function mode) not specified.');
        end
    end
    
    % Auxiliary function mode
    % multiplicative auxiliary function
    if strcmp(AuxMode,'Mul')
        if isfield(optionsStruct,'alpha')
            alpha  = optionsStruct.alpha;
        else
            if ~ND
                display('AuxMode:Mul, optionsStruct.alpha not specified for AuxMode:Mul.');
            end
        end
        AuxFuncParam = alpha;
        AuxFunc = inline('pathErr.* AuxFuncParam^(Kmax-PL);','pathErr','AuxFuncParam', 'Kmax', 'PL','dummy');
    end
    
    % adaptive auxiliary function
    if strcmp(AuxMode,'Adap')
        if isfield(optionsStruct,'beta')
            beta  = optionsStruct.beta;
        else
            if ~ND
                display('AuxMode:Adap,optionsStruct.beta not specified for AuxMode:Adap.');
            end
        end
        AuxFuncParam = beta;
        AuxFunc = inline('pathErr - AuxFuncParam*(oldErr-pathErr)*(Kmax-PL);','pathErr','AuxFuncParam', 'Kmax', 'PL','oldErr');
    end
    
    % additive auxiliary function
    if strcmp(AuxMode,'Add')
        if isfield(optionsStruct,'beta')
            beta  = optionsStruct.beta;
        else
            if ~ND
                display('optionsStruct.beta not specified for AuxMode:Add.');
            end
        end
        AuxFuncParam = beta*norm(y)/Kmax;
        AuxFunc = inline('pathErr - AuxFuncParam*(Kmax-PL);','pathErr','AuxFuncParam','Kmax', 'PL','oldErr');
    end
    
        % adaptive-multiplicative auxiliary function
    if strcmp(AuxMode,'AdapMul')
        if isfield(optionsStruct,'alpha')
            alpha  = optionsStruct.alpha;
        else
            if ~ND
                display('AuxMode:AdapMul,optionsStruct.alpha not specified for AuxMode:Adap.');
            end
        end
        AuxFuncParam = alpha;  
        AuxFunc = inline('pathErr.* (AuxFuncParam*pathErr/oldErr)^(Kmax-PL);','pathErr','AuxFuncParam', 'Kmax', 'PL','oldErr');     
    end
end

%% check input parameters 
if ~NC
    if(size(V,1) ~= M)
        error('Length of y does not match length of vectors in V.');
    end
    
    if( M >= noVectors )
        warning('More dimensions than the number of vectors in V!')
    end
    
    if( M >= N )
        error('More observations than # dimensions of x.')
    end

    if rem(P,1) || P<1
        error('P (# maximum paths) should be an integer greater than 1.');
    end
    if rem(B,1) || B<=1
        error('B (# branches per iteration) should be an integer <= 1.');
    end
    if rem(I,1) || I<=1
        error('I (# initial paths) should be an integer <= 1.');
    end
    if rem(Kmax,1) || Kmax<=1
        error('Kmax (sparsity) should be an integer <= 1.');
    end
%     if Kmax >= M/2
%         warning('Kmax seems to be too large, reconstruction may be poor!');
%     end
    if strcmp(AuxMode,'Mul')
        if  alpha <=0 || alpha > 1
            error('alpha should be in the interval (0,1]');
        end
    end
    if ismember(AuxMode, {'Add' 'Adap'})
        if  beta <=1
            error('beta should be greater than or equal to 1');
        end
    end
    if P<100
        warning('# maximum paths seems too low.');
    end
    
    % Random Check if V is normalised 
    randind = randi(size(V,2),1);
    nV=norm(V(:,randind));
    if abs(nV)-1 > 0.0001;
        warning('Columns of V (holographic basis) appear not normalized.');
    end
    if ~ismember(AuxMode, {'Add' 'Adap' 'Mul'})
        error('optionsStruct.AuxMode: Invalid Auxiliary Function.')
    end
end

% display parameters
if ~ND
    display(['Parameters: Kmax = ' num2str(Kmax) ',I = ' num2str(I) ',B = ' num2str(B) ',P = ' num2str(P)]);
    if strcmp(AuxMode,'Mul')
        display(['Cost function: Multiplicative, alpha = ' num2str(alpha)]);
    end
    if strcmp(AuxMode,'Adap')
        display(['Cost function: Adaptive, beta = ' num2str(beta)]);
    end
    if strcmp(AuxMode,'Add')
        display(['Cost function: Additive, beta = ' num2str(beta)]);
    end
end

%% Initialize A*OMP
normY = norm(y);

%initialize for orthogonalization
Y = V'*y;
G = V'*V;

PathLength = zeros(P,1);
PathLength(1:I) = 1;
CostVec(1:P) = normY;
ErrVec(1:P) = normY;

SearchMatrix = zeros(P,Kmax);
CoefMat = zeros(P,Kmax);
ResMat = repmat(y,1,P);


[Coefs, ind] = sort(abs(Y),'descend');
SearchMatrix(1:I,1) = ind(1:I);
CoefMat(1:I,1) = Y(ind(1:I));

for k=1:I
    ResMat(:,k) = y - V(:,SearchMatrix(k,1))*CoefMat(k,1);
    ErrVec(k) = norm(ResMat(:,k));
    CostVec(k) = AuxFunc(ErrVec(k),AuxFuncParam, Kmax, 1,normY);
end;

% initialize selected path
[val, SP] = min(CostVec(1:I));
l = PathLength(SP);
iter = 0;
%% run A* search
while l < Kmax
    iter = iter+1;
    % copy original path
    oldPath = SearchMatrix(SP,1:l);
    oldErr = ErrVec(SP);
    
    % compute expansions
    VV = V;
    VV(:,oldPath) = 0;
    CorrXBasis = (VV)'*ResMat(:,SP);
    [Coef, ExpansionList] = sort(abs(CorrXBasis), 'descend');
    
    PathLength(SP) = 0;
    CostVec(SP) = normY;
        
    worstCost = normY;
    worstInd = SP;
    newPL = l+1;
    ind = logical(PathLength >= newPL );
    SortedSM = sort(SearchMatrix(ind,1:newPL),2);
    
    % add new candidates
    for k=1:B
        newPath = [oldPath ExpansionList(k)];
        
        sortedPath = sort(newPath);
        isToBeIgnored = 0;
        for i = 1:size(SortedSM,1)
            if(isequal(SortedSM(i,:),sortedPath))
                isToBeIgnored = 1;
                break;
            end;
        end;
    
        if not(isToBeIgnored)
            %orthogonalize the selected path   
            Coefs = G(newPath,newPath)\Y(newPath);
            pathRes = y - V(:,newPath)*Coefs;
            pathErr = norm(pathRes);
            pathCost = AuxFunc(pathErr,AuxFuncParam, Kmax, newPL, oldErr);
            
            if (pathCost < worstCost )
                
                SearchMatrix(worstInd,:) = 0;
                CoefMat(worstInd,:) = 0;
                
                SearchMatrix(worstInd,1:newPL) = newPath;
                CoefMat(worstInd,1:newPL) = Coefs;
                ResMat(:,worstInd) = pathRes;
                PathLength(worstInd) = newPL;
                CostVec(worstInd) = pathCost;
                ErrVec(worstInd) = pathErr;
                [worstCost, worstInd] = max(CostVec);
            end;
        end;
    end;
    
    % Given current paths, find the most promising one
    [val, SP] = min(CostVec);
    l = PathLength(SP);
    if norm(ResMat(:,SP))/normY < eps
        break
    end
 
end

%% reconstruct x using the selected basis vectors and coefficients
BestBasisList = SearchMatrix(SP,1:l);
BestBasisCoefs = CoefMat(SP,1:l);
BestBasisCoefs(abs(BestBasisCoefs)<1e-10) = 0;
x = zeros(N,1);
x(BestBasisList) = BestBasisCoefs;
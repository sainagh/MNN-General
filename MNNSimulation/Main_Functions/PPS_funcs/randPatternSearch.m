%{
This function performs a randomized patten search 
ALL VARIABLES ARE STORED AS COLUMN VECTORs
funPack: is a cell array that contains the following functions
    {1} evaluates error at x 
        {error, posOut, debugOut } = funPack{1}(x)
    {2} Test function
        returns true if new point is acceptable
         funPack{2}(errorLast, errorNew)
    {3} pauses evaluation to cool components
          funPack{3}
    {4} plots optimal condition graphs
          funPack{4}(s1, s2, s3, s4, s5, s6)
    {5} plots error for all attempted points
          funPack{5}(s1, s2, s3, s4, s5, s6)
nOutputs: the number of outputs from the 

%}
function [xOpt, eOpt, bestPts, allPts] =  randPatternSearch( funPack, stepSize, stepFun, minStep, loopsToStep, deMin, allIter, typNoise, outlierMult, savePath, xRange, xStart)
% DO THIS OUT OF FUNCTION
% %Get resting position
%     %Turn off array
MAD = @(x) median(abs(x -median(x)));
Mi = @(x) abs(0.6745*(x - median(x))/MAD(x));
%% Start at some intial condition
%If no intitial condition is given start in the middle of the range
if isempty('xStart') 
    xStart = ones(size(xRange,1),1)*sum(xRange,2)/2;
end
%Check error for each IC
for icNum = 1:size(xStart,2)
    %Evaluate the error at each initial condition
    outputCell = funPack{1}(xStart(:,icNum));
    
    %Store error for the best IC
    icError(icNum) = outputCell{1};
    icPos(:,icNum) = reshape(outputCell{2},1,[]); %THIS FORLDS along dim1, dim2, dim3
    icDebug{icNum} = outputCell{3};
    %Plot if required
  
end
funPack{4}(icError, icPos, [], [],[],[])
funPack{3}; %pause for cool down
%Pick the best error as a start point
%Store the best values into history vectors 
icNum = find(min(icError) == icError); % This holds eval number with the best Initial condition
xLast  = xStart(:,icNum); %Stores the order list of ics run from (For backStepping) 
eHistB = icError(icNum);
xHistB = xStart(:,icNum);
xRun = xHistB; %This x is the current eval pt
tHistB = now;

pHistB = icPos(:, icNum);
debugHistB = {icDebug{icNum}};
estList = icError(icNum);

%Initialize ALL history array's with the BEST history ARRAY
eHistA = repmat(eHistB,3,1); %Stack 3 coppies for: 0, -stepSize, and +stepSize
xHistA = repmat(xHistB,3,1);
pHistA = repmat(pHistB,3,1);
debugHistA = repmat(debugHistB,3,1);
currentEstList = icError(icNum);

stepIter = 2;
allIter = 2; %initialize counter for number of eval
iter = 2; %counter for number in '_'Bh vectors
allSteps = 1;
backStepIters = 2;
estCell = {eHistB};
stepIterHist  = inf(1000,1000);
stepIterHist(1,1) = stepIter;
%% Perform Full Optimization from chosen IC
%Until step size is bellow min size
nVar = length(xRun);
while stepSize >= minStep
    loopsUnchanged = 0;
    %UNTIL nothing is changed for a few loops repeat checks at this step size
    while loopsUnchanged < loopsToStep
        loopsUnchanged = loopsUnchanged + 1;
        disp(['running ',num2str(loopsUnchanged),' loop of ', num2str(loopsToStep)])
        %FOR Randomly step through each link
        for beam = randperm( size(xRun,1))
            %repmat current state by 3 ( high;low; stay)    
            xToTry = repmat(xRun, 3, 1); 
            pattern = [0, stepSize, -stepSize];
            xToTry((0:1:2)*nVar + beam) = xToTry((0:1:2)*nVar + beam)+ pattern';
            %Threshold
            repRange = repmat(xRange, 3,1);
            try
            xToTry(xToTry > repRange(:,2)) = repRange(xToTry > repRange(:,2),2);
            xToTry(xToTry < repRange(:,1)) = repRange(xToTry < repRange(:,1),1);
            catch
               disp('threshold error') 
               
            end
            %Initialize HOLDER Matrices
            timeEval = inf(3,1);
            errorEval = timeEval;
            posEval = inf(size(icPos,1)*3,1);
            
            %FOR each (high low and current)
            for  pat = 1:length(pattern)
                patInd =(1:nVar)+(pat - 1)*nVar;
                posInd = (1:size(icPos,1)) + (pat -1)*size(icPos,1);
                %Log current time;
                timeEval(pat) = now;
                %get the error
                outputCell = funPack{1}(xToTry(patInd));
                
                errorEval(pat) = outputCell{1};
                debugOut{pat} = outputCell{3};
                %Store the position AS STACKED COLUMN VECTOR       
                posTemp = outputCell{2};
                posEval(posInd) = reshape(posTemp,1,[]);
                %cool off
                funPack{3};
            end %- end of checking each position
           
            %Build a vector containing all evals at current pt
            currentEst = 0;
            invalidEst = false;
            estList = [estList; errorEval(1)];
%             %Until this is a good read
            while currentEst == 0
                %Check new estimate against last best value
                if length(estList) < 3
                    %for low numbers of reads
                    if max(abs(diff(estList))) > typNoise
                        invalidEst = true;
                    else
                        currentEst = mean(estList) - deMin;
                        invalidEst = false;
                    end
                else
                    %for high numbers of reads
                    mseMi = Mi(estList);
                    validMse = (mseMi < outlierMult) | isnan(mseMi);
                    if sum(validMse) < 4
                        invalidEst = true;
                    else
                        %remove outliers
                        noOutliers = nonzeros((mseMi < outlierMult).*estList);
                        currentEst = mean(noOutliers) - deMin ;
                        if length(noOutliers) < 3
                            invalidEst = true;
                        else
                            invalidEst = false; 
                        end
                    end
                end %- end of outlier checking.

                %Retake Readings
                if invalidEst
                    recheckCell = funPack{1}(xRun);
                    estList = [estList; recheckCell{1}];
                    
                end
            end %- end of retesting current position
%             %Back Stepping Code
%             %IF this settled value isn't as good as it seemed
%             currentEst + deMin
%              badTest = (((currentEst + deMin) - estList(1) )>  + typNoise )&& iter > 2;
        
%               badTest = (((currentEst + deMin) - estList(1) )> deMin )&& iter > 2;
%             currentEst = estList(1)  - deMin
%              badTest = false;
%               badTest =  (mean(estList(2:end)) - estList(1)) > typNoise;
              badTest = ((currentEst + deMin) - estList(1) )>  typNoise;
%             badTest = (currentEst - estList(1)) > deMin
%             stepIter            
% estList
            if badTest
%                 iter
                if iter > 2
                    try
                        %Use the last position as your new position
                        iter = iter -1; %Overwrite the best Iter

                        xRun = xLast(:, iter); %Go to the last position
                        estList = estCell{iter - 1} %Use the estimation from the last position
        %                 stepIter = stepIter(1:bestCounter);

                        loopsUnchanged  = 0; %Flag to repeate the loops at this size
                        backStepIters= [backStepIters, allIter]; %Log this backstep
                        allSteps = [allSteps, allIter];
        %                 %Remove this from the trajectory
                        tHistB = tHistB(1:iter);
                        eHistB = eHistB(1:iter);
                        xHistB = xHistB(:,1:iter);
                        pHistB = pHistB(:,1:iter);
                        xLast = xLast(:,1:iter);
                        stepIter =stepIter(1:iter);
                        debugHistB = debugHistB(1:iter);
                    catch
                       stepIter;
                    end
                else
                    %don't change the eval position just update the
                    %comparison pt
                    newEst = currentEst + deMin;
                    estList(iter) = newEst;
                    eHistB(iter) = newEst;
                    pHistB(iter) = posEval(1);
                end
                    

                    %Store the data in the all data Vectors
                    tHistA(:,allIter) = timeEval;
                    xHistA(:,allIter) = xToTry;
                    eHistA(:,allIter) = errorEval;
                    pHistA(:,allIter) = posEval;
                    debugHistA{allIter} = debugOut;
                    stepIterHist(1:length(stepIter), allIter) = stepIter';
                    %Plot if needed
                    funPack{5}(allIter, stepIter, eHistA', xHistA', currentEstList, deMin);
                    %Save recorded Data
                    save(savePath);
                    allIter = allIter + 1;
    %                 allIter
    %                 stepIter
                    disp(["errarAtThisPoint is: ", num2str(currentEst + deMin)])
                
                %Break search at this xToTry
                break
            end
             %IF the best pt is not where you are go to the best pt
%              currentEst = estList(1);
           [~, minEval] = min([estList(1) - deMin; errorEval(2:end)]);
%            disp("These are the comparison Values")
%            [estList(1)- deMin; errorEval(2:end)]
%            disp("Last Value was: ")
%            currentEst
             if minEval ~= 1
                %Retest the values
                
                %Reset the number of times you have gone through the
                %loop                
                loopsUnchanged = 0;
                
                xLast = [xLast, xRun];
                %Store the data in optimization path
                patInd = [1:nVar]+(minEval - 1)*nVar;
                posInd = [1:size(icPos,1)] + (minEval -1)*size(icPos,1);
                xRun = xToTry( patInd );
                tHistB(:,iter) = timeEval( minEval);
                xHistB(:,iter) = xToTry( patInd );
                eHistB(:,iter) = errorEval(minEval);
               % posEval(posInd)
                pHistB(:,iter) = posEval(posInd);
                debugHistB{iter} = debugOut{minEval};
                %Update the estimation list for the debounce   
                estCell{iter -1} = estList;
                estList = errorEval(minEval); %Now you compare to the best pt
                allSteps = [allSteps, allIter];
                %Update the iteration counter
                stepIter(iter) = allIter
                stepIterHist(1:length(stepIter), allIter) = stepIter';
                iter = iter + 1
                
                %Plot the optimization path
                 funPack{4}(stepIter,stepIter,eHistB', xHistB', tHistB', pHistB')
            end
            %Print the data
            %Store all of the data in to the full storage path
            tHistA(:,allIter) = timeEval;
            xHistA(:,allIter) = xToTry;
            eHistA(:,allIter) = errorEval;
            pHistA(:,allIter) = posEval;
            debugHistA{allIter} = debugOut;
            stepIterHist(1:length(stepIter), allIter) = stepIter';

            currentEstList(allIter) = currentEst;
            %Plot if needed
            funPack{5}(allIter, stepIter, eHistA', xHistA', currentEstList, deMin);
            %Plot the optimization path
             funPack{4}(stepIter,stepIter,eHistB', xHistB', tHistB', pHistB')
            %Save recorded Data
            save(savePath);
            allIter = allIter + 1;
            
        end %- end of random beam selection
    end %- end of checking for a given step size
    %Reduce Step size using step function
    stepSize = stepFun(stepSize);
    disp(['Step size is now: ', num2str(stepSize)])
end%- end of optimization
%Pack output data and return

xOpt = xRun;
eOpt = currentEst + deMin;
bestPts.tHistB = tHistB;
bestPts.xHistB = xHistB;
bestPts.eHistB = eHistB;
bestPts.pHistB = pHistB;
bestPts.debugHistB = debugHistB;
bestPts.stepIter = stepIter;

allPts.tHistA = tHistA;
allPts.xHistA = xHistA;
allPts.eHistA = eHistA;
allPts.pHistA = pHistA;
allPts.debugHistA = debugHistA;
allPts.stepIter = stepIter;
allPts.backStepIters = backStepIters;
allPts.allSteps = allSteps;
allPts.stepIterHist = stepIterHist;
allPts.estCell =  estCell;
end


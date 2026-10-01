function outCell = errorMatlabFromDAQ_PP(x, sysPar, testPar, evalObj, dq, arrayPar, triggerPar)
 persistent testNum xHist errorHist dispHist timeHist
%if this is the first time the system has been run
%Initalize the holing matricies
%Alter history Vector State #THESE WILL BE A DIFFERENT SIZE SINCE 3X MORE
%PTS


% 3
% nReads = sysPar.nReads;


if isempty(testNum)
    %This function is going to save nPts points for every time that the
    %function is called and store them stac
    maxIters = sysPar.maxIters*3;
    %Initialize the holding matrix
    testNum = 1;
    xHist = zeros(maxIters, length(x));
    errorHist = zeros(maxIters,2);
    dispHist = zeros(maxIters, numel(testPar.zerosXY)*size(testPar.forceCases,3)); %iteration, dimension*beh
    timeHist = zeros(maxIters, length(clock));
end

%Collect the data using correct data evaluation interface
%Hardware interface

 if testPar.testMode == 0
     for ind = 1:testPar.numEvals
        disableMotors(evalObj,0)
        [mseOut(ind), maxDistOut(ind), dispOut(:,:,:,ind), voltOut] = getErrorDAQ_12Node(x, testPar, evalObj, dq, arrayPar, triggerPar);
        % only collect output nodes
%         dispOut = dispOut([4,8],:,:);
        disableMotors(evalObj,1)
        pause(.75);
        %Store the collected data in the history vectors
        xHist(testNum,:) = x';
        errorHist(testNum,:) = [mseOut, maxDistOut];
        dispHist(testNum,:) = reshape(dispOut,[],1)';
        errorOut = errorHist(testNum, sysPar.errorType);
        timeHist(testNum,:) = clock;
        testNum = testNum + 1;
     end
    outCell = {mean(mseOut), mean(dispOut,4), voltOut};
%Simulation interface
elseif testPar.testMode == 1
    [mseOut, maxDistOut, dispOut, voltOut] = getErrorSim_12Node( 1000*x, testPar, sysPar, evalObj,arrayPar, triggerPar);
    outCell = {mseOut, dispOut, voltOut};
    %Store the collected data in the history vectors
    xHist(testNum,:) = x';
    errorHist(testNum,:) = [mseOut, maxDistOut];
    dispHist(testNum,:) = reshape(dispOut,[],1)';
    errorOut = errorHist(testNum, sysPar.errorType)
    timeHist(testNum,:) = clock;
    testNum = testNum + 1;
 end



%Check to see if you need to run the "Cooling Function" 
%SAVES HISTORY VECTORS WHILE COOLING
if mod(testNum, sysPar.coolStep) == 0
        if testPar.testMode == 0
            disableMotors(evalObj,1)
        end
        
        save([sysPar.path, sysPar.file,'running.mat'] ,'testNum','xHist','errorHist','dispHist','sysPar','testPar',...
             'arrayPar','triggerPar','timeHist')
        
         if sysPar.dispText
            disp(['testNum = ',num2str(testNum)]) 
            disp(['Error = [',num2str(mseOut),' , ', num2str(maxDistOut),']']);
        end
    
        if testPar.testMode == 0
            pause(sysPar.coolTime);
            disableMotors(evalObj,0)
        end
end

%If asked for save the long data vectors into seperate files (speed up the
%larger file saves by splitting the saves up
if sysPar.debugSave == true
    save([sysPar.path,'\',num2str(testNum),'vHist','.mat'],'voltOut', 'x')
end

%Save all of the output data to file
     if sysPar.Last
         save([sysPar.path, sysPar.file,'running.mat'] ,'testNum','xHist','errorHist','dispHist','sysPar','testPar',...
             'arrayPar','triggerPar','timeHist')
     end
end
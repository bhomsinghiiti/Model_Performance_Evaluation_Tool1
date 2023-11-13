classdef btpapp1_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        ResultsdisplayPanel             matlab.ui.container.Panel
        ResultsTextArea                 matlab.ui.control.TextArea
        StatisticalPerformanceEvaluationMeasuresPanel  matlab.ui.container.Panel
        SimulatedEditField              matlab.ui.control.NumericEditField
        SimulatedEditFieldLabel         matlab.ui.control.Label
        ObservedEditField_2             matlab.ui.control.NumericEditField
        ObservedEditField_2Label        matlab.ui.control.Label
        RunButton                       matlab.ui.control.Button
        Tree                            matlab.ui.container.CheckBoxTree
        ScaledependentmeasuresNode      matlab.ui.container.TreeNode
        MeanErrorMENode                 matlab.ui.container.TreeNode
        MeanAbsoluteErrorMAENode        matlab.ui.container.TreeNode
        AbsoluteMaximumErrorAMENode     matlab.ui.container.TreeNode
        MeanSquaredErrorMSENode         matlab.ui.container.TreeNode
        RootMeanSquaredErrorRMSENode    matlab.ui.container.TreeNode
        FourthRootOfTheMeanQuadrupledErrorR4MS4ENode  matlab.ui.container.TreeNode
        MeanSquaredLogarithmicErrorMSLENode  matlab.ui.container.TreeNode
        MeanSquaredDerivativeErrorMSDENode  matlab.ui.container.TreeNode
        FuzzyMeanSquaredErrorFMSENode   matlab.ui.container.TreeNode
        PeakDifferencePDIFFNode         matlab.ui.container.TreeNode
        ScaleindependentmeasuresNode    matlab.ui.container.TreeNode
        FractionalStandardErrorFSENode  matlab.ui.container.TreeNode
        RootMeanSquaredErrortoStandardDeviationRatioRSRNode  matlab.ui.container.TreeNode
        InertiaRootMeanSquaredErrorIRMSENode  matlab.ui.container.TreeNode
        PercentagebasedmeasuresNode     matlab.ui.container.TreeNode
        MeanAbsolutePercentageErrorMAPENode  matlab.ui.container.TreeNode
        MedianAbsolutePercentageErrorMdAPENode  matlab.ui.container.TreeNode
        SymmetricMeanAbsolutePercentageErrorsMAPENode  matlab.ui.container.TreeNode
        PercentageErrorinPeakPEPNode    matlab.ui.container.TreeNode
        RunoffCoefficientErrorRRNode    matlab.ui.container.TreeNode
        RelativemeasuresNode            matlab.ui.container.TreeNode
        MeanRelativeErrorMRENode        matlab.ui.container.TreeNode
        RelativeAbsoluteErrorRAENode    matlab.ui.container.TreeNode
        RelativeVolumeErrorRVENode      matlab.ui.container.TreeNode
        MeanAbsoluteRelativeErrorMARENode  matlab.ui.container.TreeNode
        MeanSquaredRelativeErrorMSRENode  matlab.ui.container.TreeNode
        MultiplicativeBiasMBNode        matlab.ui.container.TreeNode
        Panel_2                         matlab.ui.container.Panel
        PlotButton                      matlab.ui.control.Button
        DownloadGraphButton             matlab.ui.control.Button
        UIAxes                          matlab.ui.control.UIAxes
        Panel                           matlab.ui.container.Panel
        SimulateddataEditField          matlab.ui.control.NumericEditField
        SimulateddataEditFieldLabel     matlab.ui.control.Label
        ObserveddataEditField           matlab.ui.control.NumericEditField
        ObserveddataEditFieldLabel      matlab.ui.control.Label
        TimeindexEditField              matlab.ui.control.NumericEditField
        TimeindexEditFieldLabel         matlab.ui.control.Label
        ProvidecolumnindexLabel         matlab.ui.control.Label
        Lamp                            matlab.ui.control.Lamp
        SelectplottypeDropDown          matlab.ui.control.DropDown
        SelectplottypeDropDownLabel     matlab.ui.control.Label
        UploaddataxlsxButton            matlab.ui.control.Button
    end

        properties (Access = private)
        data % Description
        check % to handle scatter
        sNodes % to handle selected nodes
        % data2 % for part 2 means right side asking for observed and simulated;
    end
    
    
    methods (Access = private)
        
     function results = ME(app, observed, simulated)
        ans=0; % var to store result 
    n=length(observed); % to get length of given vector
    for i = 1:n
        % Access the current element using indexing
        ans=ans+simulated(i);
        ans=ans-observed(i);
            % Perform operations on the current element
    end
    %to calculate ME divide the summation by n
    results=ans/n;           
   
        end
        
       function results = MAE(app,  observed, simulated)
        ans=0; % var to store ongoing calculation
        n=length(observed); % length of given data vectors
      for i = 1:n
        % Access the current element using indexing
        temp=abs(simulated(i)-observed(i)); % temp var to store absolutue value
        ans=ans+temp; % update ongoing calculation
         % Perform operations on the current element
    end
    %divide the summation by n;
    results=ans/n;
       end
   

  function results=AME(app, observed, simulated);
    ans=0; % var to store result 
    n=length(observed); % to get length of given vector
    for i = 1:n
        % Access the current element using indexing
        temp=abs(simulated(i)-observed(i));
        ans=max(ans, temp);
         % Perform operations on the current element
    end
    results=ans; % finally results=ans;
    
  end



function results=MSE(app, observed, simulated);
    ans=0;
    n=length(observed); % to get the size of vector
    for i = 1:n
        % Access the current element using indexing
        temp=abs(simulated(i)-observed(i));
        x=temp*temp;
        ans=ans+x;
         % Perform operations on the current element and doing calculations
    end
    results=ans/n; % final result acc. to formula of MSE
end
   

function results=RMSE(app, observed, simulated);
 ans=0; % to store ongoing calculation
n=length(observed); % to get the length of vector
 % iteration using for loop for calculation
for i = 1:n
    % Access the current element using indexing
    temp=abs(simulated(i)-observed(i));
    x=temp*temp;
    ans=ans+x;
     % Perform operations on the current element
end
ans=ans/n; % divide by n 
ans=sqrt(ans); % taking root acc. to formula

results=ans; % to return final res
end

 function results=R4MS4E(app, observed, simulated);
    ans=0; % to store ongoing calculations
    n=length(observed); %to get length of given vector
    % for loop for doing calculation on given length of vector
    for i = 1:n
        % Access the current element using indexing
        temp=abs(simulated(i)-observed(i));
        x=temp*temp*temp*temp;
        ans=ans+x;
         % Perform operations on the current element
    end
    ans=ans/n; % acc. to formual
    ans=power(ans, 1/4); % acc. to formula
    results=ans; % to return final answer
end
    
       
function results=MSLE(app, observed, simulated);
    ans=0; % to store ongoing calculations
    n=length(observed); % to get length of vector
    % for loop for doing calculations
    for i = 1:n
        % Access the current element using indexing
        x=log(simulated(i))-log(observed(i));
         temp=x*x;
        ans=ans+temp;
            % Perform operations on the current element
    end
    ans=ans/n; % acc. to formula

    results=ans; % to return final result
end

function results=MSDE(app, observed, simulated);
    ans=0; % to store ongoing calculations
    n=length(observed); % to get length of given vector
     % for loop for doing calculations
    for i = 2:n
        % Access the current element using indexing
        x=simulated(i)-simulated(i-1);
        y=observed(i)-observed(i-1);
        temp=y-x;
        ans=ans+temp*temp;
            % Perform operations on the current element
    end
     f=n-1; 
    ans=ans/f; % acc. to formula
    results=y; % to return final result
end


function results=PDIFF(app, observed, simulated);
    ans=0; 
    n=length(observed);
    % mx = intmin('int32'); 
    mx_simu=intmin('int32');  % to store the max value of simulated data
    mx_obs=intmin('int32');   % to store the max value of observed data
    
    for i = 1:n
        % Access the current element using indexing
        mx_simu=max(mx_simu, simulated(i));
        mx_obs=max(mx_obs, observed(i));
    end
    ans=mx_simu-mx_obs; % acc. to formula
    results=ans; % to return the final result
end


function results=FMSE(app, observed, simulated);
    ans=0; % to store ongoing calculations
    n=length(observed); % to get the length of vector
     % for loop for doing calculations
    for i = 1: n
        % Access the current element using indexing
        temp=abs(simulated(i)-observed(i));
        x=temp*temp;
        ans=ans+x;
        % Perform operations on the current element
    end
    ans=ans/n; % acc. to formula
    results=ans;
    
end


function results = FSE(app, observed, simulated)
    ans=0; % to store ongoing calcualtions
    r=RMSE(app, observed, simulated);  % calling RMSE function to use value of rmse acc. to  formula
    m=mean(observed); % to get mean of observed vector
    ans=r/m; % according to formula
    results=ans; % to return result
end

function results = RSR(app, observed, simulated)
    ans=0; % to store ongoing calcualtions
    r=RMSE(app, observed, simulated);  % calling RMSE function to use value of rmse acc. to  formula
    sd=std(observed); %computes the standard deviation of the elements in the array or vector
   
    ans=r/sd; % according to formula
    results=ans; % to return result
end


      
function results=IRMSE(app, observed, simulated);
    
   r=RMSE(app, observed, simulated); % to get RMSE 
    n=length(observed); % to get the length of vector

    x=0;     % x is used to store ongoing cal. for one part of formula
    for i = 2:n
        % Access the current element using indexing
        % x=simulated(i)-simulated(i-1);
        temp=observed(i)-observed(i-1);
        x=x+temp;
         % Perform operations on the current element
    end
    x=x/n; % acc. to formula

    den=0; % calculation for outer part of formula
    for i = 2:n
        % Access the current element using indexing
        % x=simulated(i)-simulated(i-1);
        temp=observed(i)-observed(i-1);
        temp=temp-x;
        den=den+temp*temp;
            % Perform operations on the current element
    end
    % acc. to given formula
    nn=n-1;  
    den=sqrt(den/nn);
    results=r/den; % to return final result
end
    

function results=MAPE(app, observed, simulated);
    ans=0; % to store ongoing calculations
    n=length(observed); % to get length of given vector
    % for loop for performing req. cal.
    for i = 1:n
        % Access the current element using indexing
        temp=simulated(i)-observed(i);
         temp1=observed(i);
         x=abs(temp/temp1);
        ans=ans+x;
            % Perform operations on the current element
    end
    ans=ans/n; % acc. to formula
    ans=ans*100; % acc. to formula
    results=ans; % to return final result
end
    
   
    function results=MdAPE(app, observed, simulated);
    %  creating a empty vector
     myVector = [];
     n=length(observed); % to get the length of given vector
     % x is used as variable to cal value inside the loop for each iteration
     x=0;
     for i = 1:n
        
         temp=abs(simulated(i)-observed(i));
         temp2=abs(observed(i));
         x=temp/temp2;
         x=x*100;
         % pushing value x into myvector
         myVector = [myVector x];
     end
    %  calculating meadian
     medianValue = median(myVector);
      y=medianValue;
      results=y; % to return final result

    end
     

function results=sMAPE(app, observed, simulated);
    ans=0; % to store ongoing calculations
    n=length(observed); % to get length of vector
    % for loop for doing cal.
    for i = 1:n
        % Access the current element using indexing
        temp=abs(simulated(i)-observed(i));
         temp1=abs(observed(i))+abs(simulated(i));
         x=temp/temp1;
        ans=ans+x;
            % Perform operations on the current element
    end
    ans=ans/n; % acc. to formula
    
    ans=ans*200; % acc. to formula

    results=ans; % to return final result
end

function results=PEP(app, observed, simulated);
    ans=0; % to store cal.
    mx_obs=max(observed); % to get max value of observed vector
    mx_sim=max(simulated); % to get max value of simulated vector
    ans=mx_sim-mx_obs; % diff
    ans=ans/mx_obs; % acc. to formula
    ans=ans*100; % acc. to formula
    results=ans; % to return final result
end


function results=RR(app, observed, simulated);
    ans=0; % to store ongoing cal.
    n=length(observed);
    for i = 1:n
        % Access the current element using indexing
        temp=simulated(i)-observed(i);
       ans=ans+temp;
            % Perform operations on the current element
    end
     sum_obs=sum(observed); % to get Denominator of formula || sum of observed vector
    ans=ans/sum_obs; % acc. to formula
    ans=ans*100; % acc. to formula
    results=ans; % to return final result
end



function results=MRE(app, observed, simulated);
    ans=0; % to store ongoing cal.
    n=length(observed); % to get length of given vector
    for i = 1:n
        % Access the current element using indexing
        temp=simulated(i)-observed(i);
        temp=temp/observed(i);
        ans=ans+temp;
      % Perform operations on the current element
    end
    
    ans=ans/n;
    results=ans; % to return result

end


function results=RAE(app, observed, simulated);
    n1=0; % to store Numerator cal.
    n2=0; % to store Denominator cal.
    n=length(observed); % to get the length of given vector
    m=mean(observed); % mean of observed vector
    for i = 1:n
        % Access the current element using indexing
        temp=abs(simulated(i)-observed(i));
        n1=n1+temp;
        temp2=abs(observed(i)-m);
        n2=n2+temp2;
            % Perform operations on the current element
    end
    y=n1/n2;
    results=y; % to return final result
end


function results=RVE( app, observed, simulated);
    n1=0; % to store Numerator cal.
    n2=0; % to store Denominator cal.
    n=length(observed); % to get the length of given vector
        for i = 1:n
        % Access the current element using indexing
        n1=n1+simulated(i)-observed(i);
        n2=n2+observed(i);
            % Perform operations on the current element
    end
    ans=n1/n2; % acc. to formula
    y=ans*100; % acc. to formula
    results=y; % to return final result
end


function results=MARE(app, observed, simulated);
    ans=0; % to store ongoing calculations
    n=length(observed); % to get length of given vector
    % doing cal. using for loop acc. to formula
    for i = 1:n
        % Access the current element using indexing
        n1=abs(simulated(i)-observed(i));
        n1=n1/observed(i);
        ans=ans+n1;
         % Perform operations on the current element
    end
    y=ans/n; % acc. to formula
    results=y; % to return final result

end



function results=MSRE(app, observed, simulated);
    ans=0; % to store ongoing cal.
    n=length(observed); % to get length of vector
    for i = 1:n
        % Access the current element using indexing
        n1=simulated(i)-observed(i);
        n1=n1/observed(i);
        n1=n1*n1;
        ans=ans+n1;
            % Perform operations on the current element
    end
    y=ans/n;
    results=y; % to return results
end



function results=MB(app, observed, simulated);
    sum_ob=sum(observed); % to get sum of observed vector
    sum_sim=sum(simulated); % to get sum of simulated vector
    y=sum_sim/sum_ob; % acc. to formula
    results=y; % to return final result
end


    
        
    
    
    
    
        
    
    
    
    
        
    
    
    
     
     
     

    
    



    



    
    





   




    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: UploaddataxlsxButton
        function UploaddataxlsxButtonPushed(app, event)
            [file, path] = uigetfile('*.xlsx', 'Select Excel File');
            
            if file ~= 0
                % Read the Excel file and store the data in the app's property
                % app.ExcelData = readtable(fullfile(path, file));
                app.data = xlsread(fullfile(path, file));
                %update lamp to ensure that your data is uploded
                app.Lamp.Color='g';
                % now read two column
              % x=app.TimeindexEditField.Value;
              % y=app.ObserveddataEditField.Value;
          % column1 = app.data(:, x); % Extract the first column (assuming it exists)
          % column3 = app.data(:, y);
                 % app.UITable.Data=readtable(fullfile(path, file));
                % Display the file name
                % app.FileLabel.Text = ['File: ', file];
                % 
                % % Update the table data
                % app.DataTable.Data = app.ExcelData;
            end

        end

        % Button pushed function: PlotButton
        function PlotButtonPushed(app, event)
               cla(app.UIAxes);
              
              % column for x axis
               x=app.TimeindexEditField.Value;
               %data from time index can be used as x axis
               x1 = app.data(:, x);
               % column for obs data 
              vy1=app.ObserveddataEditField.Value;
              obs_data=app.data(:, vy1);
                % column for sim data
              vy2=app.SimulateddataEditField.Value;
              sim_data=app.data(:, vy2);


          %var p to check which plot
          hold(app.UIAxes, 'on');
          p=app.SelectplottypeDropDown.Value;
          if strcmp(p, 'Scatter plot')
            % vy1=0;
            % vy2=0;
            cla(app.UIAxes);
            % c1 and c2 are column number of excel sheet for y1 and y2;
            % c1=app.y1EditField.Value;
            % c2=app.y2sEditField.Value;
            % c1=vy1;
            % c2=vy2;
            %here you can consider x1 as y1 and y1 as y2;
            % x1=app.data(:, c1);
            % y1=app.data(:, c2);
            %to plot scatter between x1 and y1 means y1 and y2;
            % scatter(app.UIAxes, obs_data, sim_data, 'red', 'filled' ,'DisplayName','');
           
           scatter(app.UIAxes, obs_data, sim_data, 'red', 'filled' );
           % plot(app.UIAxes, refline, 'k--');

             %below  lines are code to plot 45 degree line in scatter for
             %better understanding
             st=min(obs_data);
             st=max(st-10, 0);
             var1=max(sim_data);
             var2=max(obs_data);
             last=max(var2, var1);
             last=last+20;
             line_range = [st, last];
             plot(app.UIAxes, line_range, line_range, 'b--', 'LineWidth', 1.5);

             % refline(1,0);
             % axis equal;
           % axis(app.UIAxes, 'equal');
          xlabel(app.UIAxes, 'Observed data');
          ylabel(app.UIAxes, 'Simulated data');
          legend(app.UIAxes, 'off');

          % plot(app.UIAxes, line_range, line_range, 'b--', 'LineWidth', 2, 'DisplayName','');
          % Set the legend visibility to 'off'
          % set(app.UIAxes, 'LegendVisible', 'off');

          elseif strcmp(p,'Line plot')
              % cla(app.UIAxes,'reset');
              % hold(app.UIAxes, 'on');
              %first plot

              plot(app.UIAxes, x1, obs_data, 'b-', 'LineWidth', 2, 'DisplayName','Observed');
              plot(app.UIAxes, x1, sim_data, 'r--', 'LineWidth', 2 ,'DisplayName','Simulated');
                xlabel(app.UIAxes, 'Time');
               ylabel(app.UIAxes, 'Observed and Simulated data');
            % x1 = app.data(:, x); % Extract the first column (assuming it exists)
                % Set the legend visibility to 'on'
             % set(app.UIAxes, 'LegendVisible', 'on');

             legend(app.UIAxes, 'Location', 'best');
          else
              % cla(app.UIAxes,'reset');
              % hold(app.UIAxes, 'on');
              bar(app.UIAxes, x1, obs_data, 'b', 'DisplayName','Observed');
              bar(app.UIAxes, x1, sim_data, 'red','DisplayName','Simulated');
                         xlabel(app.UIAxes, 'Time');
               ylabel(app.UIAxes, 'Observed and Simulated data');
               % Set the legend visibility to 'on'
             % set(app.UIAxes, 'LegendVisible', 'on');

              legend(app.UIAxes, 'Location', 'best');
          end   
          
           % if vy1>0
           % y1 = app.data(:, vy1);  % extract data for y1 column vector to plot 
           %   %check it is line or bar
           % if strcmp(p,'Line plot')
           %  plot(app.UIAxes, x1, y1, 'b-', 'LineWidth', 2, 'DisplayName','Y1');
           %      else
           %     bar(app.UIAxes, x1, y1, 'b', 'DisplayName','Y1');
           % end  
           % end
           % 
           % if vy2>0
           %    y2=app.data(:, vy2); % extract data for y2
           %    % check it is line or bar
           %   if strcmp(p,'Line plot')
           %  plot(app.UIAxes, x1, y2, 'r--', 'LineWidth', 2 ,'DisplayName','Y2');
           %   else
           %     bar(app.UIAxes, x1, y2, 'red','DisplayName','Y2');
           %   end 
           % end
           % 
           % hold(app.UIAxes, 'off');
           % % check if we are plotting one of the bar of line graph then
           % % scatter input y1 and y2 should not be visible , set them off;
           % if (vy1>0 || vy2>0) 
           % app.y1EditField.Visible="off";
           % app.y2sEditField.Visible="off";
           % end

           % legend(app.UIAxes, 'Location', 'best');
          

        end

        % Value changed function: SelectplottypeDropDown
        function SelectplottypeDropDownValueChanged(app, event)
            value = app.SelectplottypeDropDown.Value;
            app.check="Scatter plot";
            
            % if strcmp(value, app.check) 
            % app.y1EditField.Visible="on";
            % app.y2sEditField.Visible="on";
            % else
            %   app.check="-1";
            % end

        end

        % Button pushed function: DownloadGraphButton
        function DownloadGraphButtonPushed(app, event)
         %centralized image
    % Open a file save dialog for the user to choose a location and filename
[filename, pathname] = uiputfile('*.tif', 'Save Figure As');

% Check if the user canceled the dialog
if isequal(filename, 0) || isequal(pathname, 0)
    disp('User canceled the operation.');
    return;
end

% Combine the chosen filename and pathname to get the full file path
filePath = fullfile(pathname, filename);

% Create a new figure and copy the entire UIAxes to the new figure
newFig = figure();

% Set the desired size for the new figure (larger than UIAxes)
newFigWidth = 1200;
newFigHeight = 800;
set(newFig, 'Position', [100, 100, newFigWidth, newFigHeight]);

% Calculate the position to center the copied UIAxes within the new figure
uiAxesPosition = app.UIAxes.Position;
uiAxesWidth = uiAxesPosition(3);
uiAxesHeight = uiAxesPosition(4);
newX = (newFigWidth - uiAxesWidth) / 2;
newY = (newFigHeight - uiAxesHeight) / 2;

% Copy the UIAxes contents to the new figure and set its position
ax = copyobj(app.UIAxes, newFig);
ax.Position = [newX, newY, uiAxesWidth, uiAxesHeight];

% Set the font size to ensure axis labels and tick labels are visible
ax.FontSize = 12; % Adjust the font size as needed

% Recreate the legend in the new figure based on DisplayName property
% lines = findall(ax, 'Type', 'line');
% legendEntries = {};
% for i = 1:numel(lines)
%     displayName = get(lines(i), 'DisplayName');
%     if ~isempty(displayName)
%         legendEntries = [legendEntries, displayName];
%     end
% end

% if ~isempty(legendEntries)
%     newLegend = legend(ax, legendEntries, 'Location', 'best');
% end

% % Save the new figure as an image (TIFF) to the chosen location
% drawnow; % Ensure the figure is fully updated before saving
% Create a new figure and copy the contents of the existing figure (newFig)
enlargedFig = figure('Visible', 'off'); % Create a new figure (not visible)
copyobj(allchild(newFig), enlargedFig); % Copy contents of newFig to enlargedFig

% Double the size of the enlargedFig
enlargedFig.Position = [100, 100, newFig.Position(3)*2, newFig.Position(4)*2];

% Center the figure on the screen
movegui(enlargedFig, 'center');

% Save the new figure as an image (TIFF) to the chosen location
drawnow; % Ensure the figure is fully updated before saving
print(enlargedFig, filePath, '-dtiff', '-r300');

% Close the new figure
close(enlargedFig);

% print(newFig, filePath, '-dtiff', '-r300');

% % Close the new figure
% close(newFig);
      
      
      

   


        end

        % Selection changed function: Tree
        function TreeSelectionChanged(app, event)
            app.sNodes = app.Tree.SelectedNodes;
        end

        % Callback function: Tree
        function TreeCheckedNodesChanged(app, event)
            checkedNodes = app.Tree.CheckedNodes;
            
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
          checkedNodes = app.Tree.CheckedNodes;

    % Initialize an empty cell array to store the text/names of checked nodes
    checkedNodeTexts = {};

    % Iterate through the checked nodes and extract the text/name of each node
         for i = 1:numel(checkedNodes)
        checkedNodeTexts{i} = checkedNodes(i).Text;
       end
     % Assuming you have an array of strings called checkedNodeTexts

    % data extraction for obs and simulated column vector form excel sheet
             sim=app.SimulatedEditField.Value;
             obs=app.ObservedEditField_2.Value;
            %here you can consider x1 as y1 and y1 as y2;
            simulated=app.data(:, sim);
            observed=app.data(:, obs);


resultVector = {};  % Initialize an empty cell array
   
for i = 1:length(checkedNodeTexts)
    if strcmp(checkedNodeTexts{i}, 'Mean Error (ME)')
         result = ME(app, observed, simulated);  % Call your function

        resultVector{end+1} = ['Mean Error (ME): ' num2str(result)];

     elseif strcmp(checkedNodeTexts{i}, 'Mean Absolute Error (MAE)')

        result = MAE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Mean Absolute Error (MAE): ' num2str(result)];

     elseif strcmp(checkedNodeTexts{i}, 'Absolute Maximum Error (AME)')

        result = AME(app, observed, simulated);  % Call your function
      
       resultVector{end+1} = ['Absolute Maximum Error (AME): ' num2str(result)];
   
 elseif strcmp(checkedNodeTexts{i}, 'Mean Squared Error (MSE)')

        result = MSE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Mean Squared Error (MSE): ' num2str(result)];

  elseif strcmp(checkedNodeTexts{i}, 'Root Mean Squared Error (RMSE)')

        result = RMSE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Root Mean Squared Error (RMSE): ' num2str(result)];

elseif strcmp(checkedNodeTexts{i}, 'Fourth Root Of The Mean Quadrupled Error (R4MS4E)')

        result = R4MS4E(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Fourth Root Of The Mean Quadrupled Error (R4MS4E): ' num2str(result)];

elseif strcmp(checkedNodeTexts{i}, 'Mean Squared Logarithmic Error (MSLE)')

        result = MSLE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Mean Squared Logarithmic Error (MSLE): ' num2str(result)];

elseif strcmp(checkedNodeTexts{i}, 'Mean Squared Derivative Error (MSDE)')

        result = MSDE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Mean Squared Derivative Error (MSDE): ' num2str(result)];

   elseif strcmp(checkedNodeTexts{i}, 'Fuzzy Mean Squared Error (FMSE)')

        result = FMSE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Fuzzy Mean Squared Error (FMSE): ' num2str(result)];

    elseif strcmp(checkedNodeTexts{i}, 'Peak Difference (PDIFF)')

        result = PDIFF(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Peak Difference (PDIFF): ' num2str(result)];

    elseif strcmp(checkedNodeTexts{i}, 'Fractional Standard Error (FSE)')

        result = FSE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Fractional Standard Error (FSE): ' num2str(result)];

     elseif strcmp(checkedNodeTexts{i}, 'Root Mean Squared Error to Standard Deviation Ratio (RSR)')

        result = RSR(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Root Mean Squared Error to Standard Deviation Ratio (RSR): ' num2str(result)];
   
elseif strcmp(checkedNodeTexts{i}, 'Inertia Root Mean Squared Error (IRMSE)')

        result = IRMSE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Inertia Root Mean Squared Error (IRMSE): ' num2str(result)];

elseif strcmp(checkedNodeTexts{i}, 'Mean Absolute Percentage Error (MAPE)')

        result = MAPE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Mean Absolute Percentage Error (MAPE): ' num2str(result)];

     elseif strcmp(checkedNodeTexts{i}, 'Median Absolute Percentage Error (MdAPE)')

        result = MdAPE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Median Absolute Percentage Error (MdAPE): ' num2str(result)];

      elseif strcmp(checkedNodeTexts{i}, 'Symmetric Mean Absolute Percentage Error (sMAPE)')

        result = sMAPE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Symmetric Mean Absolute Percentage Error (sMAPE): ' num2str(result)];

      elseif strcmp(checkedNodeTexts{i}, 'Percentage Error in Peak (PEP)')

        result = PEP(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Percentage Error in Peak (PEP): ' num2str(result)];

       elseif strcmp(checkedNodeTexts{i}, 'Runoff Coefficient Error (RR)')

        result = RR(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Runoff Coefficient Error (RR): ' num2str(result)];

     elseif strcmp(checkedNodeTexts{i}, 'Mean Relative Error (MRE)')
 
        result = MRE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Mean Relative Error (MRE): ' num2str(result)];

      elseif strcmp(checkedNodeTexts{i}, 'Relative Absolute Error (RAE)')

        result = RAE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Relative Absolute Error (RAE): ' num2str(result)];  

          elseif strcmp(checkedNodeTexts{i}, 'Relative Volume Error (RVE)')

        result = RVE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Relative Volume Error (RVE): ' num2str(result)];

          elseif strcmp(checkedNodeTexts{i}, 'Mean Absolute Relative Error (MARE)')

        result = MARE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Mean Absolute Relative Error (MARE): ' num2str(result)];

         elseif strcmp(checkedNodeTexts{i}, 'Mean Squared Relative Error (MSRE)')
       
        result = MSRE(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Mean Squared Relative Error (MSRE): ' num2str(result)];

     
         elseif strcmp(checkedNodeTexts{i}, 'Multiplicative Bias (MB)')

        result = MB(app, observed, simulated);  % Call your function
      
        resultVector{end+1} = ['Multiplicative Bias (MB): ' num2str(result)];   

    end
    
end

% combinedResultText = '';
% for i = 1:length(resultVector)
%     combinedResultText = [combinedResultText resultVector{i} '   '];
% end
output='';
for i=1: length(resultVector)
    output=sprintf('%s%s\n', output, resultVector{i});
end
app.ResultsTextArea.Value=output;

% Set the Value property of the EditField to display the combined text

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 742 626];
            app.UIFigure.Name = 'MATLAB App';

            % Create UploaddataxlsxButton
            app.UploaddataxlsxButton = uibutton(app.UIFigure, 'push');
            app.UploaddataxlsxButton.ButtonPushedFcn = createCallbackFcn(app, @UploaddataxlsxButtonPushed, true);
            app.UploaddataxlsxButton.Position = [42 580 115 23];
            app.UploaddataxlsxButton.Text = 'Upload data (.xlsx)';

            % Create SelectplottypeDropDownLabel
            app.SelectplottypeDropDownLabel = uilabel(app.UIFigure);
            app.SelectplottypeDropDownLabel.HorizontalAlignment = 'right';
            app.SelectplottypeDropDownLabel.Position = [43 392 87 22];
            app.SelectplottypeDropDownLabel.Text = 'Select plot type';

            % Create SelectplottypeDropDown
            app.SelectplottypeDropDown = uidropdown(app.UIFigure);
            app.SelectplottypeDropDown.Items = {'Line plot', 'Scatter plot', 'Bar plot'};
            app.SelectplottypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectplottypeDropDownValueChanged, true);
            app.SelectplottypeDropDown.Position = [145 392 100 22];
            app.SelectplottypeDropDown.Value = 'Line plot';

            % Create Lamp
            app.Lamp = uilamp(app.UIFigure);
            app.Lamp.Position = [177 581 20 20];
            app.Lamp.Color = [1 1 1];

            % Create Panel
            app.Panel = uipanel(app.UIFigure);
            app.Panel.TitlePosition = 'righttop';
            app.Panel.Position = [29 428 216 138];

            % Create ProvidecolumnindexLabel
            app.ProvidecolumnindexLabel = uilabel(app.Panel);
            app.ProvidecolumnindexLabel.BackgroundColor = [0.8 0.8 0.8];
            app.ProvidecolumnindexLabel.HorizontalAlignment = 'center';
            app.ProvidecolumnindexLabel.Position = [25 99 146 22];
            app.ProvidecolumnindexLabel.Text = 'Provide column index';

            % Create TimeindexEditFieldLabel
            app.TimeindexEditFieldLabel = uilabel(app.Panel);
            app.TimeindexEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeindexEditFieldLabel.Position = [25 70 63 22];
            app.TimeindexEditFieldLabel.Text = 'Time index';

            % Create TimeindexEditField
            app.TimeindexEditField = uieditfield(app.Panel, 'numeric');
            app.TimeindexEditField.Position = [119 70 49 22];

            % Create ObserveddataEditFieldLabel
            app.ObserveddataEditFieldLabel = uilabel(app.Panel);
            app.ObserveddataEditFieldLabel.HorizontalAlignment = 'right';
            app.ObserveddataEditFieldLabel.Position = [20 40 84 22];
            app.ObserveddataEditFieldLabel.Text = 'Observed data';

            % Create ObserveddataEditField
            app.ObserveddataEditField = uieditfield(app.Panel, 'numeric');
            app.ObserveddataEditField.Position = [119 40 49 23];

            % Create SimulateddataEditFieldLabel
            app.SimulateddataEditFieldLabel = uilabel(app.Panel);
            app.SimulateddataEditFieldLabel.HorizontalAlignment = 'right';
            app.SimulateddataEditFieldLabel.Position = [20 9 85 22];
            app.SimulateddataEditFieldLabel.Text = 'Simulated data';

            % Create SimulateddataEditField
            app.SimulateddataEditField = uieditfield(app.Panel, 'numeric');
            app.SimulateddataEditField.Position = [119 9 49 22];

            % Create Panel_2
            app.Panel_2 = uipanel(app.UIFigure);
            app.Panel_2.Position = [25 21 318 353];

            % Create UIAxes
            app.UIAxes = uiaxes(app.Panel_2);
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Box = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [24 90 276 250];

            % Create DownloadGraphButton
            app.DownloadGraphButton = uibutton(app.Panel_2, 'push');
            app.DownloadGraphButton.ButtonPushedFcn = createCallbackFcn(app, @DownloadGraphButtonPushed, true);
            app.DownloadGraphButton.Position = [179 53 105 23];
            app.DownloadGraphButton.Text = 'Download Graph';

            % Create PlotButton
            app.PlotButton = uibutton(app.Panel_2, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed, true);
            app.PlotButton.Position = [24 53 100 23];
            app.PlotButton.Text = 'Plot';

            % Create StatisticalPerformanceEvaluationMeasuresPanel
            app.StatisticalPerformanceEvaluationMeasuresPanel = uipanel(app.UIFigure);
            app.StatisticalPerformanceEvaluationMeasuresPanel.TitlePosition = 'centertop';
            app.StatisticalPerformanceEvaluationMeasuresPanel.Title = 'Statistical Performance Evaluation Measures';
            app.StatisticalPerformanceEvaluationMeasuresPanel.Position = [396 227 299 376];

            % Create Tree
            app.Tree = uitree(app.StatisticalPerformanceEvaluationMeasuresPanel, 'checkbox');
            app.Tree.SelectionChangedFcn = createCallbackFcn(app, @TreeSelectionChanged, true);
            app.Tree.Position = [18 52 264 249];

            % Create ScaledependentmeasuresNode
            app.ScaledependentmeasuresNode = uitreenode(app.Tree);
            app.ScaledependentmeasuresNode.Text = 'Scale-dependent measures';

            % Create MeanErrorMENode
            app.MeanErrorMENode = uitreenode(app.ScaledependentmeasuresNode);
            app.MeanErrorMENode.Text = 'Mean Error (ME)';

            % Create MeanAbsoluteErrorMAENode
            app.MeanAbsoluteErrorMAENode = uitreenode(app.ScaledependentmeasuresNode);
            app.MeanAbsoluteErrorMAENode.Text = 'Mean Absolute Error (MAE)';

            % Create AbsoluteMaximumErrorAMENode
            app.AbsoluteMaximumErrorAMENode = uitreenode(app.ScaledependentmeasuresNode);
            app.AbsoluteMaximumErrorAMENode.Text = 'Absolute Maximum Error (AME)';

            % Create MeanSquaredErrorMSENode
            app.MeanSquaredErrorMSENode = uitreenode(app.ScaledependentmeasuresNode);
            app.MeanSquaredErrorMSENode.Text = 'Mean Squared Error (MSE)';

            % Create RootMeanSquaredErrorRMSENode
            app.RootMeanSquaredErrorRMSENode = uitreenode(app.ScaledependentmeasuresNode);
            app.RootMeanSquaredErrorRMSENode.Text = 'Root Mean Squared Error (RMSE)';

            % Create FourthRootOfTheMeanQuadrupledErrorR4MS4ENode
            app.FourthRootOfTheMeanQuadrupledErrorR4MS4ENode = uitreenode(app.ScaledependentmeasuresNode);
            app.FourthRootOfTheMeanQuadrupledErrorR4MS4ENode.Text = 'Fourth Root Of The Mean Quadrupled Error (R4MS4E)';

            % Create MeanSquaredLogarithmicErrorMSLENode
            app.MeanSquaredLogarithmicErrorMSLENode = uitreenode(app.ScaledependentmeasuresNode);
            app.MeanSquaredLogarithmicErrorMSLENode.Text = 'Mean Squared Logarithmic Error (MSLE)';

            % Create MeanSquaredDerivativeErrorMSDENode
            app.MeanSquaredDerivativeErrorMSDENode = uitreenode(app.ScaledependentmeasuresNode);
            app.MeanSquaredDerivativeErrorMSDENode.Text = 'Mean Squared Derivative Error (MSDE)';

            % Create FuzzyMeanSquaredErrorFMSENode
            app.FuzzyMeanSquaredErrorFMSENode = uitreenode(app.ScaledependentmeasuresNode);
            app.FuzzyMeanSquaredErrorFMSENode.Text = 'Fuzzy Mean Squared Error (FMSE)';

            % Create PeakDifferencePDIFFNode
            app.PeakDifferencePDIFFNode = uitreenode(app.ScaledependentmeasuresNode);
            app.PeakDifferencePDIFFNode.Text = 'Peak Difference (PDIFF)';

            % Create ScaleindependentmeasuresNode
            app.ScaleindependentmeasuresNode = uitreenode(app.Tree);
            app.ScaleindependentmeasuresNode.Text = 'Scale-independent measures';

            % Create FractionalStandardErrorFSENode
            app.FractionalStandardErrorFSENode = uitreenode(app.ScaleindependentmeasuresNode);
            app.FractionalStandardErrorFSENode.Text = 'Fractional Standard Error (FSE)';

            % Create RootMeanSquaredErrortoStandardDeviationRatioRSRNode
            app.RootMeanSquaredErrortoStandardDeviationRatioRSRNode = uitreenode(app.ScaleindependentmeasuresNode);
            app.RootMeanSquaredErrortoStandardDeviationRatioRSRNode.Text = 'Root Mean Squared Error to Standard Deviation Ratio (RSR)';

            % Create InertiaRootMeanSquaredErrorIRMSENode
            app.InertiaRootMeanSquaredErrorIRMSENode = uitreenode(app.ScaleindependentmeasuresNode);
            app.InertiaRootMeanSquaredErrorIRMSENode.Text = 'Inertia Root Mean Squared Error (IRMSE)';

            % Create PercentagebasedmeasuresNode
            app.PercentagebasedmeasuresNode = uitreenode(app.Tree);
            app.PercentagebasedmeasuresNode.Text = 'Percentage-based measures';

            % Create MeanAbsolutePercentageErrorMAPENode
            app.MeanAbsolutePercentageErrorMAPENode = uitreenode(app.PercentagebasedmeasuresNode);
            app.MeanAbsolutePercentageErrorMAPENode.Text = 'Mean Absolute Percentage Error (MAPE)';

            % Create MedianAbsolutePercentageErrorMdAPENode
            app.MedianAbsolutePercentageErrorMdAPENode = uitreenode(app.PercentagebasedmeasuresNode);
            app.MedianAbsolutePercentageErrorMdAPENode.Text = 'Median Absolute Percentage Error (MdAPE)';

            % Create SymmetricMeanAbsolutePercentageErrorsMAPENode
            app.SymmetricMeanAbsolutePercentageErrorsMAPENode = uitreenode(app.PercentagebasedmeasuresNode);
            app.SymmetricMeanAbsolutePercentageErrorsMAPENode.Text = 'Symmetric Mean Absolute Percentage Error (sMAPE)';

            % Create PercentageErrorinPeakPEPNode
            app.PercentageErrorinPeakPEPNode = uitreenode(app.PercentagebasedmeasuresNode);
            app.PercentageErrorinPeakPEPNode.Text = 'Percentage Error in Peak (PEP)';

            % Create RunoffCoefficientErrorRRNode
            app.RunoffCoefficientErrorRRNode = uitreenode(app.PercentagebasedmeasuresNode);
            app.RunoffCoefficientErrorRRNode.Text = 'Runoff Coefficient Error (RR)';

            % Create RelativemeasuresNode
            app.RelativemeasuresNode = uitreenode(app.Tree);
            app.RelativemeasuresNode.Text = 'Relative measures';

            % Create MeanRelativeErrorMRENode
            app.MeanRelativeErrorMRENode = uitreenode(app.RelativemeasuresNode);
            app.MeanRelativeErrorMRENode.Text = 'Mean Relative Error (MRE)';

            % Create RelativeAbsoluteErrorRAENode
            app.RelativeAbsoluteErrorRAENode = uitreenode(app.RelativemeasuresNode);
            app.RelativeAbsoluteErrorRAENode.Text = 'Relative Absolute Error (RAE)';

            % Create RelativeVolumeErrorRVENode
            app.RelativeVolumeErrorRVENode = uitreenode(app.RelativemeasuresNode);
            app.RelativeVolumeErrorRVENode.Text = 'Relative Volume Error (RVE)';

            % Create MeanAbsoluteRelativeErrorMARENode
            app.MeanAbsoluteRelativeErrorMARENode = uitreenode(app.RelativemeasuresNode);
            app.MeanAbsoluteRelativeErrorMARENode.Text = 'Mean Absolute Relative Error (MARE)';

            % Create MeanSquaredRelativeErrorMSRENode
            app.MeanSquaredRelativeErrorMSRENode = uitreenode(app.RelativemeasuresNode);
            app.MeanSquaredRelativeErrorMSRENode.Text = 'Mean Squared Relative Error (MSRE)';

            % Create MultiplicativeBiasMBNode
            app.MultiplicativeBiasMBNode = uitreenode(app.RelativemeasuresNode);
            app.MultiplicativeBiasMBNode.Text = 'Multiplicative Bias (MB)';

            % Assign Checked Nodes
            app.Tree.CheckedNodesChangedFcn = createCallbackFcn(app, @TreeCheckedNodesChanged, true);

            % Create RunButton
            app.RunButton = uibutton(app.StatisticalPerformanceEvaluationMeasuresPanel, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Position = [98 9 100 23];
            app.RunButton.Text = 'Run';

            % Create ObservedEditField_2Label
            app.ObservedEditField_2Label = uilabel(app.StatisticalPerformanceEvaluationMeasuresPanel);
            app.ObservedEditField_2Label.HorizontalAlignment = 'right';
            app.ObservedEditField_2Label.Position = [20 313 60 22];
            app.ObservedEditField_2Label.Text = 'Observed ';

            % Create ObservedEditField_2
            app.ObservedEditField_2 = uieditfield(app.StatisticalPerformanceEvaluationMeasuresPanel, 'numeric');
            app.ObservedEditField_2.Position = [90 313 46 22];

            % Create SimulatedEditFieldLabel
            app.SimulatedEditFieldLabel = uilabel(app.StatisticalPerformanceEvaluationMeasuresPanel);
            app.SimulatedEditFieldLabel.HorizontalAlignment = 'right';
            app.SimulatedEditFieldLabel.Position = [152 313 58 22];
            app.SimulatedEditFieldLabel.Text = 'Simulated';

            % Create SimulatedEditField
            app.SimulatedEditField = uieditfield(app.StatisticalPerformanceEvaluationMeasuresPanel, 'numeric');
            app.SimulatedEditField.Position = [219 313 46 22];

            % Create ResultsdisplayPanel
            app.ResultsdisplayPanel = uipanel(app.UIFigure);
            app.ResultsdisplayPanel.Title = 'Results display';
            app.ResultsdisplayPanel.Position = [416 32 260 172];

            % Create ResultsTextArea
            app.ResultsTextArea = uitextarea(app.ResultsdisplayPanel);
            app.ResultsTextArea.Position = [11 8 234 130];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = btpapp1_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
%% --------------------------------------------------------------------

%%  File     : wumpus.pl
%%  Author   : Shajid Mohammad (shajidm@student.unimelb.edu.au)
%%  Origin   : Thu May 24, 2018
%%  Purpose  : Assignment4 Project

%% --------------------------------------------------------------------

%% INTRODUCTION
%% --------------
%% Wumpus is a planning problem. 
%% You need to find and kill a Wumpus hiding in an unknown maze. 
%% The player sends in a series of disposable robots each with a 
%% fixed list of instructions to follow. 
%% Each robot follows its instructions until it is destroyed 
%% or it finishes the instructions, or it kills the Wumpus.
%% After a robot is finished the player gets feedback on what 
%% the robot sensed following the instructions. 
%% The aim is to explore the maze and find and shoot the wumpus.

%% --------------------------------------------------------------------

%% Module declaration
%% A Prolog module is a collection of predicates which defines a 
%% public interface by means of a set of provided predicates and
%% operators, here we declare initialState, guess and updateState
%% functions, so it can be called outside the scope.

:- module(wumpus,[initialState/5, guess/3, updateState/4]).

%% --------------------------------------------------------------------

%% initialState(+NR, +NC, +XS, +YS, -State0)
%% Takes as input the number of rows NR and number of columns NC 
%% in the game, and the starting position (XS, Y S) and outputs  
%% State0 an initial state representation for the game.
%% State0 stores all the VisitedCells and all the available vantagepoints
%% to soot the wumpus and also the danger points which are required to be
%% avoided by the robot, also the CurrentGameState is used to get the
%% Wumpus location, initial position and also the boundary points.

initialState(NR, NC, XS, YS, State0):-
    generateNewPairs(NR,NC,KeyValuePairs), 
    newfactsGenerator(KeyValuePairs,NR,NC),
    VisitedCells = [(XS,YS)], WumpusLocation = empty, 
    CurrentGameState = [(NR,NC),(XS,YS),WumpusLocation],
    DangerZones = [(XS,YS)], VantagePoints = [], 
    State0 = (VisitedCells,CurrentGameState,VantagePoints,DangerZones).

%% --------------------------------------------------------------------

%% generateNewPairs is used to create NoOfRowsInMap X NoOfColsInMap pairs
%% succ(Low, NextRowNo) returns True if NextRowNo = Low + 1 and Low >= 0. 

generateNewPairs(NoOfRowsInMap, NoOfColsInMap, KeyValuePair) :-
    NextRowNo is NoOfRowsInMap+1,succ(Low, NextRowNo),
    NextColNo is NoOfColsInMap+1,succ(High, NextColNo),
    bagof(Temppair, generatePairs(Low, High, Temppair), KeyValuePair).

%% --------------------------------------------------------------------

%% pair is used to create pairs with values from a range of low to high
%% between(Low, High, Pair), returns Pair which is successively bound 
%% to all integers between Low and High.

generatePairs(Low, High, (X,Y)) :-
between(1, Low, X),
between(1, High, Y).

%% --------------------------------------------------------------------

%% newfactsGenerator is used to Create new facts or predicates from known
%% calculated relations and the Cut operator (!) in Prolog is used
%% to stop backtracking, to increase performance.

newfactsGenerator([],_,_).
newfactsGenerator([(X,Y)|Points],NoofRows,NoofCol):-
    WestofX is X - 1,
    EastofX is X + 1,
    SouthOfY is Y + 1,
    NorthOfY is Y - 1,
    (  WestofX > 0 -> rowsCoordinatesMapper((X,Y),(WestofX,Y))
    ;! 
    ),
    (  EastofX =< NoofRows -> rowsCoordinatesMapper((X,Y),(EastofX,Y))
    ;! 
    ),
    (  NorthOfY > 0 -> columnsCoordinatesMapper((X,Y),(X,NorthOfY))
    ;! 
    ),
    (  SouthOfY =< NoofCol -> columnsCoordinatesMapper((X,Y),(X,SouthOfY))
    ;!
    ),
    newfactsGenerator(Points,NoofRows,NoofCol).

%% --------------------------------------------------------------------

%% guess(+State0, -State, -Guess) 
%% Given the current state State0 returns a new state 
%% and a Guess which is a list of north, east, south, west, shoot 
%% which are instructions for the robot. 
%% The total energy of the Guess must be at most 100 energy
%% each move costs 1 energy, and each shoot requires 5 energy.
%% Here we take State0 as the input and result a response of a new 
%% Guess and an update new state. If the WumpusLocation is equal to
%% empty then we update the state and give a new guess as a plan else
%% we append a shoot instruction to the KillPlan and give a Guess.

guess(State0, State, Guess):-
    (VisitedCells,CurrentGameState,VantagePoints,DangerZones) = State0,
    [BoundaryPoints,InitialLocation,WumpusLocation] = CurrentGameState, 
    (   WumpusLocation == empty ->
    	BoundaryPoints = (NoofRows,NoofCols),
        generateNewPairs(NoofRows,NoofCols,KeyValuePairs),
        subtract(KeyValuePairs,VisitedCells,UnVisitedCells),
        detectPath(InitialLocation,UnVisitedCells,DangerZones,Guess),
        State = State0 
    ;
        bestVantagePathToWumpus(InitialLocation,VantagePoints,
            WumpusLocation,DangerZones,KillPlan),
        append(KillPlan,[shoot],Guess),
        State = State0
    ).

%% --------------------------------------------------------------------

%% bestVantagePathToWumpus is used to get the vantage Path from a set of 
%% vantagepoints to kill the wumpus.

bestVantagePathToWumpus(InitialLocation,VantagePoints,(WX,WY),
	DangerZones,KillPlan):-
	detectPath(InitialLocation,VantagePoints,DangerZones,KillPlan),
	find(InitialLocation,(SX,SY),KillPlan),
	validateKillPlan(SX,WX,SY,WY,KillPlan).

%% --------------------------------------------------------------------

%% detectPath is used to detect untraversed paths by utilziing
%% the UnVisitedCells list.

detectPath(InitialLocation,UnVisitedCells,DangerZones,Guess):-
    member(Cell,UnVisitedCells),
    find(InitialLocation,Cell,DangerZones,Guess).

%% --------------------------------------------------------------------

%%  validateKillPlan is used to validate the Kill Plan.

validateKillPlan(SX,WX,SY,WY,P):-    
    last(P,Move),
    (   SX =:= WX -> (SY > WY -> Move = north; Move = south);
        SY =:= WY -> (SX < WX -> Move = east; Move = west)).

%% --------------------------------------------------------------------

%% find finds a simple Path from Start to End
%% If already visited Previous then we 
%% dont visit that places once again.

find(Start, End, Path) :-
	find(Start, End, [Start], Path).

find(Start, Start, _, []).
find(Start, End, Previous, [Direction|Path]) :-
    edge(Start, Direction, Med),
    \+ member(Med, Previous), 
    find(Med, End, [Med|Previous], Path).

%% --------------------------------------------------------------------

%% updateState(+State0, +Guess, +Feedback, -State)
%% It takes as input the state State0, 
%% the previous guess Guess and the feedback from the
%% guess Feedback and returns a new updated state State.
%% If the WumpusLocation is empty, it checks whether wumpus is returned the 
%% feedback and if it is, it appends the wumpus position in the DangerZones
%% list and if a pit is returned in the feedback it appends this as well to
%% the DangerZones list qnd visited cells list and if wumpus location is 
%% not empty then it maps the poistions and shoots to kill the wumpus.

updateState(State0, Guess, Feedback, State):-
    (VisitedCells,CurrentGameState,VantagePoints,DangerZones) = State0,
    [BoundaryPoints,InitialLocation,WumpusLocation] = CurrentGameState,
    (       WumpusLocation == empty ->
            (
            member(wumpus,Feedback) ->
            length(Feedback,WP),
            getFirstNElements(WP,Guess,PathToWumpus),
            find(InitialLocation,WumpusPosition,PathToWumpus), 
            NewGameState = [BoundaryPoints,InitialLocation,WumpusPosition],
            getVantagePoints(BoundaryPoints,WumpusPosition,ShootPositions), 
            append([WumpusPosition],DangerZones,NewDangerZones),
            State = (VisitedCells,NewGameState,ShootPositions,NewDangerZones)
    		;
            (member(pit,Feedback)) ->
                length(Feedback,Pit),
                getFirstNElements(Pit,Guess,PathToPit),                
                find(InitialLocation,PitPosition,PathToPit),
                append([PitPosition],DangerZones,NewDangerZones),
                savePath(InitialLocation,PathToPit,UpdatedVisited,[]),
                append(UpdatedVisited,VisitedCells,NewVisited1),
                sort(NewVisited1,NewVisited),
                State = (NewVisited,CurrentGameState,ShootPositions,
                    NewDangerZones);
                savePath(InitialLocation,Guess,UpdatedVisited,[]),
                append(UpdatedVisited,VisitedCells,NewVisited),
                State = (NewVisited,CurrentGameState,ShootPositions,
                    DangerZones)
            )
    ;
            VantagePoints = [_|RestShootPos],
            State = (VisitedCells,CurrentGameState,RestShootPos,DangerZones)
    ).

%% --------------------------------------------------------------------


%% savePath is used to record untraversed paths.

savePath(_,[],VisitedCells,VisitedCells).
savePath(InitialLocation,[Loc|NewGuess],VisitedCells,Path):-
    find(InitialLocation,Last,[Loc]),
    append([Last],Path,TempPath),
    savePath(Last,NewGuess,VisitedCells,TempPath).

%% --------------------------------------------------------------------

%% getVantagePoints is used to get all the available vantagepoints
%% to shoot the wumpus.

getVantagePoints((NR,NC),(X,Y),VantagePoints):-
    vantagePointsForRows(NR,NC,(X,Y),XCoords,[]),
    vantagePointsForColumns(NR,NC,(X,Y),YCoords,[]),
    append(XCoords,YCoords,Coords), 
    OldPairs = [(X,Y),(X,1),(X,NC),(NR,Y),(1,Y)],
    subtract(Coords,OldPairs,VantagePoints).

%% --------------------------------------------------------------------

%% vantagePointsForRows is used as an iterator to get 
%% all the available vantagepoints from East to West.

vantagePointsForRows(NR,NC,(X,Y),VantagePoints,Point):-
    (   
    NR =:= 0 -> VantagePoints = Point
   	;
    NR > 0 -> append([(NR,Y)],Point,TempPoint),
              NNR is NR - 1,
              vantagePointsForRows(NNR,NC,(X,Y),VantagePoints,TempPoint)
    ).

%% --------------------------------------------------------------------

%% vantagePointsForColumns is used as an iterator to get 
%% all the available vantagepoints from North to South.

vantagePointsForColumns(NR,NC,(X,Y),VantagePoints,Point):-
    (   
    NC =:= 0 ->VantagePoints = Point
    ;
    NC > 0 -> append([(X,NC)],Point,TempPoint),
              NNC is NC - 1,
              vantagePointsForColumns(NR,NNC,(X,Y),VantagePoints,TempPoint)
    ).

%% --------------------------------------------------------------------

%% rowsCoordinatesMapper is used to assert pairs with respect to the
%% X direction i.e. from East to West.

rowsCoordinatesMapper(Pair1,Pair2):-
    Pair1 = (X,Y),
    Pair2 = (Pair2X,Y),
    (
    Pair2X - X > 0 ->
    assert(edge(Pair1,east,Pair2))
    ;
    assert(edge(Pair1,west,Pair2))
    ).

%% --------------------------------------------------------------------

%% columnsCoordinatesMapper is used to assert pairs with respect 
%% to the Y direction i.e. from North to South.

columnsCoordinatesMapper(Pair1,Pair2):-
    Pair1 = (X,Y),
    Pair2 = (X,Pair2Y),
    (
    Pair2Y - Y > 0 ->
    assert(edge(Pair1,south,Pair2))
    ;
    assert(edge(Pair1,north,Pair2))
    ).

%% --------------------------------------------------------------------

%% getFirstNElements is used to get the first N elements from a list.

getFirstNElements(Number,List,Startingelements) :- 
    findall(Temp, (nth1(I,List,Temp), I =< Number), Startingelements).
%% --------------------------------------------------------------------
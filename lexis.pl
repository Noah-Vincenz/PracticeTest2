%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Question 1 (50%)                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (a) eval(+E, +Env, -V)

% Addition case
eval(E1 + E2, Env, V) :-
  eval(E1, Env, V1),
  eval(E2, Env, V2),
  V is V1 + V2.

% Multiplication case
eval(E1 * E2, Env, V) :-
  eval(E1, Env, V1),
  eval(E2, Env, V2),
  V is V1 * V2.

% Minus case
eval(- E, Env, V) :-
  eval(E, Env, V1),
  V is -V1.

% sin case
eval(sin(E), Env, V) :-
  eval(E, Env, V1),
  V is sin(V1).

% cos case
eval(cos(E), Env, V) :-
  eval(E, Env, V1),
  V is cos(V1).

% Number case
eval(C, Env, C) :-
  float(C),!.

% Variable case
eval(Var, [(Var,X)|Tail], X).

eval(Var, [X|Tail], V) :-
  eval(Var, Tail, V).


% (b) commutes(+E1, +E2)
commutes(E,E).

commutes(sin(E1), sin(E2)) :-
  commutes(E1,E2).

commutes(cos(E1), cos(E2)) :-
  commutes(E1,E2).

commutes(E1+E2,E_1+E_2) :-
  commutes(E1,E_1),
  commutes(E2,E_2).

commutes(E1+E2,E_1+E_2) :-
  commutes(E1,E_2),
  commutes(E2,E_1).

commutes(E1*E2,E_1*E_2) :-
  commutes(E1,E_1),
  commutes(E2,E_2).

commutes(E1*E2,E_1*E_2) :-
  commutes(E1,E_2),
  commutes(E2,E_1).


% (c) diff(+E, +V, -D)
:- ensure_loaded(support).
% case 4
diff(E1 + E2, V, D) :-
  !,
  diff(E1, V, D1),
  diff(E2, V, D2),
  simplify(D1 + D2, D).

% case 5
diff(E1 * E2, V, D) :-
  !,
  diff(E1, V, D1),
  diff(E2, V, D2),
  simplify(D1, D_1),
  simplify(D2, D_2),
  simplify(E1 * D_2, D3),
  simplify(E2 * D_1, D4),
  simplify(D3 + D4, D).

% case 6
diff(sin(E), V, D) :-
  !,
  diff(E, V, D1),
  simplify(cos(E) * D1, D).

% case 7
diff(cos(E), V, D) :-
  !,
  diff(E, V, D1),
  simplify(-sin(E) * D1, D).

% case 2
diff(E, V, 1.0) :-
  E == V,!.

% special case
diff(-V, V, -1.0) :- !.

% case 3
diff(E, V, 0.0).



% (d) maclaurin(+E, +X, +N, -V)
% (d) maclaurin(+E, +X, +N, -V)
factorial(N,Result):-
  factorial(N,1,Result).
factorial(0,Acc,Acc):-!.
factorial(N,Acc,Result):-
  N > 0,
  NewAcc is Acc*N,
  N1 is N-1,
  factorial(N1,NewAcc,Result).

exponent(Base,Ex,Result):-
  exponent(Base,Ex,1,Result).
exponent(_,0,Acc,Acc):-!.
exponent(Base,Exponent,Acc,Result):-
  Exponent1 is Exponent-1,
  Acc1 is Acc*Base,
  exponent(Base,Exponent1,Acc1,Result).

kDiff(1,E,X,D):-
  !,
  diff(E,X,Tmp),
  D = Tmp.
kDiff(K,E,X,D):-
  diff(E,X,Tmp),
  K1 is K-1,
  kDiff(K1,Tmp,X,D).

maclaurin(E,X,N,V):-
  N1 is N-1,
  maclaurin(E,X,N1,0.0,V),!.

maclaurin(E,_,0,Acc,V):-
  !,eval(E,[(x,0.0)],F0),
  V is Acc+F0.

maclaurin(E,X,N,Acc,Result):-
  factorial(N,FK),
  kDiff(N,E,x,Tmp0),
  Tmp1 = Tmp0/FK,
  exponent(X,N,RX),
  Tmp2 = Tmp1*RX,
  eval(Tmp2,[(x,0.0)],Tmp3),
  NewAcc is Acc+ Tmp3,
  N1 is N-1,
  simplify(NewAcc,Simplified),
  maclaurin(E,X,N1,Simplified,Result).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Question 2 (50%)                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% (a) b_to_left(+Pos, -H).

b_to_left(Pos, H) :-
  comp_sum(Pos, 0, 0, H).

comp_sum([b|Tail], NoOfBlackTiles, Acc, H) :-
  NewNoOfBlackTiles is NoOfBlackTiles + 1,
  comp_sum(Tail, NewNoOfBlackTiles, Acc, H).

comp_sum([w|Tail], NoOfBlackTiles, Acc, H) :-
  NewAcc is Acc + NoOfBlackTiles,
  comp_sum(Tail, NoOfBlackTiles, NewAcc, H).

comp_sum([e|Tail], NoOfBlackTiles, Acc, H) :-
  comp_sum(Tail, NoOfBlackTiles, Acc, H).

comp_sum([], NoOfBlackTiles, Acc, Acc).






% (b) move(+Pos, -NewPos)

% find solutions with 3B, 3W, 1E
% find index of E in Pos
% find index of E in NewPos
% get jump by taking away new index from old (abs)
% check if jump != 0 and jump <= 3
% check if 5 elements have stayed in same place
move(Pos,NewPos):-
  move(e,Pos,NewPos).

move(X,[H1,H2,H3,X|T],[X,H2,H3,H1|T]).
move(X,[H1,H2,X|T],[X,H2,H1|T]).
move(X,[H1,X|T],[X,H1|T]).
move(X,[X,H1,H2,H3|T],[H3,H1,H2,X|T]).
move(X,[X,H1,H2|T],[H2,H1,X|T]).
move(X,[X,H|T],[H,X|T]).
move(X,[H1,H2|T],[H1|Result]):-
    move(X,[H2|T],Result).


/*
move(Pos, NewPos) :-
  all_solutions(3, 3, 1, [], NewPos),
  five_remained_the_same(Pos, NewPos, 0),
  find_blank(Pos, 1, PosIndex),
  find_blank(NewPos, 1, NewPosIndex),
  Jump is abs(NewPosIndex - PosIndex),
  Jump \== 0,
  Jump =< 3.


all_solutions(0, 0, 0, Acc, Acc).

all_solutions(BCount, WCount, ECount, Acc, Solution) :-
  BCount > 0,
  NewBCount is BCount - 1,
  all_solutions(NewBCount, WCount, ECount, [b|Acc], Solution).

all_solutions(BCount, WCount, ECount, Acc, Solution) :-
  WCount > 0,
  NewWCount is WCount - 1,
  all_solutions(BCount, NewWCount, ECount, [w|Acc], Solution).

all_solutions(BCount, WCount, ECount, Acc, Solution) :-
  ECount > 0,
  all_solutions(BCount, WCount, 0, [e|Acc], Solution).


five_remained_the_same([], [], 5).

five_remained_the_same([X|PosTail], [Y|NewPosTail], RemainedTheSame) :-
  X == Y,!,
  NewRemainedTheSame is RemainedTheSame + 1,
  five_remained_the_same(PosTail, NewPosTail, NewRemainedTheSame).

five_remained_the_same([X|PosTail], [Y|NewPosTail], RemainedTheSame) :-
  five_remained_the_same(PosTail, NewPosTail, RemainedTheSame).


find_blank([e|Tail], Acc, Acc) :- !.
find_blank([_|Tail], Acc, BlankPos) :-
  NewAcc is Acc + 1,
  find_blank(Tail, NewAcc, BlankPos).
*/

reverse(List, NewList) :-
  reverse_list(List, [], NewList).

reverse_list([], NewList, NewList).

reverse_list([X|Tail], Acc, NewList) :-
  reverse_list(Tail, [X|Acc], NewList).



% (c) search_agenda(+Agenda, -Visited, -Final)

%search_agenda(Agenda, Visited, Final) :-

% go through all nodes in agenda until found a final node reachable from one of the nodes


search_agenda(Agenda,Visited,Final):-
  search_agenda(Agenda,-1,Visited,Final).

search_agenda([Node|_], ID, [], Node):-
  node_id(Node,ID),
  goal(Node).


search_agenda([Node|T], C, [Node|Visited], Final):-
    \+goal(Node),
    C1 is C+1,
    node_pos(Node,Pos),
    findall(N,move(Pos,N),Children),
    node_g(Node,G),
    node_id(Node,C1),
    G1 is G+1,
    makeListOfNodes(G1,C1,Children,ChildrenNodes),
    merge(T,ChildrenNodes,Merged),
    sort(Merged,Sorted),
    search_agenda(Sorted,C1,Visited,Final).

merge([],List2,List2).
merge([H|T],List2,[H|NewList]):-
  merge(T,List2,NewList).

makeListOfNodes(_,_,[],[]).
makeListOfNodes(G,Pid,[Pos|T],[Node|Nodes]):-
  b_to_left(Pos,H),
  make_node(G,H,Pos,Pid,Node),
  makeListOfNodes(G,Pid,T,Nodes).

goal(Node):-
  node_pos(Node,Pos),
  b_to_left(Pos,0).

% (d) trace_moves(+Final, +Visited, -Seq)
%trace_moves(n(10, 1, [b,b,e,b,w,w,w], _, 0),
%         [n(9, 0, [b,b,b,e,w,w,w], 0, none)], Seq).

trace_moves(Final, [], [Pos-H]):-
  node_pos(Final,Pos),
  b_to_left(Pos,H).

trace_moves(Final, [Node|T], [Pos-H|Seq]):-
  node_pos(Node,Pos),
  b_to_left(Pos,H),
  trace_moves(Final,T,Seq).



% (d) trace_moves(+Final, +Visited, -Seq)
